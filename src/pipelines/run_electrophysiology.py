import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import signal
from specparam import SpectralModel

# Configurações estéticas de publicação científica
plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
plt.rcParams['font.sans-serif'] = 'DejaVu Sans'
plt.rcParams['figure.titlesize'] = 14
plt.rcParams['axes.titlesize'] = 11

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
DATA_DIR = os.path.join(REPO_ROOT, "data", "raw") if os.path.exists(os.path.join(REPO_ROOT, "data", "raw")) else os.path.join(REPO_ROOT, "X")
OUTPUT_DIR = os.path.join(REPO_ROOT, "results")
os.makedirs(os.path.join(OUTPUT_DIR, "figures"), exist_ok=True)
os.makedirs(os.path.join(OUTPUT_DIR, "tables"), exist_ok=True)

FS = 200 # Taxa de amostragem da OpenBCI Ganglion (Hz)
CHANNELS = ['F3', 'F4', 'P3/O1', 'P4/O2'] # Mapeamento 10-20 do Ultracortex Mark IV

# Bandas canônicas de EEG
BANDS = {
    'Delta': (1.0, 4.0),
    'Theta': (4.0, 8.0),
    'Alpha': (8.0, 13.0),
    'Beta':  (13.0, 30.0),
    'Gamma': (30.0, 45.0)
}

def load_and_clean_eeg(filepath, apply_alias_notch=True):
    """
    Pipeline rigoroso de pré-processamento de fase zero:
    1. Detrend linear.
    2. Filtro Notch 60 Hz (Q=30) - Fundamental da rede elétrica brasileira.
    3. Filtro Notch 20 Hz (Q=30) - Eliminação do 3º harmônico (180 Hz) rebatido por aliasing (180 - 200 = 20 Hz).
    4. Filtro Notch 80 Hz (Q=30) - Eliminação do 2º harmônico (120 Hz) rebatido por aliasing (120 - 200 = 80 Hz).
    5. Filtro Passa-Faixa Butterworth 1-45 Hz via Seções de Segunda Ordem (SOS) de fase zero.
    """
    df = pd.read_csv(filepath, skiprows=1, header=None, comment='%', usecols=[1, 2, 3, 4])
    raw = df.to_numpy(dtype=float)
    
    # 1. Detrend
    d = signal.detrend(raw, axis=0)
    
    # 2. Notch 60 Hz
    b60, a60 = signal.iirnotch(60.0, 30.0, fs=FS)
    d = signal.filtfilt(b60, a60, d, axis=0)
    
    # 3 e 4. Notches de Aliasing (180 Hz -> 20 Hz e 120 Hz -> 80 Hz)
    if apply_alias_notch:
        b20, a20 = signal.iirnotch(20.0, 30.0, fs=FS)
        d = signal.filtfilt(b20, a20, d, axis=0)
        b80, a80 = signal.iirnotch(80.0, 30.0, fs=FS)
        d = signal.filtfilt(b80, a80, d, axis=0)
        
    # 5. Bandpass 1-45 Hz (SOS ordem 4)
    sos = signal.butter(4, [1.0, 45.0], btype='bandpass', fs=FS, output='sos')
    clean = signal.sosfiltfilt(sos, d, axis=0)
    return clean

def epoch_data_with_artifact_rejection(data, epoch_sec=2.0, threshold_uV=500.0):
    """
    Divide o sinal contínuo em épocas de 2.0 segundos e aplica rejeição automática de artefatos:
    Descarta épocas com amplitude pico-a-pico excessiva nos canais frontais (> threshold_uV, ex.: movimentos bruscos, piscadas extremas)
    ou épocas com linha plana/desconexão (< 2.0 uV).
    """
    ep_len = int(epoch_sec * FS)
    n_total_epochs = data.shape[0] // ep_len
    
    valid_epochs = []
    rejected_count = 0
    
    for i in range(n_total_epochs):
        ep = data[i*ep_len : (i+1)*ep_len, :]
        # Foco nos canais frontais F3 e F4 (índices 0 e 1) para FAA e conectividade
        p2p_frontal = np.ptp(ep[:, :2], axis=0)
        
        if np.any(p2p_frontal > threshold_uV) or np.any(p2p_frontal < 2.0):
            rejected_count += 1
        else:
            valid_epochs.append(ep)
            
    valid_data = np.array(valid_epochs) # Formato: (n_epochs, n_samples_ep, n_channels)
    return valid_data, rejected_count, n_total_epochs

def compute_epoch_psd_and_relative_faa(epochs):
    """
    Calcula potência relativa em cada época e deriva a Assimetria Alfa Frontal (FAA) relativa:
    FAA_rel = ln(Rel_Alpha_F4) - ln(Rel_Alpha_F3)
    Gera distribuição de FAA para testes estatísticos e intervalos de confiança (Bootstrap).
    """
    n_ep = epochs.shape[0]
    faa_epochs = []
    f3_rel_alpha = []
    f4_rel_alpha = []
    
    for i in range(n_ep):
        ep = epochs[i, :, :]
        f, psd = signal.welch(ep, fs=FS, nperseg=epochs.shape[1], axis=0)
        
        mask_tot = (f >= 1.0) & (f <= 45.0)
        mask_alpha = (f >= 8.0) & (f <= 13.0)
        
        trapz = getattr(np, 'trapezoid', getattr(np, 'trapz', None))
        tot_pwr = trapz(psd[mask_tot, :], f[mask_tot], axis=0)
        alpha_pwr = trapz(psd[mask_alpha, :], f[mask_alpha], axis=0)
        
        rel_alpha = alpha_pwr / (tot_pwr + 1e-12)
        
        # F3 = canal 0, F4 = canal 1
        faa = np.log(rel_alpha[1] + 1e-12) - np.log(rel_alpha[0] + 1e-12)
        faa_epochs.append(faa)
        f3_rel_alpha.append(rel_alpha[0] * 100.0)
        f4_rel_alpha.append(rel_alpha[1] * 100.0)
        
    return np.array(faa_epochs), np.array(f3_rel_alpha), np.array(f4_rel_alpha)

def bootstrap_ci(data, n_boot=2000, ci=95):
    """Calcula intervalo de confiança não paramétrico via Bootstrap (percentil)."""
    boot_means = np.empty(n_boot)
    n = len(data)
    for i in range(n_boot):
        sample = np.random.choice(data, size=n, replace=True)
        boot_means[i] = np.mean(sample)
    lower = np.percentile(boot_means, (100 - ci) / 2.0)
    upper = np.percentile(boot_means, 100 - (100 - ci) / 2.0)
    return np.mean(data), lower, upper, np.std(data)

def compute_wpli_with_surrogates(epochs, band=(8.0, 13.0), n_surrogates=500):
    """
    Calcula o wPLI (Weighted Phase Lag Index) inter-épocas entre F3 (CH1) e F4 (CH2)
    e testa a significância estatística contra uma distribuição nula de fase aleatória (Surrogate Testing).
    """
    n_ep, ep_len, n_ch = epochs.shape
    sos = signal.butter(4, list(band), btype='bandpass', fs=FS, output='sos')
    
    # Filtrar cada época na banda desejada e extrair sinal analítico via Hilbert
    im_diffs = np.empty(n_ep)
    for i in range(n_ep):
        ep_filt = signal.sosfiltfilt(sos, epochs[i, :, :], axis=0)
        analytic = signal.hilbert(ep_filt, axis=0)
        # Diferença de fase cruzada entre F3 e F4
        cross = analytic[:, 0] * np.conj(analytic[:, 1])
        im_diffs[i] = np.mean(np.imag(cross))
        
    # wPLI empírico real = |E[Im]| / E[|Im|]
    wpli_real = np.abs(np.mean(im_diffs)) / (np.mean(np.abs(im_diffs)) + 1e-12)
    
    # Teste de Surrogates (Hipótese Nula: Acoplamento de fase espúrio/aleatório)
    surr_wplis = np.empty(n_surrogates)
    for s in range(n_surrogates):
        # Inverte aleatoriamente o sinal de imaginação ou permuta épocas
        signs = np.random.choice([-1.0, 1.0], size=n_ep)
        im_surr = im_diffs * signs
        surr_wplis[s] = np.abs(np.mean(im_surr)) / (np.mean(np.abs(im_surr)) + 1e-12)
        
    # P-valor empírico e Z-score
    p_value = (np.sum(surr_wplis >= wpli_real) + 1.0) / (n_surrogates + 1.0)
    z_score = (wpli_real - np.mean(surr_wplis)) / (np.std(surr_wplis) + 1e-12)
    
    return wpli_real, p_value, z_score, np.percentile(surr_wplis, 95)

def main():
    base_file = os.path.join(DATA_DIR, "OpenBCI-RAW-2023-11-28_15-56-59.txt")
    long_file = os.path.join(DATA_DIR, "OpenBCI-RAW-2023-11-28_16-04-54.txt")
    s3_file   = os.path.join(DATA_DIR, "OpenBCI-RAW-2023-11-28_16-36-02.txt")
    s4_file   = os.path.join(DATA_DIR, "OpenBCI-RAW-2023-11-28_16-44-15.txt")
    
    # Carregar dados com duplo notch (60 Hz + 20 Hz de aliasing)
    clean_base = load_and_clean_eeg(base_file)
    clean_long = load_and_clean_eeg(long_file)
    clean_s3   = load_and_clean_eeg(s3_file)
    clean_s4   = load_and_clean_eeg(s4_file)
    
    block_samples = 5 * 60 * FS
    conditions = [
        ("Linha de Base (15:56)", clean_base[:block_samples, :]),
        ("BB Bloco 1 (0-5 min)", clean_long[0:block_samples, :]),
        ("BB Bloco 2 (5-10 min)", clean_long[block_samples:2*block_samples, :]),
        ("BB Bloco 3 (10-15 min)", clean_long[2*block_samples:3*block_samples, :]),
        ("BB Bloco 4 (15-20 min)", clean_long[3*block_samples:4*block_samples, :]),
        ("Pós-Estimulação 1 (16:36)", clean_s3[:block_samples, :]),
        ("Pós-Estimulação 2 (16:44)", clean_s4[:block_samples, :]),
    ]
    
    results = []
    all_faa_data = []
    
    for label, raw_block in conditions:
        epochs, n_rej, n_tot = epoch_data_with_artifact_rejection(raw_block, epoch_sec=2.0, threshold_uV=150.0)
        
        # 1. FAA Relativo e Bootstrap CI
        faa_eps, f3_alphas, f4_alphas = compute_epoch_psd_and_relative_faa(epochs)
        mean_faa, ci_low, ci_high, std_faa = bootstrap_ci(faa_eps)
        
        # 2. wPLI na banda Alfa com teste de Surrogates
        wpli_val, p_val, z_val, null_thresh_95 = compute_wpli_with_surrogates(epochs, band=(8.0, 13.0), n_surrogates=1000)
        
        # 3. Potência Relativa Beta (após remoção do aliasing de 20 Hz)
        f_ep, psd_all = signal.welch(epochs, fs=FS, nperseg=epochs.shape[1], axis=1)
        mean_psd = np.mean(psd_all, axis=0) # Média entre épocas: (freqs, channels)
        
        mask_tot = (f_ep >= 1.0) & (f_ep <= 45.0)
        mask_beta = (f_ep >= 13.0) & (f_ep <= 30.0)
        trapz = getattr(np, 'trapezoid', getattr(np, 'trapz', None))
        tot_pwr = trapz(mean_psd[mask_tot, :], f_ep[mask_tot], axis=0)
        beta_pwr = trapz(mean_psd[mask_beta, :], f_ep[mask_beta], axis=0)
        rel_beta_ch4 = (beta_pwr[3] / tot_pwr[3]) * 100.0 # CH4 (P4/O2)
        rel_beta_f4 = (beta_pwr[1] / tot_pwr[1]) * 100.0  # F4
        
        results.append({
            'Condicao': label,
            'Epocas_Validas': len(epochs),
            'Epocas_Rejeitadas': n_rej,
            'FAA_Rel_Mean': mean_faa,
            'FAA_CI_Lower': ci_low,
            'FAA_CI_Upper': ci_high,
            'FAA_Std': std_faa,
            'wPLI_Alpha_F3_F4': wpli_val,
            'wPLI_pValue': p_val,
            'wPLI_Zscore': z_val,
            'wPLI_Null_95': null_thresh_95,
            'Rel_Beta_CH4_Clean': rel_beta_ch4,
            'Rel_Beta_F4_Clean': rel_beta_f4
        })
        
        all_faa_data.append((label, faa_eps))
        
    df_res = pd.DataFrame(results)
    
    # Calcular Effect Size (Cohen's d) do FAA com desvio padrão amostral rigoroso (ddof=1)
    base_faa = all_faa_data[0][1]
    n_base = len(base_faa)
    base_mean = np.mean(base_faa)
    base_std = np.std(base_faa, ddof=1)
    
    cohens_d = []
    for label, faa_eps in all_faa_data:
        n_cur = len(faa_eps)
        m_diff = np.mean(faa_eps) - base_mean
        cur_std = np.std(faa_eps, ddof=1)
        pooled_std = np.sqrt(((n_base - 1) * (base_std**2) + (n_cur - 1) * (cur_std**2)) / (n_base + n_cur - 2))
        d_val = m_diff / pooled_std
        cohens_d.append(d_val)
        
    df_res['Cohens_d_vs_Baseline'] = cohens_d
    
    # Salvar tabela com rigor estatístico
    csv_path = os.path.join(OUTPUT_DIR, "tables", "tabela_estatistica_alta_confiabilidade.csv")
    df_res.to_csv(csv_path, index=False)
    print("Resultados salvos em:", csv_path)
    print(df_res[['Condicao', 'FAA_Rel_Mean', 'FAA_CI_Lower', 'FAA_CI_Upper', 'Cohens_d_vs_Baseline', 'wPLI_Alpha_F3_F4', 'wPLI_pValue', 'Rel_Beta_CH4_Clean']].to_string())
    
    # -------------------------------------------------------------
    # Geração dos Gráficos com Barras de Erro e Nível Nulo
    # -------------------------------------------------------------
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))
    x = np.arange(len(df_res))
    cond_labels = [c.replace(" (15:56)", "").replace(" (16:36)", "").replace(" (16:44)", "") for c in df_res['Condicao']]
    
    # Painel 1: Assimetria Alfa Frontal Relativa com IC 95% (Bootstrap)
    yerr_lower = df_res['FAA_Rel_Mean'] - df_res['FAA_CI_Lower']
    yerr_upper = df_res['FAA_CI_Upper'] - df_res['FAA_Rel_Mean']
    
    bars1 = axes[0].bar(x, df_res['FAA_Rel_Mean'], yerr=[yerr_lower, yerr_upper], 
                        capsize=5, color=['#718096', '#4299e1', '#3182ce', '#2b6cb0', '#2c5282', '#9f7aea', '#805ad5'],
                        edgecolor='black', linewidth=0.8, alpha=0.9)
    axes[0].axhline(0, color='black', linestyle='--', linewidth=1, label='Simetria Hemisférica Perfeita (0.0)')
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(cond_labels, rotation=35, ha='right', fontsize=9)
    axes[0].set_title('A. Assimetria Alfa Frontal Relativa (FAA)\nIC 95% Bootstrap (Normalizada Intra-Canal)', fontweight='bold')
    axes[0].set_ylabel('FAA Relativa = ln(Rel_Alpha F4) - ln(Rel_Alpha F3)')
    axes[0].grid(True, linestyle='--', alpha=0.5, axis='y')
    axes[0].legend(loc='lower left', fontsize=8.5)
    
    # Painel 2: Conectividade Funcional wPLI F3-F4 vs. Distribuição Nula de Surrogates
    axes[1].bar(x, df_res['wPLI_Alpha_F3_F4'], color='#dd6b20', edgecolor='black', linewidth=0.8, alpha=0.85, label='wPLI Observado (F3-F4 Alfa)')
    axes[1].plot(x, df_res['wPLI_Null_95'], color='red', linestyle='--', marker='o', linewidth=2, label='Limiar Nulo P95 (1000 Surrogates)')
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(cond_labels, rotation=35, ha='right', fontsize=9)
    axes[1].set_title('B. Sincronização de Fase wPLI (F3–F4 Alfa)\nTeste de Significância contra 1000 Surrogates', fontweight='bold')
    axes[1].set_ylabel('wPLI (Inter-Épocas)')
    axes[1].grid(True, linestyle='--', alpha=0.5, axis='y')
    axes[1].legend(loc='upper left', fontsize=8.5)
    
    # Anotação de significância estatística
    for i, p in enumerate(df_res['wPLI_pValue']):
        txt = "***" if p < 0.001 else ("**" if p < 0.01 else ("*" if p < 0.05 else "n.s."))
        axes[1].text(i, df_res['wPLI_Alpha_F3_F4'][i] + 0.03, txt, ha='center', fontweight='bold', color='black')
        
    # Painel 3: Banda Beta Limpa (Sem Aliasing de 180 Hz)
    axes[2].bar(x, df_res['Rel_Beta_CH4_Clean'], color='#38a169', edgecolor='black', linewidth=0.8, alpha=0.85, label='CH4 (P4/O2) - Limpo')
    axes[2].plot(x, df_res['Rel_Beta_F4_Clean'], color='#2f855a', marker='s', linewidth=2, label='F4 (CH2) - Limpo')
    axes[2].axhline(6.0, color='gray', linestyle=':', label='Patamar Basal Fisiológico (~6%)')
    axes[2].set_xticks(x)
    axes[2].set_xticklabels(cond_labels, rotation=35, ha='right', fontsize=9)
    axes[2].set_title('C. Potência Relativa da Banda Beta (%)\nApós Remoção do Aliasing da Rede (180 Hz -> 20 Hz)', fontweight='bold')
    axes[2].set_ylabel('Potência Beta Relativa (%)')
    axes[2].grid(True, linestyle='--', alpha=0.5, axis='y')
    axes[2].legend(loc='upper right', fontsize=8.5)
    
    plt.tight_layout()
    plot_path = os.path.join(OUTPUT_DIR, "figures", "validacao_estatistica_alta_confiabilidade.png")
    plt.savefig(plot_path, dpi=180, bbox_inches='tight')
    plt.close()
    print("Figura de validação estatística salva em:", plot_path)

if __name__ == "__main__":
    main()
