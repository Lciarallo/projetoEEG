import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import signal
from specparam import SpectralModel

# Configurações do matplotlib para relatórios científicos
plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
plt.rcParams['font.sans-serif'] = 'DejaVu Sans'
plt.rcParams['figure.titlesize'] = 14
plt.rcParams['axes.titlesize'] = 12

DATA_DIR = "/home/lciarallo/projetoEEG/X"
OUTPUT_DIR = "/home/lciarallo/projetoEEG/resultados"
os.makedirs(OUTPUT_DIR, exist_ok=True)

SAMPLING_RATE = 200 # Hz (OpenBCI Ganglion)
CHANNELS = ['CH1', 'CH2', 'CH3', 'CH4']

BANDS = {
    'Delta': (1.0, 4.0),
    'Theta': (4.0, 8.0),
    'Alpha': (8.0, 13.0),
    'Beta': (13.0, 30.0),
    'Gamma': (30.0, 45.0)
}

def load_openbci_raw(filepath):
    """Carrega dados brutos do OpenBCI ignorando comentários e selecionando colunas EXG 0 a 3."""
    # Ler arquivo pulando cabeçalho
    df = pd.read_csv(filepath, skiprows=1, header=None, comment='%', usecols=[1, 2, 3, 4])
    data = df.to_numpy(dtype=float)
    return data

def preprocess_eeg(data, fs=200):
    """
    Pré-processamento de alta fidelidade para EEG do OpenBCI:
    1. Detrend linear para remover deriva de linha de base.
    2. Filtro Notch em 60 Hz (Q=30) de fase zero para rede elétrica brasileira.
    3. Filtro Passa-Faixa Butterworth 1-45 Hz usando Seções de Segunda Ordem (SOS) e sosfiltfilt (fase zero).
    """
    # 1. Detrend
    data_detrended = signal.detrend(data, axis=0)
    
    # 2. Notch 60 Hz
    b_notch, a_notch = signal.iirnotch(w0=60.0, Q=30, fs=fs)
    data_notched = signal.filtfilt(b_notch, a_notch, data_detrended, axis=0)
    
    # 3. Bandpass 1-45 Hz (SOS - 4ª ordem)
    sos = signal.butter(N=4, Wn=[1.0, 45.0], btype='bandpass', fs=fs, output='sos')
    data_clean = signal.sosfiltfilt(sos, data_notched, axis=0)
    
    return data_clean

def compute_welch_psd(data, fs=200, nperseg=400):
    """Calcula Densidade Espectral de Potência usando janela de Welch de 2 segundos (400 pts a 200 Hz)."""
    freqs, psd = signal.welch(data, fs=fs, nperseg=nperseg, noverlap=nperseg//2, axis=0, window='hann')
    return freqs, psd

def extract_band_metrics(freqs, psd):
    """Calcula potência absoluta, relativa e razões funcionais para cada canal."""
    metrics = []
    # Banda total considerada: 1 a 45 Hz
    total_mask = (freqs >= 1.0) & (freqs <= 45.0)
    
    for ch_idx, ch_name in enumerate(CHANNELS):
        ch_psd = psd[:, ch_idx]
        total_pwr = np.trapezoid(ch_psd[total_mask], freqs[total_mask]) if hasattr(np, 'trapezoid') else np.trapz(ch_psd[total_mask], freqs[total_mask])
        
        row = {'Canal': ch_name, 'Total_Power': total_pwr}
        for b_name, (low, high) in BANDS.items():
            b_mask = (freqs >= low) & (freqs <= high)
            trapz_func = getattr(np, 'trapezoid', getattr(np, 'trapz', None))
            abs_pwr = trapz_func(ch_psd[b_mask], freqs[b_mask])
            row[f'Abs_{b_name}'] = abs_pwr
            row[f'Rel_{b_name}'] = (abs_pwr / total_pwr) * 100.0
            
        row['Alpha_Beta_Ratio'] = row['Abs_Alpha'] / (row['Abs_Beta'] + 1e-12)
        row['Theta_Beta_Ratio'] = row['Abs_Theta'] / (row['Abs_Beta'] + 1e-12)
        metrics.append(row)
        
    return pd.DataFrame(metrics)

def fit_specparam_fooof(freqs, psd_ch):
    """
    Decomposição aperiódica 1/f vs periódica via specparam (FOOOF).
    Retorna o expoente aperiódico chi (índice do balanço Excitação/Inibição) e r2 do ajuste.
    """
    fm = SpectralModel(peak_width_limits=[1.0, 8.0], max_n_peaks=4, min_peak_height=0.1, verbose=False)
    # Ajuste na faixa 2 a 40 Hz para evitar bordas extremas
    mask = (freqs >= 2.0) & (freqs <= 40.0)
    try:
        fm.fit(freqs[mask], psd_ch[mask])
        ap = fm.results.get_params('aperiodic')
        aperiodic_offset = ap[0]
        aperiodic_exp = ap[1]
        r2 = fm.results.get_metrics('gof_rsquared')
        peaks = fm.results.get_params('peak')
    except Exception as e:
        print(f"Erro no ajuste specparam: {e}")
        aperiodic_exp, aperiodic_offset, r2, peaks = np.nan, np.nan, np.nan, []
    return aperiodic_exp, aperiodic_offset, r2, peaks

def compute_wpli(data, fs=200, band=(8.0, 13.0)):
    """
    Calcula weighted Phase Lag Index (wPLI) entre todos os pares dos 4 canais.
    wPLI é estritamente imune à condução de volume (fase zero eliminada).
    """
    # Filtra estreito na banda de interesse
    sos = signal.butter(N=4, Wn=list(band), btype='bandpass', fs=fs, output='sos')
    filtered = signal.sosfiltfilt(sos, data, axis=0)
    
    # Transformada de Hilbert para obter fase instantânea
    analytic = signal.hilbert(filtered, axis=0)
    
    n_channels = data.shape[1]
    pairs = {}
    for i in range(n_channels):
        for j in range(i + 1, n_channels):
            # Diferença de fase cruzada
            cross_spec = analytic[:, i] * np.conj(analytic[:, j])
            im = np.imag(cross_spec)
            # wPLI = |E[Im]| / E[|Im|]
            wpli = np.abs(np.mean(im)) / (np.mean(np.abs(im)) + 1e-12)
            pairs[f"{CHANNELS[i]}-{CHANNELS[j]}"] = wpli
    return pairs

def main():
    files = sorted(glob.glob(os.path.join(DATA_DIR, "OpenBCI-RAW-*.txt")))
    print(f"Arquivos encontrados ({len(files)}):")
    for f in files:
        print(" -", os.path.basename(f))
        
    all_summary = []
    
    # 1. Processar cada sessão individualmente
    fig_psd, axes_psd = plt.subplots(len(files), 1, figsize=(11, 3.5 * len(files)), sharex=True)
    if len(files) == 1:
        axes_psd = [axes_psd]
        
    for idx, filepath in enumerate(files):
        fname = os.path.basename(filepath)
        tag = fname.replace("OpenBCI-RAW-2023-11-28_", "").replace(".txt", "")
        
        raw = load_openbci_raw(filepath)
        clean = preprocess_eeg(raw, fs=SAMPLING_RATE)
        freqs, psd = compute_welch_psd(clean, fs=SAMPLING_RATE)
        df_bands = extract_band_metrics(freqs, psd)
        
        # Calcular wPLI na banda Alfa (8-13 Hz) e Beta (13-30 Hz)
        wpli_alpha = compute_wpli(clean, fs=SAMPLING_RATE, band=(8.0, 13.0))
        wpli_beta = compute_wpli(clean, fs=SAMPLING_RATE, band=(13.0, 30.0))
        
        # Ajuste SpecParam (FOOOF) por canal
        exponents = []
        for ch in range(4):
            exp_val, off_val, r2_val, _ = fit_specparam_fooof(freqs, psd[:, ch])
            exponents.append(exp_val)
            
        df_bands['Aperiodic_Exponent_1f'] = exponents
        df_bands['Session'] = tag
        all_summary.append(df_bands)
        
        # Plot do PSD da sessão
        ax = axes_psd[idx]
        for ch in range(4):
            ax.semilogy(freqs, psd[:, ch], label=f"{CHANNELS[ch]} (1/f exp={exponents[ch]:.2f})", linewidth=1.5)
        ax.set_title(f"Sessão: {tag} ({raw.shape[0]/SAMPLING_RATE:.1f} s)", fontsize=11, fontweight='bold')
        ax.set_ylabel(r"PSD ($\mu V^2 / Hz$)")
        ax.set_xlim(1, 45)
        ax.axvspan(8, 13, color='orange', alpha=0.15, label='Banda Alfa (8-13 Hz)' if ch==0 else "")
        ax.axvspan(13, 30, color='green', alpha=0.15, label='Banda Beta (13-30 Hz)' if ch==0 else "")
        ax.grid(True, linestyle='--', alpha=0.5)
        ax.legend(loc='upper right', fontsize=9)
        
    axes_psd[-1].set_xlabel("Frequência (Hz)", fontsize=11)
    fig_psd.tight_layout()
    psd_plot_path = os.path.join(OUTPUT_DIR, "comparativo_psd_sessoes.png")
    fig_psd.savefig(psd_plot_path, dpi=150)
    plt.close(fig_psd)
    print(f"Gráfico de PSD salvo em: {psd_plot_path}")
    
    # Salvar tabela com resumo de todas as sessões
    df_all = pd.concat(all_summary, ignore_index=True)
    summary_csv = os.path.join(OUTPUT_DIR, "resumo_metricas_todas_sessoes.csv")
    df_all.to_csv(summary_csv, index=False)
    print(f"Resumo de métricas salvo em: {summary_csv}")
    
    # 2. Análise detalhada da Sessão Longa de 20 minutos (16-04-54)
    # Dividindo em 4 blocos de 5 minutos para simular as 4 fases/estimulações (como Gao et al., 2014)
    long_file = [f for f in files if "16-04-54" in f][0]
    raw_long = load_openbci_raw(long_file)
    clean_long = preprocess_eeg(raw_long, fs=SAMPLING_RATE)
    
    n_samples = clean_long.shape[0]
    block_len = 5 * 60 * SAMPLING_RATE # 5 min = 60000 amostras
    n_blocks = min(4, n_samples // block_len)
    
    block_records = []
    fig_blocks, ax_blocks = plt.subplots(figsize=(10, 5))
    
    for b in range(n_blocks):
        start_idx = b * block_len
        end_idx = (b + 1) * block_len
        block_data = clean_long[start_idx:end_idx, :]
        f_b, psd_b = compute_welch_psd(block_data, fs=SAMPLING_RATE)
        df_b = extract_band_metrics(f_b, psd_b)
        df_b['Bloco'] = f"Min {b*5}-{(b+1)*5}"
        
        # Média dos canais para curva média do bloco
        mean_psd = np.mean(psd_b, axis=1)
        ax_blocks.semilogy(f_b, mean_psd, label=f"Bloco {b+1} (Min {b*5}-{(b+1)*5})", linewidth=1.8)
        block_records.append(df_b)
        
    ax_blocks.set_title("Evolução Temporal da Sessão Longa (Blocos de 5 min) - Média dos 4 Canais", fontsize=12, fontweight='bold')
    ax_blocks.set_xlabel("Frequência (Hz)")
    ax_blocks.set_ylabel(r"PSD Médio ($\mu V^2 / Hz$)")
    ax_blocks.set_xlim(1, 45)
    ax_blocks.axvspan(8, 13, color='orange', alpha=0.15, label='Alfa (8-13 Hz)')
    ax_blocks.axvspan(13, 30, color='green', alpha=0.15, label='Beta (13-30 Hz)')
    ax_blocks.legend()
    ax_blocks.grid(True, linestyle='--', alpha=0.5)
    fig_blocks.tight_layout()
    blocks_plot_path = os.path.join(OUTPUT_DIR, "evolucao_blocos_sessao_longa.png")
    fig_blocks.savefig(blocks_plot_path, dpi=150)
    plt.close(fig_blocks)
    print(f"Gráfico de blocos salvo em: {blocks_plot_path}")
    
    df_blocks = pd.concat(block_records, ignore_index=True)
    blocks_csv = os.path.join(OUTPUT_DIR, "resumo_blocos_sessao_longa.csv")
    df_blocks.to_csv(blocks_csv, index=False)
    print("Processamento concluído com sucesso!")

if __name__ == "__main__":
    main()
