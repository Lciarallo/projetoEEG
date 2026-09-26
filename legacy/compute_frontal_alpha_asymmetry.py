import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import signal

DATA_DIR = "/home/lciarallo/projetoEEG/X"
OUTPUT_DIR = "/home/lciarallo/projetoEEG/resultados"
os.makedirs(OUTPUT_DIR, exist_ok=True)

FS = 200
# Mapeamento com base nas fotos do Ultracortex Mark IV
# CH1: F3 (Frontal Esquerdo)
# CH2: F4 (Frontal Direito)
# CH3: P3 / O1 (Posterior Esquerdo)
# CH4: P4 / O2 (Posterior Direito)

def load_and_preprocess(filepath):
    df = pd.read_csv(filepath, skiprows=1, header=None, comment='%', usecols=[1, 2, 3, 4])
    raw = df.to_numpy(dtype=float)
    # 1. Detrend
    d = signal.detrend(raw, axis=0)
    # 2. Notch 60 Hz
    b, a = signal.iirnotch(60.0, 30.0, fs=FS)
    d = signal.filtfilt(b, a, d, axis=0)
    # 3. Bandpass 1-45 Hz (SOS de fase zero)
    sos = signal.butter(4, [1.0, 45.0], btype='bandpass', fs=FS, output='sos')
    clean = signal.sosfiltfilt(sos, d, axis=0)
    return clean

def compute_alpha_power(data_segment, fs=200):
    """Calcula potência absoluta na banda Alfa (8-13 Hz) usando Welch."""
    f, psd = signal.welch(data_segment, fs=fs, nperseg=min(len(data_segment), 400), axis=0)
    alpha_mask = (f >= 8.0) & (f <= 13.0)
    trapz_func = getattr(np, 'trapezoid', getattr(np, 'trapz', None))
    alpha_pwr = trapz_func(psd[alpha_mask, :], f[alpha_mask], axis=0)
    return alpha_pwr

def compute_continuous_faa(clean_data, fs=200, window_sec=6.0, step_sec=1.5):
    """
    Calcula Assimetria Alfa Frontal (FAA) contínua em janelas deslizantes:
    FAA = ln(Alpha_F4) - ln(Alpha_F3)
    F3 = CH1 (índice 0), F4 = CH2 (índice 1)
    """
    win_len = int(window_sec * fs)
    step_len = int(step_sec * fs)
    n_samples = clean_data.shape[0]
    
    times = []
    faa_vals = []
    alpha_f3_list = []
    alpha_f4_list = []
    
    for start in range(0, n_samples - win_len, step_len):
        end = start + win_len
        seg = clean_data[start:end, :]
        pwr = compute_alpha_power(seg, fs=fs)
        
        # CH1 = F3, CH2 = F4
        p_f3 = pwr[0]
        p_f4 = pwr[1]
        
        # Evitar log de zero
        if p_f3 > 1e-6 and p_f4 > 1e-6:
            faa = np.log(p_f4) - np.log(p_f3)
            times.append((start + end) / (2.0 * fs)) # centro da janela em segundos
            faa_vals.append(faa)
            alpha_f3_list.append(p_f3)
            alpha_f4_list.append(p_f4)
            
    return np.array(times), np.array(faa_vals), np.array(alpha_f3_list), np.array(alpha_f4_list)

def main():
    base_file = os.path.join(DATA_DIR, "OpenBCI-RAW-2023-11-28_15-56-59.txt")
    long_file = os.path.join(DATA_DIR, "OpenBCI-RAW-2023-11-28_16-04-54.txt")
    s3_file   = os.path.join(DATA_DIR, "OpenBCI-RAW-2023-11-28_16-36-02.txt")
    s4_file   = os.path.join(DATA_DIR, "OpenBCI-RAW-2023-11-28_16-44-15.txt")
    
    clean_base = load_and_preprocess(base_file)
    clean_long = load_and_preprocess(long_file)
    clean_s3   = load_and_preprocess(s3_file)
    clean_s4   = load_and_preprocess(s4_file)
    
    # 1. FAA em Janela Contínua para a Sessão Longa de Batimentos Binaurais (16:04:54)
    t_long, faa_long, a_f3_long, a_f4_long = compute_continuous_faa(clean_long, fs=FS, window_sec=8.0, step_sec=2.0)
    
    # Também para a Linha de Base (15:56:59)
    t_base, faa_base, _, _ = compute_continuous_faa(clean_base, fs=FS, window_sec=8.0, step_sec=2.0)
    
    # 2. FAA Médio por Condição e por Bloco de 5 min
    block_len = 5 * 60 * FS
    conditions = [
        ("Linha de Base (15:56)", clean_base),
        ("BB Bloco 1 (0-5 min)", clean_long[0:block_len, :]),
        ("BB Bloco 2 (5-10 min)", clean_long[block_len:2*block_len, :]),
        ("BB Bloco 3 (10-15 min)", clean_long[2*block_len:3*block_len, :]),
        ("BB Bloco 4 (15-20 min)", clean_long[3*block_len:4*block_len, :]),
        ("Sessão Pós 1 (16:36)", clean_s3),
        ("Sessão Pós 2 (16:44)", clean_s4),
    ]
    
    faa_summary = []
    for label, d in conditions:
        pwr = compute_alpha_power(d, fs=FS)
        p_f3 = pwr[0]
        p_f4 = pwr[1]
        faa = np.log(p_f4) - np.log(p_f3)
        faa_summary.append({
            'Condicao': label,
            'Alpha_F3': p_f3,
            'Alpha_F4': p_f4,
            'FAA': faa,
            'Estado': 'Predomínio Esquerdo (Aproximação/Relaxamento)' if faa > 0 else 'Predomínio Direito (Evitação/Alerta)'
        })
        
    df_faa = pd.DataFrame(faa_summary)
    csv_out = os.path.join(OUTPUT_DIR, "tabela_assimetria_alfa_frontal.csv")
    df_faa.to_csv(csv_out, index=False)
    print("Tabela FAA salva em:", csv_out)
    print(df_faa[['Condicao', 'Alpha_F3', 'Alpha_F4', 'FAA', 'Estado']].to_string())
    
    # 3. Gerar Figura com Análise Completa de FAA
    fig = plt.figure(figsize=(15, 9))
    gs = fig.add_gridspec(2, 2, hspace=0.35, wspace=0.25)
    
    # Painel 1: Trajetória Contínua de FAA durante os 20 min de Batimentos Binaurais
    ax1 = fig.add_subplot(gs[0, :])
    ax1.plot(t_long / 60.0, faa_long, color='#3182ce', alpha=0.35, label='FAA Janela 8s')
    # Média móvel suave
    kernel_size = 15
    faa_smooth = np.convolve(faa_long, np.ones(kernel_size)/kernel_size, mode='valid')
    t_smooth = t_long[:len(faa_smooth)] / 60.0
    ax1.plot(t_smooth, faa_smooth, color='#2b6cb0', linewidth=2.5, label='Tendência Suavizada (Média Móvel)')
    
    # Linha zero de simetria
    ax1.axhline(0, color='gray', linestyle='--', linewidth=1.2)
    ax1.axhline(np.mean(faa_base), color='#e53e3e', linestyle=':', linewidth=2, label=f'Média Linha de Base ({np.mean(faa_base):.2f})')
    
    # Sombreamento das zonas
    ax1.fill_between(t_long / 60.0, 0, np.maximum(faa_long, 0), color='#48bb78', alpha=0.12, label='Zona Positiva (Relaxamento/Aproximação)')
    ax1.fill_between(t_long / 60.0, np.minimum(faa_long, 0), 0, color='#f56565', alpha=0.12, label='Zona Negativa (Evitação/Alerta)')
    
    ax1.set_title('Evolução Temporal Contínua da Assimetria Alfa Frontal (F4 vs. F3) durante Estimulação', fontsize=12, fontweight='bold')
    ax1.set_xlabel('Tempo (minutos)')
    ax1.set_ylabel('FAA = ln(Alfa F4) - ln(Alfa F3)')
    ax1.set_xlim(0, t_long[-1]/60.0)
    ax1.legend(loc='upper right', fontsize=9, ncol=2)
    ax1.grid(True, linestyle='--', alpha=0.5)
    
    # Painel 2: Gráfico de Barras Comparativo entre Todas as Fases
    ax2 = fig.add_subplot(gs[1, 0])
    cores = ['#718096', '#4299e1', '#3182ce', '#2b6cb0', '#2c5282', '#9f7aea', '#805ad5']
    bars = ax2.bar(range(len(df_faa)), df_faa['FAA'], color=cores, width=0.6, edgecolor='black', linewidth=0.8)
    ax2.axhline(0, color='black', linestyle='--', linewidth=1)
    ax2.set_xticks(range(len(df_faa)))
    ax2.set_xticklabels(df_faa['Condicao'], rotation=30, ha='right', fontsize=9)
    ax2.set_title('FAA Médio por Condição e Bloco', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Índice FAA de Davidson')
    ax2.grid(True, linestyle='--', alpha=0.4, axis='y')
    
    # Valores em cima das barras
    for bar in bars:
        h = bar.get_height()
        va = 'bottom' if h >= 0 else 'top'
        ax2.annotate(f"{h:.2f}",
                     xy=(bar.get_x() + bar.get_width() / 2, h),
                     xytext=(0, 3 if h >= 0 else -10),
                     textcoords="offset points",
                     ha='center', va=va, fontsize=8.5, fontweight='bold')
                     
    # Painel 3: Comparação de Potência Absoluta Alfa em F3 vs. F4
    ax3 = fig.add_subplot(gs[1, 1])
    x = np.arange(len(df_faa))
    w = 0.35
    ax3.bar(x - w/2, df_faa['Alpha_F3'], width=w, label='F3 (Esquerdo - CH1)', color='#ed8936', edgecolor='black', linewidth=0.8)
    ax3.bar(x + w/2, df_faa['Alpha_F4'], width=w, label='F4 (Direito - CH2)', color='#4299e1', edgecolor='black', linewidth=0.8)
    ax3.set_xticks(x)
    ax3.set_xticklabels(df_faa['Condicao'], rotation=30, ha='right', fontsize=9)
    ax3.set_title('Potência Absoluta Alfa (F3 vs. F4)', fontsize=12, fontweight='bold')
    ax3.set_ylabel(r'Potência Alfa ($\mu V^2$)')
    ax3.legend(fontsize=9)
    ax3.grid(True, linestyle='--', alpha=0.4, axis='y')
    
    plot_path = os.path.join(OUTPUT_DIR, "assimetria_alfa_frontal_faa.png")
    plt.savefig(plot_path, dpi=160, bbox_inches='tight')
    plt.close()
    print("Gráfico de FAA salvo em:", plot_path)

if __name__ == "__main__":
    main()
