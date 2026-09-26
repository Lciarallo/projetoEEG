import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import signal

DATA_DIR = "/home/lciarallo/projetoEEG/X"
OUTPUT_DIR = "/home/lciarallo/projetoEEG/resultados"
os.makedirs(OUTPUT_DIR, exist_ok=True)

FS = 200
CHANNELS = ['CH1', 'CH2', 'CH3', 'CH4']

def load_clean_data(filepath):
    df = pd.read_csv(filepath, skiprows=1, header=None, comment='%', usecols=[1, 2, 3, 4])
    raw = df.to_numpy(dtype=float)
    # Detrend
    d = signal.detrend(raw, axis=0)
    # Notch 60 Hz
    b, a = signal.iirnotch(60.0, 30.0, fs=FS)
    d = signal.filtfilt(b, a, d, axis=0)
    # Bandpass 1-45 Hz
    sos = signal.butter(4, [1.0, 45.0], btype='bandpass', fs=FS, output='sos')
    clean = signal.sosfiltfilt(sos, d, axis=0)
    return clean

def compute_wpli_matrix(data, fs=200, band=(8.0, 13.0)):
    sos = signal.butter(4, list(band), btype='bandpass', fs=fs, output='sos')
    filtered = signal.sosfiltfilt(sos, data, axis=0)
    analytic = signal.hilbert(filtered, axis=0)
    
    n_ch = data.shape[1]
    mat = np.zeros((n_ch, n_ch))
    for i in range(n_ch):
        for j in range(n_ch):
            if i == j:
                mat[i, j] = 1.0
            elif j > i:
                cross = analytic[:, i] * np.conj(analytic[:, j])
                im = np.imag(cross)
                wpli = np.abs(np.mean(im)) / (np.mean(np.abs(im)) + 1e-12)
                mat[i, j] = wpli
                mat[j, i] = wpli
    return mat

def main():
    base_file = os.path.join(DATA_DIR, "OpenBCI-RAW-2023-11-28_15-56-59.txt")
    long_file = os.path.join(DATA_DIR, "OpenBCI-RAW-2023-11-28_16-04-54.txt")
    
    clean_base = load_clean_data(base_file)
    clean_long = load_clean_data(long_file)
    
    # Blocos de 5 min da sessão longa
    block_len = 5 * 60 * FS
    blocks = [
        ("Linha de Base (15:56)", clean_base[:block_len, :]),
        ("BB Bloco 1 (0-5 min)", clean_long[0:block_len, :]),
        ("BB Bloco 2 (5-10 min)", clean_long[block_len:2*block_len, :]),
        ("BB Bloco 3 (10-15 min)", clean_long[2*block_len:3*block_len, :]),
        ("BB Bloco 4 (15-20 min)", clean_long[3*block_len:4*block_len, :]),
    ]
    
    # 1. Matrizes de Conectividade wPLI na Banda Alfa
    fig, axes = plt.subplots(1, 5, figsize=(18, 4))
    for idx, (label, b_data) in enumerate(blocks):
        wpli_mat = compute_wpli_matrix(b_data, fs=FS, band=(8.0, 13.0))
        im = axes[idx].imshow(wpli_mat, cmap='YlOrRd', vmin=0, vmax=0.6)
        axes[idx].set_title(label, fontsize=11, fontweight='bold')
        axes[idx].set_xticks(range(4))
        axes[idx].set_yticks(range(4))
        axes[idx].set_xticklabels(CHANNELS)
        axes[idx].set_yticklabels(CHANNELS)
        for i in range(4):
            for j in range(4):
                if i != j:
                    axes[idx].text(j, i, f"{wpli_mat[i, j]:.2f}", ha='center', va='center', color='black', fontsize=9)
                    
    fig.subplots_adjust(right=0.9)
    cbar_ax = fig.add_axes([0.92, 0.2, 0.015, 0.6])
    fig.colorbar(im, cax=cbar_ax, label='wPLI (Imune a Condução de Volume)')
    fig.suptitle('Evolução da Conectividade Funcional wPLI (Banda Alfa 8-13 Hz): Repouso vs. Batimento Binaural', fontsize=13, y=1.03)
    
    out_heatmaps = os.path.join(OUTPUT_DIR, "conectividade_wpli_evolucao.png")
    plt.savefig(out_heatmaps, dpi=160, bbox_inches='tight')
    plt.close()
    print("Salvo:", out_heatmaps)
    
    # 2. Gráfico de Dinâmica de Potência Relativa (Beta e Theta) e Razão Alfa/Beta
    df_dyn = []
    for label, b_data in blocks:
        f, psd = signal.welch(b_data, fs=FS, nperseg=400, noverlap=200, axis=0)
        tot_mask = (f >= 1.0) & (f <= 45.0)
        tot_pwr = np.trapezoid(psd[tot_mask, :], f[tot_mask], axis=0) if hasattr(np, 'trapezoid') else np.trapz(psd[tot_mask, :], f[tot_mask], axis=0)
        
        trapz_func = getattr(np, 'trapezoid', getattr(np, 'trapz', None))
        theta_pwr = trapz_func(psd[(f >= 4.0) & (f <= 8.0), :], f[(f >= 4.0) & (f <= 8.0)], axis=0)
        alpha_pwr = trapz_func(psd[(f >= 8.0) & (f <= 13.0), :], f[(f >= 8.0) & (f <= 13.0)], axis=0)
        beta_pwr = trapz_func(psd[(f >= 13.0) & (f <= 30.0), :], f[(f >= 13.0) & (f <= 30.0)], axis=0)
        
        for c in range(4):
            df_dyn.append({
                'Condicao': label,
                'Canal': CHANNELS[c],
                'Rel_Theta': (theta_pwr[c] / tot_pwr[c]) * 100,
                'Rel_Alpha': (alpha_pwr[c] / tot_pwr[c]) * 100,
                'Rel_Beta': (beta_pwr[c] / tot_pwr[c]) * 100,
                'Alpha_Beta': alpha_pwr[c] / (beta_pwr[c] + 1e-12),
                'Theta_Beta': theta_pwr[c] / (beta_pwr[c] + 1e-12)
            })
    df_dyn = pd.DataFrame(df_dyn)
    
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    condicoes = [b[0] for b in blocks]
    
    # Subplot 1: Potência Relativa Beta
    for c in CHANNELS:
        sub = df_dyn[df_dyn['Canal'] == c]
        axes[0].plot(sub['Condicao'], sub['Rel_Beta'], marker='o', label=c, linewidth=2)
    axes[0].set_title('Potência Relativa Beta (%)', fontweight='bold')
    axes[0].set_ylabel('Beta (%)')
    axes[0].tick_params(axis='x', rotation=25)
    axes[0].grid(True, linestyle='--', alpha=0.5)
    axes[0].legend()
    
    # Subplot 2: Razão Alfa/Beta
    for c in CHANNELS:
        sub = df_dyn[df_dyn['Canal'] == c]
        axes[1].plot(sub['Condicao'], sub['Alpha_Beta'], marker='s', label=c, linewidth=2)
    axes[1].set_title('Razão Alfa/Beta (Índice de Alerta)', fontweight='bold')
    axes[1].set_ylabel('Alfa / Beta')
    axes[1].tick_params(axis='x', rotation=25)
    axes[1].grid(True, linestyle='--', alpha=0.5)
    axes[1].legend()
    
    # Subplot 3: Razão Teta/Beta (TBR)
    for c in CHANNELS:
        sub = df_dyn[df_dyn['Canal'] == c]
        axes[2].plot(sub['Condicao'], sub['Theta_Beta'], marker='^', label=c, linewidth=2)
    axes[2].set_title('Razão Teta/Beta (TBR - Fadiga/Desatenção)', fontweight='bold')
    axes[2].set_ylabel('Teta / Beta')
    axes[2].tick_params(axis='x', rotation=25)
    axes[2].grid(True, linestyle='--', alpha=0.5)
    axes[2].legend()
    
    fig.suptitle('Trajetória Neurofisiológica durante Estimulação com Batimentos Binaurais', fontsize=13, y=1.02)
    plt.tight_layout()
    out_trajectory = os.path.join(OUTPUT_DIR, "trajetoria_neurofisiologica_bb.png")
    plt.savefig(out_trajectory, dpi=160, bbox_inches='tight')
    plt.close()
    print("Salvo:", out_trajectory)

if __name__ == "__main__":
    main()
