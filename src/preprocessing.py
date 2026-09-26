"""
Módulo de Pré-Processamento e Rejeição de Artefatos de Sinais de EEG.
Projetado para dados OpenBCI Ganglion (4 canais, 200 Hz) com Ultracortex Mark IV.
"""

import numpy as np
import pandas as pd
from scipy import signal

FS = 200 # Taxa de amostragem padrão da Ganglion (Hz)
CHANNELS = ['F3', 'F4', 'P3/O1', 'P4/O2']

def load_raw_openbci(filepath):
    """Carrega as colunas EXG 0 a 3 de arquivos CSV salvos pelo OpenBCI GUI."""
    df = pd.read_csv(filepath, skiprows=1, header=None, comment='%', usecols=[1, 2, 3, 4])
    return df.to_numpy(dtype=float)

def apply_filters(raw_data, fs=FS, notch_grid=60.0, filter_alias_harmonics=True, bp_band=(1.0, 45.0)):
    """
    Aplica cadeia de filtragem de fase zero (bidirecional):
    1. Detrend linear para eliminação de offsets DC e drift de polarização.
    2. Filtro Notch na frequência da rede elétrica (60 Hz no Brasil, Q=30).
    3. Filtros Notch de Aliasing:
       - 180 Hz (3º harmônico) dobra em |180 - 200| = 20 Hz sob fs=200 Hz.
       - 120 Hz (2º harmônico) dobra em |120 - 200| = 80 Hz sob fs=200 Hz.
    4. Filtro Passa-Faixa Butterworth 1-45 Hz de 4ª ordem em Seções de Segunda Ordem (SOS).
    """
    # 1. Detrend linear
    d = signal.detrend(raw_data, axis=0)
    
    # 2. Notch 60 Hz
    b60, a60 = signal.iirnotch(notch_grid, 30.0, fs=fs)
    d = signal.filtfilt(b60, a60, d, axis=0)
    
    # 3. Notches de Aliasing (rebatimento de harmônicos superiores à Nyquist de 100 Hz)
    if filter_alias_harmonics:
        b20, a20 = signal.iirnotch(20.0, 30.0, fs=fs)
        d = signal.filtfilt(b20, a20, d, axis=0)
        b80, a80 = signal.iirnotch(80.0, 30.0, fs=fs)
        d = signal.filtfilt(b80, a80, d, axis=0)
        
    # 4. Passa-faixa 1-45 Hz (SOS fase zero)
    sos = signal.butter(4, list(bp_band), btype='bandpass', fs=fs, output='sos')
    clean = signal.sosfiltfilt(sos, d, axis=0)
    return clean

def segment_epochs(clean_data, fs=FS, epoch_sec=2.0, threshold_frontal_p2p=500.0, min_p2p=2.0):
    """
    Segmenta dados contínuos em épocas não-sobrepostas e rejeita épocas com artefatos extremos:
    - Rejeita se amplitude pico-a-pico frontal (F3 ou F4) > threshold_frontal_p2p (ex.: espasmo muscular, movimento brusco).
    - Rejeita se amplitude pico-a-pico frontal < min_p2p (desconexão ou flatline).
    """
    ep_len = int(epoch_sec * fs)
    n_epochs = clean_data.shape[0] // ep_len
    
    valid_epochs = []
    valid_indices = []
    
    for i in range(n_epochs):
        ep = clean_data[i*ep_len : (i+1)*ep_len, :]
        p2p_frontal = np.ptp(ep[:, :2], axis=0)
        
        if np.any(p2p_frontal > threshold_frontal_p2p) or np.any(p2p_frontal < min_p2p):
            continue
            
        valid_epochs.append(ep)
        valid_indices.append(i)
        
    return np.array(valid_epochs), valid_indices
