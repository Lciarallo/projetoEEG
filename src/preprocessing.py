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
    columns = pd.read_csv(filepath, comment='%', nrows=0).columns
    normalized = {column.strip(): column for column in columns}
    selected = []
    for index in range(4):
        name = next((normalized[candidate] for candidate in
                     (f'EXG Channel {index}', f'EXG {index}') if candidate in normalized), None)
        if name is None:
            raise ValueError(f'Coluna EXG {index} ausente: {filepath}')
        selected.append(name)
    df = pd.read_csv(filepath, comment='%', usecols=selected)
    raw = df[selected].to_numpy(dtype=float)
    if raw.ndim != 2 or raw.shape[1] != 4 or len(raw) == 0 or not np.isfinite(raw).all():
        raise ValueError(f"Dados EXG vazios ou não finitos: {filepath}")
    return raw

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
    raw_data = np.asarray(raw_data, dtype=float)
    if raw_data.ndim != 2 or len(raw_data) < 32 or not np.isfinite(raw_data).all():
        raise ValueError("A filtragem requer uma matriz finita com pelo menos 32 amostras")
    if not 0 < bp_band[0] < bp_band[1] < fs / 2 or not 0 < notch_grid < fs / 2:
        raise ValueError("Frequências de filtragem fora do intervalo permitido")
    # 1. Detrend linear
    d = signal.detrend(raw_data, axis=0)
    
    # 2. Notch 60 Hz
    b60, a60 = signal.iirnotch(notch_grid, 30.0, fs=fs)
    d = signal.filtfilt(b60, a60, d, axis=0)
    
    # 3. Notches de Aliasing (rebatimento de harmônicos superiores à Nyquist de 100 Hz)
    if filter_alias_harmonics:
        # Hipóteses de aliasing, não identificação da origem do pico.
        for harmonic in (2, 3):
            alias = abs((harmonic * notch_grid + fs / 2) % fs - fs / 2)
            if 0 < alias < fs / 2:
                b, a = signal.iirnotch(alias, 30.0, fs=fs)
                d = signal.filtfilt(b, a, d, axis=0)
        
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
    clean_data = np.asarray(clean_data, dtype=float)
    ep_len = int(epoch_sec * fs)
    if clean_data.ndim != 2 or clean_data.shape[1] < 2 or ep_len < 1:
        raise ValueError("Dimensões ou duração de época inválidas")
    if not 0 <= min_p2p < threshold_frontal_p2p:
        raise ValueError("Limiares de rejeição inválidos")
    n_epochs = clean_data.shape[0] // ep_len
    
    valid_epochs = []
    valid_indices = []
    
    for i in range(n_epochs):
        ep = clean_data[i*ep_len : (i+1)*ep_len, :]
        p2p_frontal = np.ptp(ep[:, :2], axis=0)
        
        if not np.isfinite(ep).all() or np.any(p2p_frontal > threshold_frontal_p2p) or np.any(p2p_frontal < min_p2p):
            continue
            
        valid_epochs.append(ep)
        valid_indices.append(i)
        
    return np.asarray(valid_epochs).reshape(-1, ep_len, clean_data.shape[1]), np.asarray(valid_indices, dtype=int)
