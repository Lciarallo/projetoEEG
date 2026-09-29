"""
Índice de defasagem ponderado com menor sensibilidade a componentes de fase zero.
O teste exploratório usa inversões de sinal, não randomização de fases de Fourier.
"""

import numpy as np
from scipy import signal

def compute_wpli_epoch_based(epochs, ch1_idx=0, ch2_idx=1, fs=200, band=(8.0, 13.0)):
    """
    Calcula o Weighted Phase Lag Index (wPLI) inter-épocas entre dois canais (Vinck et al., 2011).
    Variante baseada na média imaginária de Hilbert por época; não é wPLI espectral debiased.
    """
    if epochs.ndim != 3 or len(epochs) < 2 or not np.isfinite(epochs).all():
        raise ValueError("wPLI requer pelo menos duas épocas finitas")
    n_ep, ep_len, n_ch = epochs.shape
    sos = signal.butter(4, list(band), btype='bandpass', fs=fs, output='sos')
    
    im_diffs = np.empty(n_ep)
    for i in range(n_ep):
        ep_filt = signal.sosfiltfilt(sos, epochs[i, :, :], axis=0)
        analytic = signal.hilbert(ep_filt, axis=0)
        cross = analytic[:, ch1_idx] * np.conj(analytic[:, ch2_idx])
        im_diffs[i] = np.mean(np.imag(cross))
        
    wpli_obs = np.abs(np.mean(im_diffs)) / (np.mean(np.abs(im_diffs)) + 1e-12)
    return wpli_obs, im_diffs

def wpli_surrogate_test(im_diffs, n_surrogates=1000, seed=42):
    """
    Gera uma distribuição nula não-paramétrica invertendo aleatoriamente os sinais imaginários
    (destrói o acoplamento de fase consistente preservando a amplitude).
    Retorna: p-valor empírico, Z-score e percentil 95 da distribuição nula.
    """
    im_diffs = np.asarray(im_diffs, dtype=float)
    if im_diffs.ndim != 1 or len(im_diffs) < 2 or not np.isfinite(im_diffs).all() or n_surrogates < 1:
        raise ValueError("Amostras ou número de substitutos inválidos")
    rng = np.random.RandomState(seed)
    n_ep = len(im_diffs)
    wpli_real = np.abs(np.mean(im_diffs)) / (np.mean(np.abs(im_diffs)) + 1e-12)
    
    surr_wplis = np.empty(n_surrogates)
    for s in range(n_surrogates):
        signs = rng.choice([-1.0, 1.0], size=n_ep)
        im_surr = im_diffs * signs
        surr_wplis[s] = np.abs(np.mean(im_surr)) / (np.mean(np.abs(im_surr)) + 1e-12)
        
    p_value = (np.sum(surr_wplis >= wpli_real) + 1.0) / (n_surrogates + 1.0)
    z_score = (wpli_real - np.mean(surr_wplis)) / (np.std(surr_wplis) + 1e-12)
    null_p95 = np.percentile(surr_wplis, 95.0)
    
    return p_value, z_score, null_p95, surr_wplis
