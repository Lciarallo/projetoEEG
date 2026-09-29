"""Descritores espectrais e amplitude imaginária (esta última não é wPLI)."""
import numpy as np
from scipy import signal
from src.preprocessing import FS, CHANNELS


def extract_neuro_features(epochs):
    # Entrada: épocas x tempo x canais, sem z-score.
    f, psd = signal.welch(epochs, fs=FS, nperseg=epochs.shape[1], axis=1)
    integrate = getattr(np, 'trapezoid', None) or np.trapz
    total_mask = (f >= 1) & (f <= 45)
    total = integrate(psd[:, total_mask], f[total_mask], axis=1) + 1e-12
    bands = [('Delta',1,4), ('Theta',4,8), ('Alpha',8,13), ('Beta',13,30), ('Gamma',30,45)]
    relative, names = [], []
    for name, low, high in bands:
        mask = (f >= low) & (f <= high)
        relative.append(integrate(psd[:, mask], f[mask], axis=1)/total)
        names.extend(f'{name}_{ch}' for ch in CHANNELS)
    faa = np.log(relative[2][:,1]+1e-12)-np.log(relative[2][:,0]+1e-12)
    ratios = relative[1][:,:2]/(relative[3][:,:2]+1e-12)
    normalized = (epochs-epochs.mean(axis=1,keepdims=True))/(epochs.std(axis=1,keepdims=True)+1e-8)
    sos = signal.butter(4, [8,13], btype='bandpass', fs=FS, output='sos')
    analytic = signal.hilbert(signal.sosfiltfilt(sos, normalized, axis=1), axis=1)
    imaginary = np.abs(np.imag(analytic[:,:,0]*np.conj(analytic[:,:,1]))).mean(axis=1)
    names += ['FAA_rel', 'TBR_F3', 'TBR_F4', 'MeanAbsImagAlpha_normalized']
    return np.column_stack([*relative, faa, ratios, imaginary]), names
