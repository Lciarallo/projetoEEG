"""
Módulo de Assimetria Alfa Frontal (FAA) Normalizada por Potência Relativa
com Intervalos de Confiança Não-Paramétricos via Bootstrap.
"""

import numpy as np
from scipy import signal

def compute_relative_faa_epochs(epochs, fs=200, f3_idx=0, f4_idx=1, alpha_band=(8.0, 13.0), tot_band=(1.0, 45.0)):
    """
    Calcula FAA relativa; cancela ganho multiplicativo, mas não corrige todo artefato de contato:
    FAA_rel = ln(RelAlpha_F4) - ln(RelAlpha_F3)
    """
    n_ep = epochs.shape[0]
    faa_vals = []
    rel_alphas_f3 = []
    rel_alphas_f4 = []
    
    trapz = getattr(np, 'trapezoid', getattr(np, 'trapz', None))
    
    for i in range(n_ep):
        ep = epochs[i, :, :]
        f, psd = signal.welch(ep, fs=fs, nperseg=epochs.shape[1], axis=0)
        
        mask_tot = (f >= tot_band[0]) & (f <= tot_band[1])
        mask_a   = (f >= alpha_band[0]) & (f <= alpha_band[1])
        
        tot_pwr = trapz(psd[mask_tot, :], f[mask_tot], axis=0) + 1e-12
        a_pwr   = trapz(psd[mask_a, :], f[mask_a], axis=0)
        
        rel_a = a_pwr / tot_pwr
        
        faa = np.log(rel_a[f4_idx] + 1e-12) - np.log(rel_a[f3_idx] + 1e-12)
        faa_vals.append(faa)
        rel_alphas_f3.append(rel_a[f3_idx] * 100.0)
        rel_alphas_f4.append(rel_a[f4_idx] * 100.0)
        
    return np.array(faa_vals), np.array(rel_alphas_f3), np.array(rel_alphas_f4)

def bootstrap_ci(data, n_boot=2000, ci=95.0, seed=42):
    """Calcula intervalo de confiança percentil não-paramétrico via Bootstrap."""
    data = np.asarray(data, dtype=float)
    if data.ndim != 1 or len(data) < 2 or not np.isfinite(data).all() or n_boot < 1 or not 0 < ci < 100:
        raise ValueError("Bootstrap requer ao menos duas observações finitas e parâmetros válidos")
    rng = np.random.RandomState(seed)
    n = len(data)
    boot_means = np.empty(n_boot)
    for i in range(n_boot):
        sample = rng.choice(data, size=n, replace=True)
        boot_means[i] = np.mean(sample)
    lower = np.percentile(boot_means, (100.0 - ci) / 2.0)
    upper = np.percentile(boot_means, 100.0 - (100.0 - ci) / 2.0)
    return np.mean(data), lower, upper, np.std(data, ddof=1)

def compute_cohens_d(group_exp, group_base):
    """Calcula o tamanho de efeito de Cohen (d) com variância amostral conjunta estrita (ddof=1)."""
    n1, n2 = len(group_exp), len(group_base)
    if min(n1, n2) < 2:
        raise ValueError("Cohen d requer pelo menos duas observações por grupo")
    m1, m2 = np.mean(group_exp), np.mean(group_base)
    s1, s2 = np.std(group_exp, ddof=1), np.std(group_base, ddof=1)
    pooled_std = np.sqrt(((n1 - 1)*(s1**2) + (n2 - 1)*(s2**2)) / (n1 + n2 - 2))
    return (m1 - m2) / pooled_std if pooled_std > 0 else np.nan
