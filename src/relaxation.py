"""Descritores exploratórios comuns a EEGMAT e ao piloto; não escala de relaxamento."""
from fractions import Fraction
import numpy as np
import pandas as pd
from scipy import signal, stats
from src.preprocessing import segment_epochs

METRICS = ['alpha_rel', 'theta_rel', 'log_alpha_lowbeta', 'faa']
SEED = 20260928


def frontal_features(raw_uv, fs, threshold=150., notch20=False):
    """F3/F4 em µV, filtro por janela, bordas 2 s removidas, épocas de 2 s.

    Retorna índices nominais desde início da janela (primeira época inicia em 2s).
    Não aplica notch a 20 Hz na análise principal.
    """
    raw_uv = np.asarray(raw_uv, dtype=float)
    if raw_uv.ndim != 2 or raw_uv.shape[1] != 2 or not np.isfinite(raw_uv).all():
        raise ValueError('Esperados F3/F4 finitos, em µV, com forma tempo x 2')
    if fs <= 60 or len(raw_uv)/fs < 10:
        raise ValueError('Taxa de amostragem ou duração insuficiente')
    d = signal.detrend(raw_uv, axis=0)
    if notch20:
        b,a = signal.iirnotch(20,30,fs=fs)
        d = signal.filtfilt(b,a,d,axis=0)
    sos = signal.butter(4,[1,30],btype='bandpass',fs=fs,output='sos')
    d = signal.sosfiltfilt(sos,d,axis=0)
    ratio = Fraction(200/float(fs)).limit_denominator(1000)
    d = signal.resample_poly(d,ratio.numerator,ratio.denominator,axis=0)
    d = d[400:-400]
    epochs, indices = segment_epochs(d,threshold_frontal_p2p=threshold)
    candidate = len(d)//400
    if len(epochs) < 2:
        return pd.DataFrame(columns=['epoch_index',*METRICS]),dict(candidate=candidate,accepted=len(epochs))
    f, psd = signal.welch(epochs,fs=200,nperseg=400,axis=1)
    def power(low,high):
        mask = (f>=low)&(f<=high)
        return np.trapezoid(psd[:,mask],f[mask],axis=1)
    total,alpha,theta,lowbeta = power(1,30),power(8,13),power(4,8),power(13,18)
    a = alpha/(total+1e-12)
    frame = pd.DataFrame(dict(epoch_index=indices+1,alpha_rel=a.mean(axis=1),
        theta_rel=(theta/(total+1e-12)).mean(axis=1),
        log_alpha_lowbeta=np.log((alpha+1e-12)/(lowbeta+1e-12)).mean(axis=1),
        faa=np.log(a[:,1]+1e-12)-np.log(a[:,0]+1e-12)))
    return frame,dict(candidate=candidate,accepted=len(epochs))


def paired_summary(differences,seed=SEED,n_boot=5000):
    """IC da média de diferenças, reamostrando pessoas (ou blocos explicitados)."""
    d = np.asarray(differences,dtype=float)
    d = d[np.isfinite(d)]
    if len(d)<3:
        raise ValueError('Menos de três pares válidos')
    rng = np.random.default_rng(seed)
    means = rng.choice(d,size=(n_boot,len(d)),replace=True).mean(axis=1)
    p = 1. if np.allclose(d,0) else stats.wilcoxon(d,alternative='two-sided',method='auto').pvalue
    return dict(n=len(d),difference=float(d.mean()),ci_low=float(np.quantile(means,.025)),
                ci_high=float(np.quantile(means,.975)),p=float(p),
                dz=float(d.mean()/d.std(ddof=1)) if d.std(ddof=1)>0 else np.nan,
                positive_fraction=float(np.mean(d>0)))


def local_block_difference(base,condition,metric,block_s=30,n_boot=5000):
    """IC exploratório condicional ao registro: bootstrap de blocos cronológicos.

    Blocos sem ao menos metade das épocas não são representados. Não é IC entre sujeitos.
    """
    groups=[]
    for frame in (base,condition):
        assigned=(frame.epoch_index*2//block_s).astype(int)
        grouped=frame.groupby(assigned)[metric].agg(['mean','count'])
        groups.append(grouped.loc[grouped['count']>=block_s/4,'mean'].to_numpy())
    a,b=groups
    if min(len(a),len(b))<3:
        return dict(base_blocks=len(a),condition_blocks=len(b),difference=np.nan,ci_low=np.nan,ci_high=np.nan)
    rng=np.random.default_rng(SEED)
    draws=rng.choice(b,(n_boot,len(b))).mean(axis=1)-rng.choice(a,(n_boot,len(a))).mean(axis=1)
    return dict(base_blocks=len(a),condition_blocks=len(b),difference=float(b.mean()-a.mean()),
                ci_low=float(np.quantile(draws,.025)),ci_high=float(np.quantile(draws,.975)))
