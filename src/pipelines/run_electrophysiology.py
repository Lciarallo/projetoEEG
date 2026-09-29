"""Reanálise reproduzível e sensibilidade a notch/limiar."""
import sys
from pathlib import Path
if __package__ in (None, ''):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import signal
from src.preprocessing import load_raw_openbci, apply_filters, segment_epochs, FS
from src.asymmetry import compute_relative_faa_epochs, bootstrap_ci, compute_cohens_d
from src.connectivity import compute_wpli_epoch_based, wpli_surrogate_test
from src.validation import holm_adjust
from src.study import SESSIONS, RESULTS, session_path, provenance


def main():
    RESULTS.mkdir(parents=True, exist_ok=True)
    quality, rows = [], []
    for label, stamp in SESSIONS:
        path = session_path(stamp)
        raw = load_raw_openbci(path)
        timestamps = pd.read_csv(path, comment='%', usecols=[13]).iloc[:, 0].to_numpy(float)
        delta = np.diff(timestamps)
        quality.append(dict(sessao=label, arquivo=path.name, amostras=len(raw),
                            duracao_nominal_s=len(raw)/FS, duracao_timestamp_s=timestamps[-1]-timestamps[0],
                            timestamps_repetidos=int(np.sum(delta == 0)),
                            timestamps_regressivos=int(np.sum(delta < 0)),
                            saltos_maiores_100ms=int(np.sum(delta > .1)),
                            maior_salto_s=float(delta.max()), epocas_completas=len(raw)//400,
                            amostras_residuais=len(raw)%400))
        for notch in (True, False):
            clean = apply_filters(raw, filter_alias_harmonics=notch)
            blocks = [(f'BB {b+1}', b*60000, (b+1)*60000) for b in range(4)] if label == 'BB' else [(label, 0, 60000)]
            for threshold in (150., 500.):
                for condition, start, stop in blocks:
                    block = clean[start:stop]
                    epochs, indices = segment_epochs(block, threshold_frontal_p2p=threshold)
                    if len(epochs) < 2:
                        raise ValueError(f'Épocas insuficientes: {condition}')
                    faa, _, _ = compute_relative_faa_epochs(epochs)
                    mean, low, high, std = bootstrap_ci(faa, seed=42)
                    wpli, imaginary = compute_wpli_epoch_based(epochs)
                    p, _, _, _ = wpli_surrogate_test(imaginary, seed=42)
                    f, psd = signal.welch(epochs, fs=FS, nperseg=400, axis=1)
                    psd = psd.mean(axis=0)
                    total, beta = (f >= 1) & (f <= 45), (f >= 13) & (f <= 30)
                    integrate = getattr(np, 'trapezoid', None) or np.trapz
                    beta_pct = 100*integrate(psd[beta, 3], f[beta])/integrate(psd[total, 3], f[total])
                    posterior_bad = np.any((np.ptp(epochs[:, :, 2:], axis=1) > threshold) | (np.ptp(epochs[:, :, 2:], axis=1) < 2), axis=1)
                    rows.append(dict(condicao=condition, notch_alias=notch, limiar_uv=threshold,
                                     epocas_total=len(block)//400, epocas_validas=len(epochs),
                                     rejeicao_pct=100*(1-len(epochs)/(len(block)//400)),
                                     posteriores_fora_limiar=int(posterior_bad.sum()),
                                     faa=mean, faa_ci_low=low, faa_ci_high=high, faa_std=std,
                                     wpli=wpli, p_signflip=p, beta_ch4_pct=beta_pct, _faa=faa))
    frame = pd.DataFrame(rows)
    for _, group in frame.groupby(['notch_alias', 'limiar_uv']):
        frame.loc[group.index, 'p_holm_7'] = holm_adjust(group.p_signflip)
        base = group.loc[group.condicao == 'Base', '_faa'].iloc[0]
        frame.loc[group.index, 'cohens_d'] = [compute_cohens_d(a, base) for a in group._faa]
    frame = frame.drop(columns='_faa')
    frame.to_csv(RESULTS / 'sensibilidade_eletrofisiologia.csv', index=False)
    primary = frame[(frame.notch_alias) & (frame.limiar_uv == 150)].copy()
    primary.to_csv(RESULTS / 'eletrofisiologia.csv', index=False)
    pd.DataFrame(quality).to_csv(RESULTS / 'qualidade_aquisicao.csv', index=False)
    fig, axes = plt.subplots(3, 1, figsize=(9, 10), layout='constrained')
    x = np.arange(len(primary))
    axes[0].bar(x, primary.faa, yerr=[primary.faa-primary.faa_ci_low, primary.faa_ci_high-primary.faa], capsize=4)
    axes[0].set_title('FAA relativa: IC 95% bootstrap por época (exploratório)')
    axes[1].bar(x, primary.wpli, color='#b36c25')
    axes[1].set_title('Índice de defasagem por época; p ajustado por Holm')
    for i, (_, row) in enumerate(primary.iterrows()):
        axes[1].text(i, row.wpli+.01, f'p={row.p_holm_7:.3f}', ha='center', fontsize=9)
    axes[1].set_ylim(0, .55)
    for notch, marker in [(True, 'o'), (False, 's')]:
        group = frame[(frame.notch_alias == notch) & (frame.limiar_uv == 150)]
        axes[2].plot(x, group.beta_ch4_pct, marker=marker, label='Com notch 20/80 Hz' if notch else 'Sem notch 20/80 Hz')
    axes[2].set_title('Sensibilidade da potência beta CH4 à filtragem (limiar 150 µV)')
    axes[2].set_ylabel('Potência relativa (%)')
    axes[2].legend()
    for ax in axes:
        ax.set_xticks(x, primary.condicao)
        ax.grid(axis='y', alpha=.2)
    fig.savefig(RESULTS / 'eletrofisiologia.png', dpi=160)
    plt.close(fig)
    provenance('eletrofisiologia', dict(fs=200, epoch_sec=2, seed=42, n_boot=2000,
               n_signflips=1000, primary_threshold_uv=150, sensitivity_threshold_uv=500,
               primary_alias_notch=True, window_seconds=300, rejection_channels=[0,1],
               ci_assumption='iid epochs; exploratory', multiple_testing='Holm, 7 conditions per configuration'))
    print(primary.to_string(index=False))


if __name__ == '__main__':
    main()
