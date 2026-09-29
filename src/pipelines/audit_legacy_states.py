"""Reproduce the legacy FAA state labels from the original OpenBCI files.

The labels are historical software output, not measurements of subjective
relaxation. The eight-second windows use the calculation of the legacy curve,
restricted to each reported condition; the aggregate FAA reproduces its table.
"""

import hashlib
import json

import numpy as np
import pandas as pd
from scipy import signal

from src.study import ROOT, SESSIONS, session_path


OUT = ROOT / 'results' / 'critical_2026_09_29'
LEGACY_TABLE = ROOT / 'results' / 'tables' / 'tabela_assimetria_alfa_frontal.csv'
CONDITIONS = ('Base', 'BB 1', 'BB 2', 'BB 3', 'BB 4', 'Pós 1', 'Pós 2')
POSITIVE = 'Predomínio Esquerdo (Aproximação/Relaxamento)'
NEGATIVE = 'Predomínio Direito (Evitação/Alerta)'


def alpha_power(data):
    frequencies, psd = signal.welch(data, fs=200, nperseg=min(len(data), 400), axis=0)
    alpha = (frequencies >= 8) & (frequencies <= 13)
    integrate = getattr(np, 'trapezoid', getattr(np, 'trapz', None))
    return integrate(psd[alpha], frequencies[alpha], axis=0)


def faa(power):
    if np.any(power[:2] <= 0) or not np.isfinite(power[:2]).all():
        raise ValueError('Potência alfa frontal inválida')
    return float(np.log(power[1]) - np.log(power[0]))


def preprocess_like_legacy(path):
    raw = pd.read_csv(path, skiprows=1, header=None, comment='%', usecols=[1, 2, 3, 4]).to_numpy(float)
    if raw.ndim != 2 or raw.shape[1] != 4 or not np.isfinite(raw).all():
        raise ValueError(f'Dados OpenBCI inválidos: {path}')
    detrended = signal.detrend(raw, axis=0)
    b, a = signal.iirnotch(60.0, 30.0, fs=200)
    notched = signal.filtfilt(b, a, detrended, axis=0)
    sos = signal.butter(4, [1.0, 45.0], btype='bandpass', fs=200, output='sos')
    return signal.sosfiltfilt(sos, notched, axis=0)


def legacy_windows(data):
    """Same eight-second and two-second-step rule, within one condition."""
    values = []
    for start in range(0, len(data) - 1600, 400):
        power = alpha_power(data[start:start + 1600])
        if power[0] > 1e-6 and power[1] > 1e-6:
            values.append(faa(power))
    if not values:
        raise ValueError('Nenhuma janela FAA válida')
    return np.asarray(values)


def audit():
    rows = []
    for label, stamp in SESSIONS:
        clean = preprocess_like_legacy(session_path(stamp))
        sections = ([(f'BB {index + 1}', clean[index * 60000:(index + 1) * 60000])
                     for index in range(4)] if label == 'BB' else [(label, clean)])
        for condition, segment in sections:
            if len(segment) < 1600:
                raise ValueError(f'Segmento incompleto: {condition}')
            power = alpha_power(segment)
            aggregate = faa(power)
            windows = legacy_windows(segment)
            rows.append(dict(condition=condition, alpha_f3=power[0], alpha_f4=power[1],
                             faa_aggregate=aggregate,
                             legacy_state=POSITIVE if aggregate > 0 else NEGATIVE,
                             windows_positive=int(np.count_nonzero(windows > 0)),
                             windows_total=len(windows),
                             windows_positive_pct=100 * np.mean(windows > 0),
                             faa_window_median=float(np.median(windows))))
    result = pd.DataFrame(rows)
    if result.condition.tolist() != list(CONDITIONS):
        raise ValueError('Ordem das condições inesperada')

    historical = pd.read_csv(LEGACY_TABLE)
    if len(historical) != len(result):
        raise ValueError('Tabela histórica incompleta')
    for old, new in zip(historical.itertuples(index=False), result.itertuples(index=False)):
        if old.Estado != new.legacy_state or not np.allclose(
                [old.Alpha_F3, old.Alpha_F4, old.FAA],
                [new.alpha_f3, new.alpha_f4, new.faa_aggregate], rtol=1e-9, atol=1e-9):
            raise ValueError(f'Tabela histórica não reproduzida: {old.Condicao}')
    return result


def main():
    result = audit()
    OUT.mkdir(parents=True, exist_ok=True)
    result.to_csv(OUT / 'legacy_state_audit.csv', index=False)
    inputs = [session_path(stamp) for _, stamp in SESSIONS]
    inputs += [LEGACY_TABLE, ROOT / 'legacy' / 'compute_frontal_alpha_asymmetry.py',
               ROOT / 'src' / 'pipelines' / 'audit_legacy_states.py']
    manifest = dict(parameters=dict(fs=200, alpha_hz=[8, 13], notch_hz=60,
                                    bandpass_hz=[1, 45], window_s=8, step_s=2,
                                    aggregate='log of whole-segment Welch alpha-power ratio',
                                    state_threshold=0),
                    input_sha256={str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                                  for path in inputs})
    (OUT / 'legacy_state_audit_manifest.json').write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
    print(result[['condition', 'faa_aggregate', 'legacy_state', 'windows_positive',
                  'windows_total']].to_string(index=False))


if __name__ == '__main__':
    main()
