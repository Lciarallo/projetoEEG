"""Secondary, exploratory checks of timing, selection and external validity.

These analyses describe the available recordings. They do not estimate a causal
effect of binaural audio or a population effect from the local participant.
"""

import ast
from datetime import datetime, timezone
import hashlib
from itertools import permutations
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from src.relaxation import METRICS, SEED, frontal_features, paired_summary
from src.preprocessing import load_raw_openbci
from src.study import ROOT, SESSIONS, session_path
from src.validation import holm_adjust


OUT = ROOT / 'results' / 'critical_2026_09_29'
EXTERNAL = ROOT / 'data' / 'external'
PHASES = ('Bruit_1', 'Stim1', 'Stim2', 'Stim3')


def acquisition_timeline(sessions=SESSIONS):
    """Read recorded Unix timestamps; gaps do not establish audio offset."""
    rows = []
    previous_end = None
    for label, stamp in sessions:
        path = session_path(stamp)
        timestamps = pd.read_csv(path, usecols=[13]).iloc[:, 0].to_numpy(float)
        if len(timestamps) < 2 or not np.isfinite(timestamps).all() or np.any(np.diff(timestamps) < 0):
            raise ValueError(f'Timestamps inválidos: {path}')
        start, end = float(timestamps[0]), float(timestamps[-1])
        if previous_end is not None and start <= previous_end:
            raise ValueError(f'Sessões sobrepostas ou fora de ordem: {path}')
        rows.append(dict(condition=label, file=path.name, samples=len(timestamps),
                         start_unix_s=start, end_unix_s=end,
                         recorded_duration_s=end-start,
                         gap_from_previous_recording_s=np.nan if previous_end is None else start-previous_end,
                         repeated_timestamps=int(np.count_nonzero(np.diff(timestamps) == 0))))
        previous_end = end
    return pd.DataFrame(rows)


def local_primary_epochs():
    """Recompute the primary local epochs directly from the four raw files."""
    frames = []
    for label, stamp in SESSIONS:
        raw = load_raw_openbci(session_path(stamp))[:, :2]
        windows = [(f'BB {index+1}', index*60000, (index+1)*60000)
                   for index in range(4)] if label == 'BB' else [(label, 0, 60000)]
        for condition, start, stop in windows:
            if len(raw[start:stop]) != 60000:
                raise ValueError(f'Janela incompleta: {condition}')
            frame, _ = frontal_features(raw[start:stop], 200, threshold=150, notch20=False)
            frame['condition'] = condition
            frame['threshold_uv'] = 150
            frame['notch20'] = False
            frames.append(frame)
    return pd.concat(frames, ignore_index=True)


def local_quality_sensitivity(epochs, cutoffs=(0.5, 0.8)):
    """Mean of 30 s blocks selected by *their own* epoch retention.

    The 2 s edge removal gives 14 candidate epochs in blocks 0/9 and 15 in
    the others. This is a sensitivity analysis, not a matched comparison.
    """
    required = {'condition', 'epoch_index', 'alpha_rel', 'threshold_uv', 'notch20'}
    if not required.issubset(epochs):
        raise ValueError(f'Colunas ausentes: {required-set(epochs)}')
    data = epochs[(epochs.threshold_uv == 150) & (epochs.notch20 == False)].copy()  # noqa: E712
    if data.empty or not data.alpha_rel.notna().all():
        raise ValueError('Épocas primárias ausentes ou não finitas')
    data['block'] = data.epoch_index.mul(2).floordiv(30).astype(int)
    if data.block.min() < 0 or data.block.max() > 9:
        raise ValueError('Índice de época fora da janela de cinco minutos')
    if data.duplicated(['condition', 'epoch_index']).any():
        raise ValueError('Época duplicada')
    blocks = data.groupby(['condition', 'block']).alpha_rel.agg(['mean', 'count']).reset_index()
    blocks['candidate'] = np.where(blocks.block.isin((0, 9)), 14, 15)
    blocks['retention'] = blocks['count'] / blocks['candidate']
    if (blocks['retention'] > 1).any():
        raise ValueError('Mais épocas aceitas que candidatas')
    conditions = ['Base', 'BB 1', 'BB 2', 'BB 3', 'BB 4', 'Pós 1', 'Pós 2']
    if set(data.condition) != set(conditions):
        raise ValueError('Condições locais incompletas')
    rows = []
    for cutoff in cutoffs:
        if not 0 < cutoff <= 1:
            raise ValueError('Retenção deve estar em (0, 1]')
        kept = blocks[blocks.retention >= cutoff]
        base = kept[kept.condition == 'Base']['mean'].mean()
        for condition in conditions:
            values = kept[kept.condition == condition]
            mean = values['mean'].mean()
            rows.append(dict(condition=condition, minimum_retention=cutoff,
                             retained_blocks=len(values), alpha_mean_pct=100*mean,
                             delta_vs_base_pp=100*(mean-base)))
    return pd.DataFrame(rows), blocks


def sudre_eeg_by_electrodes(frame, electrodes):
    """Paired C1-C0 alpha after subtracting each condition's noise baseline."""
    electrodes = tuple(electrodes)
    if not electrodes or len(set(electrodes)) != len(electrodes):
        raise ValueError('Conjunto de eletrodos inválido')
    eeg = frame.drop_duplicates()
    if eeg.duplicated(['subject', 'conditions', 'epoch', 'electrode']).any():
        raise ValueError('Duplicatas conflitantes de EEG')
    roi = eeg[eeg.electrode.isin(electrodes) & eeg.conditions.isin(['C0', 'C1'])].copy()
    roi['phase'] = roi.epoch.str.extract(r'^(Bruit_1|Stim1|Stim2|Stim3)(?:_|$)', expand=False)
    roi = roi.dropna(subset=['phase'])
    if not np.isfinite(roi.relativeAlpha).all():
        raise ValueError('Alfa não finito')
    channels = roi.groupby(['subject', 'conditions', 'phase', 'electrode']).relativeAlpha.mean().reset_index()
    phases = channels.groupby(['subject', 'conditions', 'phase']).relativeAlpha.agg(['mean', 'count']).reset_index()
    phases = phases[phases['count'] == len(electrodes)]
    wide = phases.pivot(index=['subject', 'conditions'], columns='phase', values='mean')
    complete = wide.reindex(columns=PHASES).dropna()
    delta = complete[['Stim1', 'Stim2', 'Stim3']].mean(axis=1) - complete.Bruit_1
    paired = delta.unstack('conditions').reindex(columns=['C0', 'C1']).dropna()
    return paired.C1 - paired.C0


def sudre_sensor_sensitivity():
    rows = []
    for group in ('G6', 'G40'):
        path = EXTERNAL / 'sudre2024' / f'eeg_{group}.csv'
        eeg = pd.read_csv(path, sep=';')
        for label, electrodes in [('F3/Fz/F4', ('F3', 'Fz', 'F4')), ('Fz', ('Fz',))]:
            differences = sudre_eeg_by_electrodes(eeg, electrodes)
            result = paired_summary(differences)
            rows.append(dict(group=group, sensors=label, subject_ids=','.join(map(str, differences.index)),
                             **result))
    return pd.DataFrame(rows)


def sudre_order_sensitivity():
    """Descriptive stratification; order was randomized in the source study."""
    rows = []
    for group in ('G6', 'G40'):
        path = EXTERNAL / 'sudre2024' / f'relax_{group}.csv'
        frame = pd.read_csv(path, sep=';')
        frame.columns = frame.columns.str.strip()
        if frame.groupby('subject').order.nunique().max() != 1:
            raise ValueError('Ordem inconsistente por participante')
        order = frame.drop_duplicates('subject').set_index('subject').order.apply(ast.literal_eval)
        if not order.apply(lambda sequence: len(sequence) == 4 and set(sequence) == {'C0','C1','C2','C3'}).all():
            raise ValueError('Ordem de condições inválida')
        after = order.apply(lambda sequence: sequence.index('C1') > sequence.index('C0'))
        for metric in ('PA', 'CT', 'PT'):
            wide = frame.pivot(index='subject', columns='Condition', values=metric)
            difference = (wide.C1 - wide.C0).dropna()
            for value, label in [(False, 'C1 antes de C0'), (True, 'C1 depois de C0')]:
                subset = difference[after.reindex(difference.index) == value]
                rows.append(dict(group=group, metric=metric, order=label, n=len(subset),
                                 mean_difference=float(subset.mean()),
                                 median_difference=float(subset.median())))
    return pd.DataFrame(rows)


def paired_signflip_mean(differences, draws=99999, seed=SEED):
    """Two-sided paired-label symmetry test for the mean difference.

    Exact through 17 pairs; larger samples use deterministic Monte Carlo.
    Exchangeability remains a scientific assumption, not a software property.
    """
    d = np.asarray(differences, dtype=float)
    if d.ndim != 1 or len(d) < 3 or not np.isfinite(d).all() or draws < 1:
        raise ValueError('Pares ou número de sorteios inválidos')
    observed = abs(d.mean())
    if len(d) <= 17:
        indices = np.arange(1 << len(d), dtype=np.uint32)[:, None]
        bits = (indices >> np.arange(len(d), dtype=np.uint32)) & 1
        null = abs(((2*bits.astype(float)-1)*d).mean(axis=1))
        p = np.mean(null >= observed - 1e-12)
        method = 'exact'
    else:
        rng = np.random.default_rng(seed)
        signs = 2*rng.integers(0, 2, size=(draws, len(d)), dtype=np.int8)-1
        null = abs((signs*d).mean(axis=1))
        p = (1 + np.count_nonzero(null >= observed - 1e-12))/(draws+1)
        method = f'monte_carlo_{draws}'
    return dict(n=len(d), mean_difference=float(d.mean()), p_signflip=float(p), method=method)


def exact_spearman_permutation(x, y):
    """Exact two-sided person-level association test for small paired samples."""
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    if x.ndim != 1 or y.shape != x.shape or not 3 <= len(x) <= 8 or not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError('Spearman exato requer 3 a 8 pares finitos')
    xr, yr = stats.rankdata(x), stats.rankdata(y)
    xr, yr = xr-xr.mean(), yr-yr.mean()
    denominator = np.linalg.norm(xr)*np.linalg.norm(yr)
    if denominator == 0:
        raise ValueError('Correlação indefinida para vetor constante')
    rho = float(xr @ yr / denominator)
    assignments = np.fromiter((xr @ yr[list(order)] / denominator
                               for order in permutations(range(len(y)))), dtype=float)
    return dict(n=len(x), rho=rho, p_exact=float(np.mean(abs(assignments) >= abs(rho)-1e-12)))


def sudre_exact_correlations():
    rows = []
    for group in ('G6', 'G40'):
        eeg = pd.read_csv(EXTERNAL / 'sudre2024' / f'eeg_{group}.csv', sep=';')
        difference = sudre_eeg_by_electrodes(eeg, ('F3', 'Fz', 'F4'))
        questionnaire = pd.read_csv(EXTERNAL / 'sudre2024' / f'relax_{group}.csv', sep=';')
        questionnaire.columns = questionnaire.columns.str.strip()
        ct = questionnaire.pivot(index='subject', columns='Condition', values='CT')
        joint = pd.concat([difference.rename('alpha'), (ct.C1-ct.C0).rename('ct')], axis=1).dropna()
        rows.append(dict(group=group, **exact_spearman_permutation(joint.alpha, joint.ct)))
    result = pd.DataFrame(rows)
    result['p_holm_2'] = holm_adjust(result.p_exact)
    return result


def eegmat_window_sensitivity():
    """Compare the first and last rest minutes with the same task minute."""
    import mne

    rows = []
    for subject in range(36):
        features = {}
        for code in (1, 2):
            path = EXTERNAL / 'eegmat' / f'Subject{subject:02d}_{code}.edf'
            recording = mne.io.read_raw_edf(path, verbose=False)
            fs = recording.info['sfreq']
            samples = round(60*fs)
            if recording.n_times < samples:
                raise ValueError(f'Registro curto: {path}')
            starts = [('rest_first', 0), ('rest_last', recording.n_times-samples)] if code == 1 else [('task_first', 0)]
            for name, start in starts:
                raw = recording.get_data(picks=['EEG F3', 'EEG F4'], start=start, stop=start+samples).T*1e6
                frame, quality = frontal_features(raw, fs)
                medians = frame[METRICS].median() if quality['accepted'] >= 14 else pd.Series(np.nan, index=METRICS)
                for metric in METRICS:
                    features[f'{name}_{metric}'] = medians[metric]
        rows.append(dict(subject=subject, **features))
    subjects = pd.DataFrame(rows)
    summary = []
    for name in ('rest_first', 'rest_last'):
        differences = (subjects[f'{name}_alpha_rel'] - subjects.task_first_alpha_rel).dropna()
        if len(differences) < 3:
            raise ValueError('Pares EEGMAT insuficientes')
        summary.append(dict(rest_window=name, n=len(differences), mean_difference_pp=100*differences.mean(),
                            median_difference_pp=100*differences.median(),
                            positive_fraction=float(np.mean(differences > 0)),
                            p_wilcoxon=float(stats.wilcoxon(differences).pvalue)))
    return subjects, pd.DataFrame(summary)


def participant_signflip_tests(eegmat_subjects):
    rows = []
    for metric in METRICS:
        differences = (eegmat_subjects[f'rest_last_{metric}'] -
                       eegmat_subjects[f'task_first_{metric}']).dropna()
        rows.append(dict(family='EEGMAT rest-task', group='all', metric=metric,
                         **paired_signflip_mean(differences)))
    for group in ('G6', 'G40'):
        questionnaire = pd.read_csv(EXTERNAL / 'sudre2024' / f'relax_{group}.csv', sep=';')
        questionnaire.columns = questionnaire.columns.str.strip()
        if questionnaire.duplicated(['subject', 'Condition']).any():
            raise ValueError('Questionário duplicado')
        for metric in ('PA', 'CT', 'PT'):
            wide = questionnaire.pivot(index='subject', columns='Condition', values=metric)
            differences = (wide.C1 - wide.C0).dropna()
            rows.append(dict(family='Sudre questionnaire C1-C0', group=group, metric=metric,
                             **paired_signflip_mean(differences)))
        eeg = pd.read_csv(EXTERNAL / 'sudre2024' / f'eeg_{group}.csv', sep=';')
        for family, electrodes in [('Sudre EEG F3/Fz/F4 C1-C0', ('F3','Fz','F4')),
                                   ('Sudre EEG Fz C1-C0', ('Fz',))]:
            differences = sudre_eeg_by_electrodes(eeg, electrodes)
            rows.append(dict(family=family, group=group, metric='alpha_delta',
                             **paired_signflip_mean(differences)))
    result = pd.DataFrame(rows)
    result['p_holm_within_family'] = result.groupby('family').p_signflip.transform(holm_adjust)
    return result


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    timeline = acquisition_timeline()
    local = local_primary_epochs()
    local_summary, local_blocks = local_quality_sensitivity(local)
    sensors = sudre_sensor_sensitivity()
    order = sudre_order_sensitivity()
    subjects, windows = eegmat_window_sensitivity()
    windows['p_holm_2_exploratory'] = holm_adjust(windows.p_wilcoxon)
    signflip = participant_signflip_tests(subjects)
    correlations = sudre_exact_correlations()
    outputs = {'acquisition_timeline': timeline, 'local_primary_epochs': local,
               'local_quality_sensitivity': local_summary,
               'local_quality_blocks': local_blocks, 'sudre_sensor_sensitivity': sensors,
               'sudre_order_sensitivity': order, 'eegmat_window_subjects': subjects,
               'eegmat_window_summary': windows, 'participant_signflip': signflip,
               'sudre_exact_correlations': correlations}
    for name, frame in outputs.items():
        frame.to_csv(OUT / f'{name}.csv', index=False)
    inputs = [session_path(stamp) for _, stamp in SESSIONS]
    inputs += sorted((EXTERNAL / 'eegmat').glob('*.edf'))
    inputs += sorted((EXTERNAL / 'sudre2024').glob('*.csv'))
    inputs += [Path(__file__), ROOT / 'src' / 'relaxation.py', ROOT / 'src' / 'preprocessing.py']
    manifest = dict(generated_utc=datetime.now(timezone.utc).isoformat(),
                    interpretation='secondary exploratory sensitivity; no local causal or population inference',
                    parameters=dict(local_retention_cutoffs=[0.5, 0.8], block_s=30,
                                    local_threshold_uv=150, local_notch20=False,
                                    eegmat_rest_windows=['first 60 s', 'last 60 s'],
                                    sudre_sensors=['F3/Fz/F4', 'Fz'],
                                    signflip_draws=99999, signflip_seed=SEED,
                                    signflip_exact_through_n=17),
                    input_sha256={str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                                  for path in inputs})
    (OUT / 'manifest.json').write_text(json.dumps(manifest, indent=2, ensure_ascii=False)+'\n')
    for name in ('acquisition_timeline', 'local_quality_sensitivity', 'sudre_sensor_sensitivity',
                 'eegmat_window_summary', 'participant_signflip', 'sudre_exact_correlations'):
        print(name, '\n', outputs[name].to_string(index=False), '\n', flush=True)


if __name__ == '__main__':
    main()
