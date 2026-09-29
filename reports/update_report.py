"""Build the consolidated report from audited CSV/JSON outputs.

The template holds interpretation; numeric tables come from reproducible runs.
"""
from pathlib import Path
import json
import re

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
REPORTS = ROOT / 'reports'
AUDIT = ROOT / 'results' / 'audit_2026_09_26'
VALIDATION = ROOT / 'results' / 'validation_2026_09_28'
CRITICAL = ROOT / 'results' / 'critical_2026_09_29'
ORDER = ['Base', 'BB 1', 'BB 2', 'BB 3', 'BB 4', 'Pós 1', 'Pós 2']


def table(headers, rows):
    if any(len(row) != len(headers) for row in rows):
        raise ValueError('Tabela com número incorreto de colunas')
    return '\n'.join(
        ['| ' + ' | '.join(headers) + ' |',
         '| ' + ' | '.join(['---'] * len(headers)) + ' |']
        + ['| ' + ' | '.join(map(str, row)) + ' |' for row in rows]
    )


def require_rows(frame, expected, label):
    if len(frame) != expected:
        raise ValueError(f'{label}: esperadas {expected} linhas; encontradas {len(frame)}')
    return frame


def main():
    timeline = require_rows(pd.read_csv(CRITICAL / 'acquisition_timeline.csv'), 4, 'cronologia')
    if timeline.condition.tolist() != ['Base', 'BB', 'Pós 1', 'Pós 2']:
        raise ValueError('Cronologia em ordem inesperada')
    gap = float(timeline.loc[timeline.condition == 'Pós 1', 'gap_from_previous_recording_s'].iloc[0])
    if not 600 < gap < 700:
        raise ValueError('Lacuna BB-Pós 1 fora da faixa verificada')
    acquisition = table(
        ['Arquivo / condição', 'Amostras', 'Duração gravada', 'Lacuna anterior'],
        [(row.condition, int(row.samples), f'{row.recorded_duration_s:.1f} s',
          '-' if pd.isna(row.gap_from_previous_recording_s)
          else f'{row.gap_from_previous_recording_s:.1f} s')
         for row in timeline.itertuples()]
    )

    local = pd.read_csv(VALIDATION / 'local_features.csv')
    local = local[(local.threshold_uv == 150) & (~local.notch20)].set_index('condition').loc[ORDER]
    require_rows(local, 7, 'janelas locais')
    if int(local.accepted.sum()) != len(pd.read_csv(CRITICAL / 'local_primary_epochs.csv')):
        raise ValueError('Contagens locais divergem da reanálise dos arquivos brutos')
    local_table = table(
        ['Condição', 'Épocas válidas / 148', 'Alfa relativo', 'Theta relativo', 'ln(alfa/beta baixa)', 'FAA'],
        [(label, f'{int(row.accepted)}/148', f'{100*row.alpha_rel:.2f}%',
          f'{100*row.theta_rel:.2f}%', f'{row.log_alpha_lowbeta:.3f}', f'{row.faa:.3f}')
         for label, row in local.iterrows()]
    )

    quality = pd.read_csv(CRITICAL / 'local_quality_sensitivity.csv')
    require_rows(quality, 14, 'sensibilidade local')
    quality_rows = []
    for label in ORDER:
        low = quality[(quality.condition == label) & (quality.minimum_retention == .5)]
        high = quality[(quality.condition == label) & (quality.minimum_retention == .8)]
        if len(low) != 1 or len(high) != 1:
            raise ValueError(f'Corte de qualidade ausente: {label}')
        a, b = low.iloc[0], high.iloc[0]
        quality_rows.append((label, int(a.retained_blocks), f'{a.delta_vs_base_pp:+.2f}',
                             int(b.retained_blocks), f'{b.delta_vs_base_pp:+.2f}'))
    quality_table = table(
        ['Condição', 'Blocos ≥50%', 'Δ alfa ≥50% (p.p.)', 'Blocos ≥80%', 'Δ alfa ≥80% (p.p.)'],
        quality_rows
    )

    tests = pd.read_csv(CRITICAL / 'participant_signflip.csv')
    require_rows(tests, 14, 'testes pareados')
    eegmat = tests[tests.family == 'EEGMAT rest-task'].set_index('metric')
    eegmat_table = table(
        ['Indicador', 'Pessoas', 'Repouso - tarefa, média', 'p Holm'],
        [(label, int(eegmat.loc[metric, 'n']),
          f"{eegmat.loc[metric, 'mean_difference'] * scale:+.2f}{unit}",
          f"{eegmat.loc[metric, 'p_holm_within_family']:.5f}")
         for metric, label, scale, unit in [
             ('alpha_rel', 'Alfa relativo', 100, ' p.p.'),
             ('theta_rel', 'Theta relativo', 100, ' p.p.'),
             ('log_alpha_lowbeta', 'ln(alfa/beta baixa)', 1, ''),
             ('faa', 'FAA relativa', 1, ''),
         ]]
    )

    windows = require_rows(pd.read_csv(CRITICAL / 'eegmat_window_summary.csv'), 2, 'janelas EEGMAT')
    window_table = table(
        ['Repouso', 'Pessoas', 'Alfa repouso - tarefa', 'Direção positiva'],
        [('Primeiros 60 s' if row.rest_window == 'rest_first' else 'Últimos 60 s',
          int(row.n), f'{row.mean_difference_pp:+.2f} p.p.',
          f'{round(row.positive_fraction * row.n)}/{int(row.n)}')
         for row in windows.itertuples()]
    )
    model = json.loads((VALIDATION / 'eegmat_model.json').read_text())
    if model['n_subjects'] != 36:
        raise ValueError('Modelo EEGMAT com número inesperado de pessoas')

    questionnaire = tests[tests.family == 'Sudre questionnaire C1-C0']
    require_rows(questionnaire, 6, 'questionários Sudre')
    q_table = table(
        ['Grupo / escala', 'Pares', 'C1 - C0, média', 'p Holm'],
        [(f'{group} / {metric}', int(row.n), f'{row.mean_difference:+.2f}',
          f'{row.p_holm_within_family:.3f}')
         for group in ('G6', 'G40')
         for metric in ('PA', 'CT', 'PT')
         for row in [questionnaire[(questionnaire.group == group) &
                                   (questionnaire.metric == metric)].iloc[0]]]
    )

    sensors = pd.read_csv(CRITICAL / 'sudre_sensor_sensitivity.csv')
    require_rows(sensors, 4, 'sensibilidade de sensores')
    sensor_rows = []
    for group in ('G6', 'G40'):
        for name, family in [('F3/Fz/F4', 'Sudre EEG F3/Fz/F4 C1-C0'),
                             ('Fz', 'Sudre EEG Fz C1-C0')]:
            estimate = sensors[(sensors.group == group) & (sensors.sensors == name)]
            test = tests[(tests.group == group) & (tests.family == family)]
            if len(estimate) != 1 or len(test) != 1:
                raise ValueError(f'Par EEG ausente: {group}, {name}')
            sensor_rows.append((group, name, int(estimate.iloc[0].n),
                                f'{estimate.iloc[0].difference:+.4f}',
                                f'{test.iloc[0].p_holm_within_family:.3f}'))
    sensor_table = table(
        ['Grupo', 'Canais', 'Pares', 'Delta alfa C1 - C0', 'p Holm'],
        sensor_rows
    )

    correlations = require_rows(pd.read_csv(CRITICAL / 'sudre_exact_correlations.csv'), 2,
                                'correlações exatas').set_index('group')
    order = pd.read_csv(CRITICAL / 'sudre_order_sensitivity.csv')
    order = order[order.metric == 'CT']
    require_rows(order, 4, 'sensibilidade de ordem CT')
    order_table = table(
        ['Grupo', 'Posição de C1', 'Pessoas', 'CT C1 - C0, média'],
        [(group, 'Antes de C0' if before else 'Depois de C0',
          int(row.n), f'{row.mean_difference:+.2f}')
         for group in ('G6', 'G40')
         for before in (True, False)
         for row in [order[(order.group == group) &
                           (order.order == ('C1 antes de C0' if before else 'C1 depois de C0'))].iloc[0]]]
    )

    connectivity = require_rows(pd.read_csv(AUDIT / 'eletrofisiologia.csv'), 7, 'defasagem local')
    connectivity = connectivity.set_index('condicao').loc[ORDER]
    connectivity_table = table(
        ['Condição', 'Índice de defasagem', 'p Holm (7)'],
        [(label, f'{row.wpli:.3f}', f'{row.p_holm_7:.4f}')
         for label, row in connectivity.iterrows()]
    )
    classifiers = []
    for name, filename in [('EEGNet', 'eegnet_temporal.csv'),
                           ('Random Forest', 'random_forest_temporal.csv')]:
        frame = require_rows(pd.read_csv(AUDIT / filename), 6, name)
        classifiers.append((name, f'{100*frame.balanced_accuracy.mean():.1f}% ± '
                                  f'{100*frame.balanced_accuracy.std():.1f} p.p.',
                            f'{100*frame.null_balanced_accuracy.mean():.1f}%',
                            f'{frame.auc.mean():.3f}'))
    classifier_table = table(
        ['Modelo', 'Acurácia balanceada (média ± DP)', 'Controle diagnóstico', 'ROC-AUC'],
        classifiers
    )

    values = {
        'ACQUISITION_TABLE': acquisition,
        'GAP_BB_POST_TEXT': f'{gap:.0f} s ({int(gap // 60)} min {round(gap % 60):02d} s)',
        'LOCAL_TABLE': local_table,
        'QUALITY_TABLE': quality_table,
        'EEGMAT_TABLE': eegmat_table,
        'WINDOW_TABLE': window_table,
        'MODEL_BALANCED': f"{100*model['balanced_accuracy']:.2f}",
        'MODEL_AUC': f"{model['auc']:.3f}",
        'MODEL_PERMUTATIONS': str(model['n_permutations']),
        'MODEL_P': f"{model['permutation_p']:.3f}",
        'QUESTIONNAIRE_TABLE': q_table,
        'SENSOR_TABLE': sensor_table,
        'CORR_G6_N': str(int(correlations.loc['G6', 'n'])),
        'CORR_G6_RHO': f"{correlations.loc['G6', 'rho']:.3f}",
        'CORR_G6_P': f"{correlations.loc['G6', 'p_holm_2']:.3f}",
        'CORR_G40_N': str(int(correlations.loc['G40', 'n'])),
        'CORR_G40_RHO': f"{correlations.loc['G40', 'rho']:.3f}",
        'CORR_G40_P': f"{correlations.loc['G40', 'p_holm_2']:.3f}",
        'ORDER_TABLE': order_table,
        'CONNECTIVITY_TABLE': connectivity_table,
        'CONNECTIVITY_MIN_P': f'{connectivity.p_holm_7.min():.4f}',
        'CLASSIFIER_TABLE': classifier_table,
    }
    template = (REPORTS / 'relatorio_cientifico_template.md').read_text(encoding='utf-8')
    placeholders = set(re.findall(r'\{\{([A-Z0-9_]+)\}\}', template))
    if placeholders != set(values):
        raise ValueError(f'Placeholders incorretos: sem valor {placeholders-set(values)}; '
                         f'sem uso {set(values)-placeholders}')
    for key, value in values.items():
        template = template.replace('{{' + key + '}}', value)
    if '{{' in template:
        raise ValueError('Placeholder não resolvido')
    output = REPORTS / 'relatorio_cientifico.md'
    output.write_text(template, encoding='utf-8')
    print(f'Relatório consolidado: {output}')


if __name__ == '__main__':
    main()
