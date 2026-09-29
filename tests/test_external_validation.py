"""Checks for cross-dataset features and participant-level comparisons."""

import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd

from src.relaxation import METRICS, frontal_features, local_block_difference, paired_summary
from src.pipelines.run_external_validation import grouped_predictions, paired_contrasts
from src.pipelines.run_critical_reanalysis import (
    acquisition_timeline,
    exact_spearman_permutation,
    local_quality_sensitivity,
    paired_signflip_mean,
    sudre_eeg_by_electrodes,
)


def synthetic(fs, alpha=10, beta=3):
    t = np.arange(int(fs * 60)) / fs
    return np.column_stack([
        alpha*np.sin(2*np.pi*10*t) + beta*np.sin(2*np.pi*16*t),
        alpha*np.sin(2*np.pi*10*t+.3) + beta*np.sin(2*np.pi*16*t),
    ])


class FrontalFeatureTests(unittest.TestCase):
    def test_sampling_rates_give_comparable_features(self):
        slow, q_slow = frontal_features(synthetic(200), 200)
        fast, q_fast = frontal_features(synthetic(500), 500)
        self.assertEqual(q_slow, {'candidate': 28, 'accepted': 28})
        self.assertEqual(q_fast, q_slow)
        np.testing.assert_allclose(slow[METRICS].median(), fast[METRICS].median(), atol=.015)
        self.assertEqual(slow.epoch_index.iloc[0], 1)

    def test_relative_features_ignore_uniform_gain(self):
        a, _ = frontal_features(synthetic(200), 200)
        b, _ = frontal_features(synthetic(200)*2, 200)
        np.testing.assert_allclose(a[METRICS], b[METRICS], atol=1e-8)

    def test_more_alpha_increases_alpha_to_low_beta(self):
        low, _ = frontal_features(synthetic(200, alpha=2, beta=8), 200)
        high, _ = frontal_features(synthetic(200, alpha=8, beta=2), 200)
        self.assertGreater(high.log_alpha_lowbeta.median(), low.log_alpha_lowbeta.median())

    def test_rejected_epoch_does_not_change_later_acquisition_index(self):
        raw = synthetic(200)
        raw[4000:4400] += 500
        frame, quality = frontal_features(raw, 200)
        self.assertLess(quality['accepted'], quality['candidate'])
        self.assertGreater(frame.epoch_index.max(), quality['accepted'])

    def test_flatline_and_nonfinite_recordings(self):
        with self.assertRaises(ValueError):
            frontal_features(np.full((12000, 2), np.nan), 200)
        frame, quality = frontal_features(np.zeros((12000, 2)), 200)
        self.assertTrue(frame.empty)
        self.assertEqual(quality['accepted'], 0)


class ParticipantInferenceTests(unittest.TestCase):
    def test_exact_correlation_counts_permutations_of_people(self):
        result = exact_spearman_permutation([1, 2, 3, 4], [1, 2, 3, 4])
        self.assertAlmostEqual(result['rho'], 1)
        self.assertAlmostEqual(result['p_exact'], 2/24)
        with self.assertRaises(ValueError):
            exact_spearman_permutation([1, 1, 1], [1, 2, 3])

    def test_signflip_tests_the_paired_mean_and_is_deterministic(self):
        self.assertEqual(paired_signflip_mean([1, 1, 1])['p_signflip'], .25)
        self.assertEqual(paired_signflip_mean([-1, 0, 1])['p_signflip'], 1)
        sample = np.r_[np.ones(9), np.full(9, -1)]
        result = paired_signflip_mean(sample, draws=999)
        self.assertEqual(result, paired_signflip_mean(sample, draws=999))
        self.assertEqual(result['method'], 'monte_carlo_999')
        with self.assertRaises(ValueError):
            paired_signflip_mean([1, 2])

    def test_paired_contrast_uses_subject_ids_not_row_order(self):
        frame = pd.DataFrame({
            'subject': [1, 2, 3, 3, 1, 2],
            'condition': ['a', 'a', 'a', 'b', 'b', 'b'],
            'x': [100, 200, 300, 301, 101, 201],
        })
        result = paired_contrasts(frame, 'condition', 'b', 'a', ['x']).iloc[0]
        self.assertEqual(result['n'], 3)
        self.assertEqual(result.difference, 1)
        self.assertEqual(result.ci_low, 1)
        self.assertEqual(result.ci_high, 1)

    def test_missing_subject_pair_is_not_counted(self):
        frame = pd.DataFrame({
            'subject': [1, 2, 3, 4, 1, 2, 3],
            'condition': ['a']*4 + ['b']*3,
            'x': [1, 2, 3, 1000, 2, 3, 4],
        })
        result = paired_contrasts(frame, 'condition', 'b', 'a', ['x']).iloc[0]
        self.assertEqual(result['n'], 3)
        self.assertEqual(result.difference, 1)

    def test_duplicate_subject_condition_is_rejected(self):
        frame = pd.DataFrame({
            'subject': [1, 1, 1, 2, 2, 3, 3],
            'condition': ['a', 'a', 'b', 'a', 'b', 'a', 'b'],
            'x': [1, 1, 2, 1, 2, 1, 2],
        })
        with self.assertRaises(ValueError):
            paired_contrasts(frame, 'condition', 'b', 'a', ['x'])

    def test_paired_summary_requires_people_and_seed_is_reproducible(self):
        result = paired_summary([1, 2, 3], n_boot=200)
        self.assertEqual(result, paired_summary([1, 2, 3], n_boot=200))
        self.assertEqual(result['n'], 3)
        with self.assertRaises(ValueError):
            paired_summary([1, 2])

    def test_grouped_predictions_leave_each_person_out(self):
        groups = np.repeat(np.arange(8), 2)
        y = np.tile([0, 1], 8)
        x = np.column_stack([y*2-1, np.zeros(len(y))])
        predictions = grouped_predictions(x, y, groups)
        np.testing.assert_array_equal(predictions >= .5, y)
        np.testing.assert_array_equal(predictions, grouped_predictions(x, y, groups))


class SecondarySensitivityTests(unittest.TestCase):
    def test_timeline_uses_recorded_edges_not_filename_clock(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            header = ','.join(f'col{i}' for i in range(13)) + ',Timestamp\n'
            for stamp, times in [('a', [10, 10, 11]), ('b', [16, 17])]:
                body = ''.join(','.join(['0']*13 + [str(value)])+'\n' for value in times)
                (root / f'{stamp}.csv').write_text(header + body)
            with patch('src.pipelines.run_critical_reanalysis.session_path',
                       side_effect=lambda stamp: root / f'{stamp}.csv'):
                timeline = acquisition_timeline([('A', 'a'), ('B', 'b')])
            self.assertEqual(timeline.gap_from_previous_recording_s.iloc[1], 5)
            self.assertEqual(timeline.repeated_timestamps.iloc[0], 1)
            self.assertEqual(timeline.recorded_duration_s.iloc[0], 1)

    def test_local_blocks_keep_rejection_gaps(self):
        ids = np.r_[np.arange(15), np.arange(30, 60)]
        base = pd.DataFrame({'epoch_index': ids, 'x': np.ones(len(ids))})
        condition = base.copy()
        condition['x'] = 2
        result = local_block_difference(base, condition, 'x')
        self.assertEqual(result['base_blocks'], 3)
        self.assertEqual(result['difference'], 1)
        self.assertEqual(result['ci_low'], 1)

    def test_quality_cutoff_can_remove_a_condition_without_inventing_blocks(self):
        rows = []
        for condition in ['Base', 'BB 1', 'BB 2', 'BB 3', 'BB 4', 'Pós 1', 'Pós 2']:
            count = 5 if condition == 'BB 1' else 12
            for epoch_index in range(count):
                rows.append(dict(condition=condition, epoch_index=epoch_index,
                                 alpha_rel=.1, threshold_uv=150, notch20=False))
        summary, blocks = local_quality_sensitivity(pd.DataFrame(rows), cutoffs=(.3, .8))
        self.assertEqual(len(blocks), 7)
        low = summary[(summary.condition == 'BB 1') & (summary.minimum_retention == .3)].iloc[0]
        high = summary[(summary.condition == 'BB 1') & (summary.minimum_retention == .8)].iloc[0]
        self.assertEqual(low.retained_blocks, 1)
        self.assertEqual(high.retained_blocks, 0)
        self.assertTrue(np.isnan(high.alpha_mean_pct))

    def test_sudre_sensor_choice_changes_coverage_without_mixing_people(self):
        rows = []
        for subject in (1, 2):
            for condition in ('C0', 'C1'):
                for phase in ('Bruit_1', 'Stim1', 'Stim2', 'Stim3'):
                    for electrode in ('F3', 'Fz', 'F4'):
                        if subject == 2 and condition == 'C0' and electrode == 'F3':
                            continue
                        value = .1 if phase == 'Bruit_1' else .2
                        if condition == 'C1' and phase != 'Bruit_1':
                            value += .05
                        rows.append(dict(subject=subject, conditions=condition,
                                         epoch=phase+'_1', electrode=electrode,
                                         relativeAlpha=value))
        frame = pd.DataFrame(rows)
        strict = sudre_eeg_by_electrodes(frame, ('F3', 'Fz', 'F4'))
        fz = sudre_eeg_by_electrodes(frame, ('Fz',))
        self.assertEqual(len(strict), 1)
        self.assertEqual(len(fz), 2)
        np.testing.assert_allclose(fz, [.05, .05])
        duplicate = pd.concat([frame, frame.iloc[[0]].assign(relativeAlpha=.9)])
        with self.assertRaises(ValueError):
            sudre_eeg_by_electrodes(duplicate, ('Fz',))


if __name__ == '__main__':
    unittest.main()
