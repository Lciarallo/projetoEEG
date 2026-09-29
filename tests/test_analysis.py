"""Contracts for acquisition order, signal processing and inference."""

import tempfile
import unittest
from pathlib import Path

import numpy as np

from src.asymmetry import bootstrap_ci, compute_relative_faa_epochs
from src.connectivity import compute_wpli_epoch_based, wpli_surrogate_test
from src.features import extract_neuro_features
from src.preprocessing import apply_filters, load_raw_openbci, segment_epochs
from src.validation import holm_adjust, temporal_folds


FS = 200


def tone(frequency, seconds=6, amplitude=10):
    t = np.arange(FS * seconds) / FS
    return amplitude * np.sin(2 * np.pi * frequency * t)


class RawInputTests(unittest.TestCase):
    def test_loader_uses_channel_names_and_keeps_first_sample(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'raw.txt'
            path.write_text(
                '% OpenBCI metadata\n'
                'Sample Index, EXG Channel 2, Other, EXG Channel 0, EXG Channel 3, EXG Channel 1\n'
                '0, 3, 999, 1, 4, 2\n1, 7, 999, 5, 8, 6\n'
            )
            np.testing.assert_array_equal(load_raw_openbci(path), [[1, 2, 3, 4], [5, 6, 7, 8]])

    def test_loader_rejects_missing_channel_and_nonfinite_values(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'raw.txt'
            path.write_text('Index,EXG 0,EXG 1,EXG 2,Wrong\n0,1,2,3,4\n')
            with self.assertRaisesRegex(ValueError, 'EXG 3'):
                load_raw_openbci(path)
            path.write_text('Index,EXG 0,EXG 1,EXG 2,EXG 3\n0,1,NaN,3,4\n')
            with self.assertRaisesRegex(ValueError, 'não finitos'):
                load_raw_openbci(path)


class SignalTests(unittest.TestCase):
    def test_rejection_preserves_original_epoch_ids_and_empty_shape(self):
        raw = np.tile(tone(10, seconds=8)[:, None], (1, 4))
        raw[400:800] = np.nan
        raw[1200:1600, :2] *= 30
        epochs, ids = segment_epochs(raw, threshold_frontal_p2p=150)
        np.testing.assert_array_equal(ids, [0, 2])
        self.assertEqual(epochs.shape, (2, 400, 4))
        empty, rejected_ids = segment_epochs(np.zeros((800, 4)))
        self.assertEqual(empty.shape, (0, 400, 4))
        self.assertEqual(rejected_ids.size, 0)

    def test_rejection_checks_frontal_channels_without_compressing_time(self):
        raw = np.tile(tone(10, seconds=6)[:, None], (1, 4))
        raw[400:800, 1] = 0
        raw[800:1200, 2:] *= 100
        epochs, ids = segment_epochs(raw, threshold_frontal_p2p=150)
        np.testing.assert_array_equal(ids, [0, 2])
        self.assertGreater(np.ptp(epochs[1, :, 2]), 150)

    def test_filter_removes_twenty_hz_but_preserves_alpha(self):
        t = np.arange(10000) / FS
        data = np.tile((np.sin(2*np.pi*20*t) + np.sin(2*np.pi*10*t))[:, None], (1, 4))
        with_notch = apply_filters(data)
        without_notch = apply_filters(data, filter_alias_harmonics=False)
        interior = slice(1000, -1000)
        for frequency, maximum_ratio, minimum_ratio in [(20, .1, 0), (10, 1.2, .8)]:
            basis = np.sin(2*np.pi*frequency*t[interior])
            ratio = abs(with_notch[interior, 0] @ basis) / abs(without_notch[interior, 0] @ basis)
            self.assertLess(ratio, maximum_ratio)
            self.assertGreater(ratio, minimum_ratio)

    def test_features_locate_injected_alpha_and_are_gain_invariant(self):
        rng = np.random.default_rng(17)
        base = np.tile(tone(10, seconds=2)[:, None], (1, 4))
        epochs = np.stack([base + rng.normal(0, .01, base.shape) for _ in range(4)])
        values, names = extract_neuro_features(epochs)
        scaled, _ = extract_neuro_features(epochs * 3)
        self.assertEqual(values.shape, (4, 24))
        self.assertEqual(len(names), 24)
        self.assertTrue(np.isfinite(values).all())
        self.assertGreater(values[:, 8:12].min(), .9)
        np.testing.assert_allclose(values[:, :20], scaled[:, :20], atol=1e-6)

    def test_relative_faa_cancels_channel_gain(self):
        left = tone(10, seconds=2, amplitude=3) + tone(16, seconds=2, amplitude=2)
        right = tone(10, seconds=2, amplitude=6) + tone(16, seconds=2, amplitude=2)
        sample = np.column_stack([left, right, left, right])[None, :, :]
        faa, left_alpha, right_alpha = compute_relative_faa_epochs(sample)
        gain, _, _ = compute_relative_faa_epochs(sample * 4)
        self.assertGreater(faa[0], 0)
        self.assertGreater(right_alpha[0], left_alpha[0])
        np.testing.assert_allclose(faa, gain, atol=1e-8)


class InferenceTests(unittest.TestCase):
    def test_temporal_folds_keep_acquisition_gap_after_rejection(self):
        arrays = [np.delete(np.arange(174), [3, 6, 42, 100]),
                  np.delete(np.arange(612), np.arange(50, 90))]
        folds = list(temporal_folds(arrays))
        self.assertGreaterEqual(len(folds), 3)
        self.assertGreater(folds[-1][2][1][0] * 2, 1000)
        for train, test, windows in folds:
            for ids, tr, te, (start, stop) in zip(arrays, train, test, windows):
                self.assertTrue(np.all((ids[tr]+1 <= start-60) | (ids[tr] >= stop+60)))
                self.assertTrue(np.all((ids[te] >= start) & (ids[te] < stop)))
                self.assertFalse(set(tr) & set(te))

    def test_temporal_folds_reject_unsorted_and_duplicate_indices(self):
        for invalid in ([0, 2, 1], [0, 1, 1], []):
            with self.subTest(invalid=invalid), self.assertRaises(ValueError):
                list(temporal_folds([invalid, np.arange(200)]))

    def test_holm_matches_hand_calculation_and_permutation(self):
        p = np.array([.026, .2, .013, .021, .56, .309, .78])
        expected = np.array([.13, .8, .091, .126, 1, .927, 1])
        np.testing.assert_allclose(holm_adjust(p), expected)
        order = np.array([6, 4, 2, 0, 5, 1, 3])
        np.testing.assert_allclose(holm_adjust(p[order]), expected[order])
        for invalid in ([np.nan], [-.1], [1.1]):
            with self.assertRaises(ValueError):
                holm_adjust(invalid)

    def test_bootstrap_and_signflip_are_seeded_and_two_sided(self):
        data = np.linspace(-2, 1, 20)
        self.assertEqual(bootstrap_ci(data), bootstrap_ci(data))
        first = wpli_surrogate_test(data, n_surrogates=200)
        second = wpli_surrogate_test(data, n_surrogates=200)
        np.testing.assert_array_equal(first[3], second[3])
        self.assertGreater(first[0], 0)
        self.assertLess(wpli_surrogate_test(np.ones(40), n_surrogates=500)[0], .01)
        with self.assertRaises(ValueError):
            bootstrap_ci([])
        with self.assertRaises(ValueError):
            compute_wpli_epoch_based(np.empty((0, 400, 4)))


class ModelTests(unittest.TestCase):
    def test_eegnet_supported_lengths_forward_and_backward(self):
        import torch
        from src.models.eegnet import EEGNet

        torch.set_num_threads(2)
        for length in (400, 500, 512):
            with self.subTest(length=length):
                model = EEGNet(samples=length)
                output = model(torch.randn(2, 1, 4, length))
                self.assertEqual(tuple(output.shape), (2, 2))
                output.sum().backward()
                self.assertIsNotNone(model.conv1.weight.grad)


if __name__ == '__main__':
    unittest.main()
