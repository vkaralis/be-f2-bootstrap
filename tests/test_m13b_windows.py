"""M13B-specific boundary, plateau and replicate-window regressions."""
import unittest
import numpy as np
from be_f2_bootstrap import bootstrap_f2
from be_f2_bootstrap.core import _select_window


def profiles(values):
    return np.tile(values, (12, 1)).astype(float)


class M13BWindowTests(unittest.TestCase):
    def window(self, test, reference=None, framework="ich-m13b"):
        return _select_window(test, test if reference is None else reference,
                              np.arange(1, test.shape[1] + 1), framework)

    def test_exactly_85_is_m13b_cutoff_but_not_ema(self):
        test = profiles([10, 40, 70, 85, 92])
        columns, reason = self.window(test)
        self.assertEqual(columns.tolist(), [0, 1, 2, 3])
        self.assertEqual(reason, "85-percent-cutoff")
        self.assertEqual(self.window(test, framework="ema")[0].tolist(), [0, 1, 2, 3, 4])
        result = bootstrap_f2(test, test, [5, 10, 15, 20, 30], framework="ich-m13b", num_bootstraps=5000)
        self.assertEqual(result["selected_columns"], [0, 1, 2, 3])
        self.assertTrue(all(columns == [0, 1, 2, 3] for columns in result["bootstrap_selected_columns"]))

    def test_reference_alone_reaches_85(self):
        test = profiles([10, 40, 65, 75, 82])
        reference = profiles([10, 40, 70, 85, 92])
        self.assertEqual(self.window(test, reference)[0].tolist(), [0, 1, 2, 3])

    def test_plateau_range_exactly_five(self):
        test = profiles([10, 40, 60, 62, 65])
        columns, reason = self.window(test)
        self.assertEqual(columns.tolist(), [0, 1, 2])
        self.assertEqual(reason, "confirmed-sub85-plateau")

    def test_adjacent_differences_do_not_define_plateau(self):
        # Each adjacent change is five, but the three-point range is ten.
        self.assertEqual(self.window(profiles([10, 40, 60, 65, 70]))[0].tolist(), list(range(5)))

    def test_both_products_must_plateau(self):
        test = profiles([10, 40, 60, 62, 64, 64])
        reference = profiles([10, 30, 45, 60, 72, 80])
        self.assertEqual(self.window(test, reference)[0].tolist(), list(range(6)))

    def test_different_plateau_starts_use_common_later_start(self):
        test = profiles([10, 40, 60, 62, 64, 64])
        reference = profiles([10, 30, 45, 60, 62, 64])
        self.assertEqual(self.window(test, reference)[0].tolist(), [0, 1, 2, 3])

    def test_fast_sub85_plateau_keeps_three_points(self):
        test = profiles([10, 60, 62, 64])
        self.assertEqual(self.window(test)[0].tolist(), [0, 1, 2])

    def test_slow_release_false_plateau_is_not_used(self):
        test = profiles([10, 40, 60, 62, 63, 75, 77, 79])
        columns, reason = self.window(test)
        self.assertEqual(columns.tolist(), list(range(6)))
        self.assertEqual(reason, "confirmed-sub85-plateau")

    def test_more_than_six_raw_points_can_select_six_or_fewer(self):
        test = profiles([10, 30, 50, 70, 85, 90, 95, 99])
        result = bootstrap_f2(test, test, np.arange(1, 9), framework="ich-m13b", num_bootstraps=5000)
        self.assertEqual(result["selected_columns"], list(range(5)))

    def test_more_than_six_selected_points_reports_software_limit(self):
        with self.assertRaisesRegex(ValueError, "software restriction, not an absolute guideline prohibition"):
            self.window(profiles([10, 20, 30, 40, 50, 60, 70]))

    def test_each_replicate_reselects_cutoff_using_full_profiles(self):
        test = profiles([10, 40, 70, 84, 94, 99])
        test[:, 3] += np.tile([-6, 6], 6)
        result = bootstrap_f2(test, test, np.arange(1, 7), framework="ich-m13b", num_bootstraps=5000, seed=2)
        self.assertEqual(result["analysis_window_policy"], "per-replicate")
        self.assertEqual(result["selected_columns"], list(range(5)))
        self.assertEqual({len(cols) for cols in result["bootstrap_selected_columns"]}, {4, 5})
        # Reconstruct draws; determine expected windows independently of selector.
        rng = np.random.Generator(np.random.PCG64(2))
        for index in range(100):
            t = test[rng.integers(12, size=12)]
            r = test[rng.integers(12, size=12)]
            count = 4 if t[:, 3].mean() >= 85 or r[:, 3].mean() >= 85 else 5
            self.assertEqual(result["bootstrap_selected_columns"][index], list(range(count)))
            expected = 50 * np.log10(100 / np.sqrt(1 + np.mean((t[:, :count].mean(0) - r[:, :count].mean(0)) ** 2)))
            self.assertAlmostEqual(result["bootstrap_f2"][index], expected)

    def test_each_replicate_reselects_plateau(self):
        test = profiles([10, 40, 65, 69, 70, 70])
        test[:, 2] += np.tile([-6, 6], 6)
        result = bootstrap_f2(test, test, np.arange(1, 7), framework="ich-m13b", num_bootstraps=5000, seed=3)
        self.assertEqual(result["selected_columns"], [0, 1, 2])
        self.assertEqual({len(cols) for cols in result["bootstrap_selected_columns"]}, {3, 4})

    def test_invalid_replicate_is_not_discarded_or_redrawn(self):
        test = profiles([10, 20, 30, 40, 50, 86, 100])
        test[:, 5] += np.tile([-20, 20], 6)
        with self.assertRaisesRegex(ValueError, "M13B bootstrap replicate .*six selected points.*no replicates were discarded or redrawn"):
            bootstrap_f2(test, test, np.arange(1, 8), framework="ich-m13b", num_bootstraps=5000, seed=1)

    def test_cutoff_with_fewer_than_three_points_is_not_extended(self):
        test = profiles([40, 85, 90])
        with self.assertRaisesRegex(ValueError, "fewer than three"):
            self.window(test)

    def test_replicate_early_cutoff_is_reported(self):
        test = profiles([10, 84, 92, 99])
        test[:, 1] += np.tile([-12, 12], 6)
        with self.assertRaisesRegex(ValueError, "M13B bootstrap replicate .*fewer than three.*no replicates were discarded or redrawn"):
            bootstrap_f2(test, test, np.arange(1, 5), framework="ich-m13b", num_bootstraps=5000, seed=1)

    def test_ema_window_remains_fixed(self):
        test = profiles([10, 40, 70, 84, 94, 99])
        test[:, 3] += np.tile([-6, 6], 6)
        result = bootstrap_f2(test, test, np.arange(1, 7), framework="ema", num_bootstraps=5000, seed=2)
        self.assertEqual(result["analysis_window_policy"], "fixed-observed")
        self.assertTrue(all(cols == list(range(5)) for cols in result["bootstrap_selected_columns"]))


if __name__ == "__main__":
    unittest.main()
