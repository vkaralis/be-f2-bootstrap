"""Regression tests using Python's standard-library unittest runner."""
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
from be_f2_bootstrap import bootstrap_f2
from be_f2_bootstrap.core import _statistic


class BootstrapTests(unittest.TestCase):
    def setUp(self):
        self.test = np.tile([0, 20, 40, 60, 85, 90], (12, 1))
        self.reference = self.test + [0, 10, 10, 10, 0, 0]
        self.times = [0, 5, 10, 15, 20, 30]

    def analyze(self, test=None, reference=None, times=None, **options):
        defaults = dict(framework="ema", num_bootstraps=5000)
        defaults.update(options)
        return bootstrap_f2(self.test if test is None else test,
                            self.reference if reference is None else reference,
                            self.times if times is None else times, **defaults)

    def test_analytic_statistics(self):
        test = np.tile([20, 40, 60], (12, 1))
        self.assertAlmostEqual(_statistic(test, test + 10, "classical"), 49.891965655433935)
        self.assertEqual(_statistic(test, test, "ema"), 100)
        varying = np.tile(np.arange(1, 13)[:, None], (1, 3))
        expected = 50 * np.log10(100 / np.sqrt(1 + 2 * 13 / 12))
        self.assertAlmostEqual(_statistic(varying, varying, "ema"), expected)

    def test_selection_and_constant_interval(self):
        result = self.analyze()
        self.assertEqual(result["selected_columns"], [1, 2, 3, 4, 5])
        np.testing.assert_allclose(result["ci90"], [result["observed_statistic"]] * 2)

    def test_reference_indexing(self):
        reference = self.test.copy(); reference[:, 4] = 60
        result = self.analyze(reference=reference)
        self.assertAlmostEqual(result["observed_statistic"], _statistic(self.test[:, 1:], reference[:, 1:], "ema"))

    def test_crossing_and_m13b(self):
        test = np.tile([0, 20, 40, 86, 95, 99], (12, 1))
        result = self.analyze(test=test, reference=test, framework="ich-m13b")
        self.assertEqual(result["selected_columns"], [1, 2, 3])
        self.assertEqual(result["ci90"], [100, 100])
        self.assertTrue(result["bootstrap_criterion_met"])

    def test_whole_unit_resampling_seed_quantiles_and_criteria(self):
        for framework in ("ema", "ich-m13b"):
            with self.subTest(framework=framework):
                test = np.tile([20, 40, 60], (12, 1)) + np.arange(1, 13)[:, None]
                state = np.random.get_state()
                result = self.analyze(test, test, [5, 10, 15], framework=framework, seed=2)
                after = np.random.get_state()
                self.assertEqual(state[0], after[0]); np.testing.assert_array_equal(state[1], after[1])
                self.assertEqual(state[2:], after[2:])
                again = self.analyze(test, test, [5, 10, 15], framework=framework, seed=2)
                other = self.analyze(test, test, [5, 10, 15], framework=framework, seed=3)
                self.assertEqual(result, again)
                self.assertNotEqual(result["bootstrap_f2"], other["bootstrap_f2"])
                np.testing.assert_allclose(result["ci90"], np.quantile(result["bootstrap_f2"], [.05, .95], method="linear"))
                rng = np.random.Generator(np.random.PCG64(2))
                t = test[rng.integers(12, size=12)]; r = test[rng.integers(12, size=12)]
                self.assertAlmostEqual(result["bootstrap_f2"][0], _statistic(t, r, framework))
                criterion = result["ci90"][0] >= 50 if framework == "ema" else result["ci90"][0] >= 46 and result["bootstrap_median"] >= 50
                self.assertEqual(result["bootstrap_criterion_met"], criterion)
                json.dumps(result, allow_nan=False)

    def test_invalid_inputs(self):
        for change in ("units", "duplicate", "few_points", "nan", "negative", "complex", "shape", "framework", "count", "seed", "zero_times", "m13b_points"):
            with self.subTest(change=change):
                test, reference, times = self.test.copy(), self.reference.copy(), self.times.copy()
                options = dict(framework="ema", num_bootstraps=5000)
                if change == "units": test = test[:11]
                elif change == "duplicate": times = [0, 5, 5, 15, 20, 30]
                elif change == "few_points": test = np.full((12, 6), 90)
                elif change == "nan": test = test.astype(float); test[0, 1] = np.nan
                elif change == "negative": test[0, 1] = -1
                elif change == "complex": test = test.astype(complex)
                elif change == "shape": times = times[:-1]
                elif change == "framework": options["framework"] = "invalid"
                elif change == "count": options["num_bootstraps"] = 4999
                elif change == "seed": options["seed"] = True
                elif change == "zero_times": test = reference = np.zeros((12, 1)); times = [0]
                elif change == "m13b_points":
                    test = reference = np.tile([10, 20, 30, 40, 50, 60, 70], (12, 1)); times = list(range(1, 8)); options["framework"] = "ich-m13b"
                with self.assertRaises(ValueError):
                    bootstrap_f2(test, reference, times, **options)

    def test_unequal_units_and_zero_cv(self):
        test = np.tile([0, 30, 50], (12, 1))
        reference = np.tile([0, 30, 50], (15, 1))
        result = self.analyze(test, reference, [5, 10, 15])
        self.assertEqual(result["unit_counts"], [12, 15])
        self.assertEqual(result["cv_test"][0], 0)

    def test_cli(self):
        with tempfile.TemporaryDirectory() as directory:
            tmp = Path(directory); output = tmp / "output"
            np.savetxt(tmp / "t.csv", self.test, delimiter=",")
            np.savetxt(tmp / "r.csv", self.reference, delimiter=",")
            completed = subprocess.run([sys.executable, "-m", "be_f2_bootstrap", "--test", str(tmp / "t.csv"),
                "--reference", str(tmp / "r.csv"), "--times", *map(str, self.times), "--framework", "ema",
                "--num-bootstraps", "5000", "--output", str(output)], check=True, capture_output=True, text=True)
            result = json.loads((output / "analysis.json").read_text())
            self.assertEqual(json.loads(completed.stdout)["ci90"], result["ci90"])
            self.assertEqual(len((output / "bootstrap_f2.csv").read_text().splitlines()), 5001)
            np.testing.assert_array_equal(np.loadtxt(output / "test.csv", delimiter=","), self.test)
            self.assertEqual(json.loads((output / "times.json").read_text()), self.times)


if __name__ == "__main__":
    unittest.main()
