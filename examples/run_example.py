"""Synthetic demonstration; sampling times are illustrative."""
import json
import numpy as np
from be_f2_bootstrap import bootstrap_f2


def main():
    times = [0, 5, 10, 15, 20, 30]
    offset = np.arange(-5.5, 6)[:, None]
    test = np.array([0, 30, 55, 72, 82, 90]) + offset * [0, 1, 1, 0.5, 0.3, 0.1]
    reference = np.array([0, 32, 57, 74, 83, 91]) + offset * [0, 0.8, 0.8, 0.4, 0.2, 0.1]
    result = bootstrap_f2(test, reference, times, framework="ema", num_bootstraps=5000, seed=1)
    print(json.dumps({key: result[key] for key in ("observed_statistic", "ci90", "bootstrap_criterion_met")}, indent=2))


if __name__ == "__main__":
    main()
