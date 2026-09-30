"""Numerical methods; see docs/methodology.md for the statistical contract."""
import platform
from numbers import Integral

import numpy as np


def _statistic(test, reference, framework):
    distance = (test.mean(axis=0) - reference.mean(axis=0)) ** 2
    if framework == "ema":
        distance += test.var(axis=0, ddof=1) / len(test)
        distance += reference.var(axis=0, ddof=1) / len(reference)
    return float(50 * np.log10(100 / np.sqrt(1 + distance.mean())))


def _numeric(value, name, ndim):
    raw = np.asarray(value)
    if raw.dtype.kind not in "iuf" or raw.ndim != ndim or not raw.size:
        raise ValueError(f"{name} must be a nonempty {ndim}-dimensional real numeric array")
    array = raw.astype(np.float64)
    if not np.isfinite(array).all() or (array < 0).any():
        raise ValueError(f"{name} must contain finite nonnegative values")
    return array


def _cv(data):
    mean = data.mean(axis=0)
    sd = data.std(axis=0, ddof=1)
    return np.divide(100 * sd, mean, out=np.zeros_like(mean), where=mean != 0)


def _plateau_start(profile):
    """First confirmed sub-85 plateau; reject later departures from its band.

    Three consecutive means must span <=5 percentage points. As an explicit
    implementation policy, the entire available suffix must also span <=5,
    preventing an early slow-release artefact from truncating the profile.
    The result is a position in the nonzero-time profile, not a column index.
    """
    for start in range(len(profile) - 2):
        suffix = profile[start:]
        if np.all(suffix < 85) and np.ptp(profile[start:start + 3]) <= 5 and np.ptp(suffix) <= 5:
            return start
    return None


def _select_window(test, reference, times, framework):
    """Select shared columns using the specified framework's window policy."""
    nonzero = np.flatnonzero(times > 0)
    mean_t = test[:, nonzero].mean(axis=0)
    mean_r = reference[:, nonzero].mean(axis=0)
    crossing = np.flatnonzero(
        (mean_t >= 85) | (mean_r >= 85) if framework == "ich-m13b"
        else (mean_t > 85) | (mean_r > 85)
    )
    stop = int(crossing[0]) + 1 if crossing.size else len(nonzero)
    reason = "85-percent-cutoff" if crossing.size else "available-schedule"
    if framework == "ich-m13b":
        plateau_t, plateau_r = _plateau_start(mean_t), _plateau_start(mean_r)
        if plateau_t is not None and plateau_r is not None:
            # Both products must have plateaued. Keep the common chronological
            # prefix through the later start, and at least three points for
            # a rapidly established sub-85 plateau (M13B section 2.3.1).
            plateau_stop = max(3, max(plateau_t, plateau_r) + 1)
            if plateau_stop <= stop:
                stop, reason = plateau_stop, "confirmed-sub85-plateau"
    selected = nonzero[:stop]
    if len(selected) < 3:
        raise ValueError("fewer than three usable nonzero time points; review the protocol")
    if framework == "ich-m13b" and len(selected) > 6:
        raise ValueError(
            "M13B implementation supports at most six selected points; this is "
            "a software restriction, not an absolute guideline prohibition"
        )
    return selected, reason


def bootstrap_f2(test, reference, times, *, framework, num_bootstraps=30000, seed=1):
    """Compare complete unit profiles, returning JSON-compatible results.

    Rows are units, columns are shared times in minutes. ``framework`` is
    explicitly ``ema`` (expected f2) or ``ich-m13b`` (conventional f2).
    The 90% interval uses NumPy's linear quantiles (Hyndman-Fan type 7).
    The independent local PCG64 generator leaves global random state untouched.
    M13B reselects its analysis window from full profiles in each replicate;
    EMA holds the observed window fixed. Invalid replicate windows raise an
    error rather than dropping or redrawing replicates.
    """
    if framework not in ("ema", "ich-m13b"):
        raise ValueError("framework must be 'ema' or 'ich-m13b'")
    if isinstance(num_bootstraps, bool) or not isinstance(num_bootstraps, Integral) or num_bootstraps < 5000:
        raise ValueError("num_bootstraps must be an integer >= 5000")
    if isinstance(seed, bool) or not isinstance(seed, Integral) or not 0 <= seed <= 2**32 - 1:
        raise ValueError("seed must be an integer between 0 and 2**32-1")
    test = _numeric(test, "test", 2)
    reference = _numeric(reference, "reference", 2)
    times = _numeric(times, "times", 1)
    if test.shape[1] != len(times) or reference.shape[1] != len(times) or (np.diff(times) <= 0).any():
        raise ValueError("columns must match strictly increasing shared times")
    if len(test) < 12 or len(reference) < 12:
        raise ValueError("at least 12 complete unit profiles per formulation are required")
    selected, window_reason = _select_window(test, reference, times, framework)
    test_used, reference_used = test[:, selected], reference[:, selected]
    sd_t, sd_r = test_used.std(axis=0, ddof=1), reference_used.std(axis=0, ddof=1)
    cv_t, cv_r = _cv(test_used), _cv(reference_used)
    limits = np.where(times[selected] <= 10, 20, 10)
    high_variability = bool(np.any((cv_t > limits) | (cv_r > limits))) if framework == "ema" else bool(np.any((sd_t > 8) | (sd_r > 8)))
    rng = np.random.Generator(np.random.PCG64(int(seed)))
    values = np.empty(num_bootstraps)
    replicate_columns = []
    replicate_reasons = []
    for index in range(num_bootstraps):
        t = test[rng.integers(len(test), size=len(test))]
        r = reference[rng.integers(len(reference), size=len(reference))]
        if framework == "ich-m13b":
            try:
                columns, reason = _select_window(t, r, times, framework)
            except ValueError as exc:
                raise ValueError(f"M13B bootstrap replicate {index + 1}: {exc}; no replicates were discarded or redrawn") from exc
        else:
            columns, reason = selected, window_reason
        replicate_columns.append(columns.tolist())
        replicate_reasons.append(reason)
        values[index] = _statistic(t[:, columns], r[:, columns], framework)
    ci = np.quantile(values, [0.05, 0.95], method="linear")
    median = float(np.median(values))
    bootstrap_pass = ci[0] >= 50 if framework == "ema" else ci[0] >= 46 and median >= 50
    classical = _statistic(test_used, reference_used, "classical")
    return {
        "project_version": "0.3.0", "python_version": platform.python_version(),
        "numpy_version": np.__version__, "framework": framework,
        "estimator": "expected-f2" if framework == "ema" else "classical-f2",
        "num_bootstraps": int(num_bootstraps), "seed": int(seed), "rng": "PCG64",
        "resampling": "whole-unit vectors, independent formulations",
        "percentile_method": "Hyndman-Fan type 7", "confidence_level": 0.90,
        "selected_columns": selected.tolist(), "selected_times": times[selected].tolist(),
        "analysis_window_policy": "per-replicate" if framework == "ich-m13b" else "fixed-observed",
        "observed_window_reason": window_reason,
        "bootstrap_selected_columns": replicate_columns,
        "bootstrap_window_reasons": replicate_reasons,
        "unit_counts": [len(test), len(reference)], "classical_f2": classical,
        "observed_statistic": _statistic(test_used, reference_used, framework),
        "mean_test": test_used.mean(axis=0).tolist(), "mean_reference": reference_used.mean(axis=0).tolist(),
        "sd_test": sd_t.tolist(), "sd_reference": sd_r.tolist(),
        "cv_test": cv_t.tolist(), "cv_reference": cv_r.tolist(),
        "high_variability": high_variability, "bootstrap_f2": values.tolist(),
        "bootstrap_mean": float(values.mean()), "bootstrap_median": median,
        "bootstrap_sd": float(values.std(ddof=1)), "ci90": ci.tolist(),
        "bootstrap_criterion_met": bool(bootstrap_pass),
        "classical_criterion_met": bool(not high_variability and classical >= 50),
    }
