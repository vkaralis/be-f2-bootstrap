# Statistical specification

## Data contract

Measurements are percent dissolved; rows are complete individual unit profiles,
columns are common sampling times in minutes. At least 12 rows per product are
required. Different unit counts are supported, with product-specific variance
terms. NaN, Inf, negative measurements and duplicate times are rejected. No
monotonicity constraint is imposed on experimental measurements.

The original file contained 12 units and six columns, including zero, but no
actual times. `data/legacy_test.csv` and `data/legacy_reference.csv` preserve those values.
They cannot support CV checks based on minutes or rapid-dissolution conclusions
without the actual sampling schedule.

## Window and resampling

Both modes exclude t=0 and use a shared chronological prefix. EMA retains the
first mean >85% in either product and fixes this observed window for all draws.
M13B uses >=85% and evaluates plateaus as described below. Its observed window
is reported separately from the windows recomputed in every bootstrap replicate.

For M13B, a sub-85 plateau candidate has three consecutive mean values with
maximum minus minimum <=5 percentage points. Differences are absolute percentage
points, not relative percentages, and the full three-point range is checked.
Only a plateau in both products can terminate selection. Normally the shared
prefix ends at the later product's first plateau point. An early plateau retains
three points, matching the limited-solubility exception; an early 85% crossing
with fewer than three points still raises an error.

As a conservative **implementation policy**, a candidate is accepted only when
all later available mean values remain below 85% and the full suffix spans <=5
percentage points. This rejects apparent plateaus followed by renewed release.
It is a reproducible operational choice, not a claim that the guideline defines
this exact automatic suffix algorithm. Without a confirmed plateau or cutoff,
all available nonzero points are used. Profile completeness and slow dissolution
still need protocol review; future unmeasured release cannot be inferred.

In each M13B replicate, independently draw whole rows from the complete original
matrices, then apply the same cutoff/plateau selection to the resulting means.
Calculate f2 only on that replicate's selected columns. This permits both shorter
and longer windows than the observed window. The per-replicate policy implements
the requested analysis design; it is not presented as a unique guideline-mandated
bootstrap algorithm. EMA continues to use its fixed observed window.

The implementation requires 3–6 selected M13B points in both the observed dataset
and every replicate. The upper bound is a **software limitation**, not an absolute
M13B prohibition: the Q&A allows justified use of more than six points. More than
six raw columns are accepted when selection meets this limit. No justification
override is provided. A replicate with fewer than three or more than six selected
points aborts with its one-based replicate number. It is never silently omitted,
redrawn, extended or truncated to change the replicate distribution.

`selected_columns` and `selected_times` describe the observed window;
`bootstrap_selected_columns` and `bootstrap_window_reasons` describe every draw.
Column indices are zero-based. `analysis_window_policy` distinguishes
`fixed-observed` from `per-replicate`. Statistics such as reported means, CVs and
`observed_statistic` describe the observed window.

Lag-time adjustment, very-rapid-dissolution exemptions, sub-10% plateau exemptions
and study-specific point selection overrides are not implemented.

## Statistics

Let P be the retained point count, nT and nR the product unit counts, mT and mR
sample mean profiles, and sT² and sR² unbiased sample variances (denominator n−1).

Conventional f2:

```text
50 * log10(100 / sqrt(1 + mean((mT - mR)^2)))
```

Expected f2, used in EMA mode:

```text
50 * log10(100 / sqrt(1 + mean((mT - mR)^2 + sT²/nT + sR²/nR)))
```

The variance terms are added, not subtracted. The same estimator is recalculated
on every replicate. Identical observed mean profiles with nonzero variability
need not have expected f2 = 100. Helpers accept already validated, selected data;
call the public API for analysis.

## Interval and numerical criteria

The CI uses the 5th and 95th percentiles of at least 5000 bootstrap statistics.
Quantiles use Hyndman–Fan type 7: sorted sample x, h=1+(B−1)p, linear interpolation
between floor(h) and ceil(h). This is an interval from the bootstrap distribution,
not a t interval for its mean. Increasing B does not increase the experimental
sample size. Comparisons use unrounded values.

EMA mode: expected-f2 lower bound >=50. M13B mode: conventional-f2 lower bound
>=46 and bootstrap median >=50. Framework
selection is mandatory. The two modes are not interchangeable.

The function additionally reports variability (EMA: CV >20% through 10 min or
>10% later; M13B: SD >8 percentage points), and a separate classical criterion
flag. Bootstrap is always computed. No automatic regulatory pathway is chosen.
A zero mean with zero SD is assigned CV=0.

## Reproducibility and validation status

The function records Python and NumPy versions, project version, seed, PCG64 generator,
replicate count, quantile method, estimator, selected points and unit counts.
Save original inputs and the result together. A local random generator leaves
global NumPy random state unchanged. Selected columns use zero-based indices.
Identical settings reproduce results with the same Python/NumPy environment;
bit-for-bit agreement with MATLAB or across software versions is not asserted.

Tests include independent analytic values, whole-unit resampling reconstruction,
M13B cutoff/plateau and per-replicate window regressions, and CLI export checks. The Python suite is executable
locally and through GitHub Actions. On 2026-09-30, all 24 unittest cases
(including parameterized subcases) passed locally; the wheel and source archive
were built and a synthetic example ran from the installed wheel. Published
reference-dataset validation and study-specific regulatory qualification remain
necessary before relying on this implementation for a submission. The numerical
flags do not establish in-vivo bioequivalence or biowaiver eligibility.

## Sources and scope

The following primary references specify the regulatory methods. Jurisdictional
applicability is not inferred from the numerical criterion alone:

- [EMA Clinical pharmacology and pharmacokinetics Q&A, sections 3.11 and 3.13](https://www.ema.europa.eu/en/human-regulatory-overview/research-development/scientific-guidelines/clinical-pharmacology-pharmacokinetics-guidelines/clinical-pharmacology-pharmacokinetics-questions-answers)
- [ICH M13B Q&A, Step 5](https://www.ema.europa.eu/en/documents/comments/ich-m13b-guideline-bioequivalence-immediate-release-solid-oral-dosage-forms-additional-strengths-biowaiver-questions-answers-step-5_en.pdf)
- [Final ICH M13B guideline, Step 5, sections 2.3.1 and 2.3.2](https://www.ema.europa.eu/system/files/documents/scientific-guideline/ich-m13b-guideline-bioequivalence-immediate-release-solid-oral-dosage-forms-additional-en.pdf)

The final M13B guideline was retrieved for this revision. This software implements
the numerical comparison described here, rather than the complete regulatory
decision tree. Regional implementation dates and application scope must be
assessed separately.
