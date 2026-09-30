# be-f2-bootstrap

Python package and command-line tool for comparing Test and Reference dissolution
profiles with whole-unit bootstrap resampling and a two-sided 90% percentile CI.
Author: **Vangelis D. Karalis**.

## Install

Requires Python 3.10+ and NumPy. From the repository root:

```sh
python -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[dev]'
python examples/run_example.py
```

On Windows, activate with `.venv\Scripts\activate`. The demonstration uses
synthetic measurements and illustrative sampling times.

## Python API

```python
import numpy as np
from be_f2_bootstrap import bootstrap_f2

test = np.loadtxt('test.csv', delimiter=',')
reference = np.loadtxt('reference.csv', delimiter=',')
times = [0, 5, 10, 15, 20, 30]  # replace with actual protocol times
result = bootstrap_f2(test, reference, times,
                      framework='ema', num_bootstraps=30000, seed=1)
print(result['ci90'], result['bootstrap_criterion_met'])
```

Rows represent complete individual unit trajectories; columns represent shared
sampling times in minutes. At least 12 units per product are required. CSV files
must contain numeric measurements without headers. Missing, negative and
nonfinite values are rejected. Values above 100% are retained without clipping.

| Parameter | Value | Purpose |
| --- | --- | --- |
| `framework` | Required: `'ema'` or `'ich-m13b'` | Explicit estimator and criterion selection |
| `num_bootstraps` | Default 30000; minimum 5000 | Number of bootstrap replicates |
| `seed` | Default 1; integer 0 to 2^32-1 | Local PCG64 generator seed |

EMA mode uses expected f2 with lower CI bound >=50. M13B mode uses conventional
f2 with lower bound >=46 and bootstrap median >=50.
These numerical flags do not determine regulatory eligibility. Read
[methodology](docs/methodology.md) for assumptions and validation scope.

Both modes exclude zero and resample whole unit profiles independently.
EMA uses a fixed observed window through the first mean >85%. M13B uses >=85%
and confirmed sub-85 plateaus, and reselects its window from full data in every
replicate. Plateau detection requires three successive mean values spanning
at most five percentage points and, as an implementation policy, no later
departure from that band. The shared window ends at the later product's plateau
start, retaining at least three points for an early sub-85 plateau.

M13B supports at most six **selected** points, including in each replicate;
more raw sampling columns are accepted if selection reduces them to six or fewer.
This is an implementation restriction. The guideline permits more points with
adequate justification. Invalid replicate windows raise an error with the
replicate number; no replicates are discarded or redrawn. See
[methodology](docs/methodology.md) for the full selection policy.

## Command line

After installation:

```sh
be-f2-bootstrap --test test.csv --reference reference.csv \
  --times 0 5 10 15 20 30 --framework ema \
  --num-bootstraps 30000 --seed 1 --output results
```

Alternatively use `python -m be_f2_bootstrap` with the same arguments. Replace
the illustrative times with the actual sampling schedule. The output directory
contains `analysis.json`, `bootstrap_f2.csv`, both input matrices and `times.json`.
An existing output directory's files with these names are overwritten.

Results include raw statistics, mean/SD/CV profiles, zero-based selected columns,
selected times, software versions, seed and analysis settings.
`selected_columns` describes the observed dataset. `bootstrap_selected_columns`
and `bootstrap_window_reasons` record each replicate's actual analysis window. The CI uses
Hyndman–Fan type 7 quantiles. Threshold comparisons use unrounded values.
Global NumPy random state is left unchanged. Reproduction requires the same
inputs, settings and software versions; MATLAB random sequences are not matched.

## Tests and build

```sh
python -m unittest discover -s tests -v
python -m build
```

GitHub Actions runs tests, the example and package builds across Python
3.10–3.14. Tests cover analytic formulas, Reference indexing, point selection,
whole-unit draws, repeatability, invalid inputs, both criteria and CLI exports.
M13B regressions cover exactly 85%, plateau boundaries and artefacts, dynamic
replicate selection and more than six raw or selected points.

## Files

- `be_f2_bootstrap/`: Python API and CLI.
- `examples/run_example.py`: runnable synthetic demonstration.
- `data/legacy_test.csv`, `data/legacy_reference.csv`: original measurements.
- `tests/`: regression tests.
- `.gitignore`, `.gitattributes`: excluded artifacts and text line endings.
- `.github/workflows/python.yml`: Python CI workflow.
- `pyproject.toml`, `MANIFEST.in`: package metadata and source archive contents.
- `docs/methodology.md`: statistical specification and references.

Actual sampling times were absent from the original dataset and must be supplied.
No open-source license has been granted for this repository.
