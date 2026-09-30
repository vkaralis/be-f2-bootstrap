# Changelog

## 0.3.0 — 2026-09-30

- Use >=85% for M13B cutoff; retain EMA's >85% cutoff.
- Add confirmed sub-85 plateau selection with a documented suffix check to
  reject apparent plateaus followed by further release.
- Reselect M13B windows from full profiles within every bootstrap draw and
  record actual columns and reasons for each draw. EMA remains fixed.
- Report invalid replicate windows without discarding or redrawing replicates.
- Describe six selected M13B points as a software limit, not a guideline ban.
- Add cutoff, plateau, per-replicate selection and >6-point regression cases.
- Verify ignore rules, repository filenames, CI files and upload archive.

## 0.2.0 — 2026-09-30

- Replace MATLAB deliverables with an installable Python package using NumPy.
- Add Python API, CSV command-line interface and JSON/CSV exports.
- Preserve original measurements as numeric CSV files.
- Use an isolated PCG64 generator and record Python and NumPy versions.
- Add executable Python regression tests and Python GitHub Actions builds.
- Keep MATLAB files only as ignored local references.

## 0.1.0 — 2026-09-30

- Initial MATLAB refactor: corrected Reference indexing, whole-unit resampling,
  fixed analysis window, expected-f2 EMA mode and percentile 90% intervals.

Legacy results are not reproduced: the original script contained data indexing,
point selection and interval calculation defects.
