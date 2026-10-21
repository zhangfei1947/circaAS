# circaAS 0.1.0

**circaAS** is the GitHub project name; the Python distribution and command-line tool are named **circadian-as**.

This initial release provides an installable Python package and CLI for exploratory circadian alternative-splicing analysis. The primary engine uses a beta-binomial (BB) model with bias-corrected null bootstrap calibration. An interface for continuous binomial–logit-normal (BLN) sensitivity analysis is reserved but not yet implemented. This release does not provide certified genome-wide false discovery rate (FDR) control.

The code is organized into input validation, coverage diagnostics, model engines, execution, and provenance modules. A future R interface can reuse the same core after the input/output schema stabilizes. Retaining the current numerical implementation avoids introducing optimizer and random-number differences during an immediate migration to R.

## Installation and quick start

Requires Python 3.10+, NumPy, and SciPy. Plotting additionally requires matplotlib.

```bash
git clone https://github.com/zhangfei1947/circaAS.git
cd circaAS
python -m venv .venv
source .venv/bin/activate
python -m pip install .
circadian-as --version
circadian-as engines
```

Run a minimal example without the original research data. These entirely synthetic counts demonstrate software behavior and are not a performance benchmark:

```bash
python examples/make_demo.py --out demo.jsonl
circadian-as run --input demo.jsonl --out demo-run --bootstrap 199 --workers 2
```

Analyze your own raw rMATS counts:

```bash
circadian-as prepare-rmats \
  --counts-dir /path/to/rmats/output \
  --samples /path/to/samples.tsv \
  --bam-list-1 /path/to/all_bams.txt \
  --dataset my_dataset --mode JC --out prepared

circadian-as coverage --input prepared/records.jsonl --out coverage-run

circadian-as run --input prepared/records.jsonl \
  --out analysis-run --bootstrap 199 --workers 4 --seed 20261020

# Use the event JSON path listed in the artifact column of results.tsv.
python -m pip install '.[plot]'
circadian-as plot --event analysis-run/events/EVENT_HASH.json --out event.png
```

See the [input specification](docs/input-format.md) for `samples.tsv` fields and BAM-order mapping requirements. Output directories must not already exist, to prevent overwriting previous results. The CLI sets one BLAS thread per process; `--workers` controls parallel execution across events.

## Model and coverage

Each biological condition is analyzed separately using raw inclusion/skipping counts and effective lengths. The model uses a fixed 24-hour period, experiment-specific fixed intercepts, and sine/cosine terms on the logit PSI scale. The observed inclusion-count probability includes a `log(LI/LS)` offset. The BB distribution accounts for sampling variation and overdispersion. Zero-depth observations contribute no likelihood information, and their PSI values are not imputed.

The primary test retains the V2 likelihood-ratio statistic. A pilot of 99 null simulations estimates bias in the null dispersion, followed by a parametric bootstrap from the corrected null. The default final bootstrap budget is 199. Failed simulations are conservatively counted in the upper tail; no primary p-value is reported when the pilot or final failure rate exceeds the predefined threshold. The historical numerical implementation was retained with only relative-import changes, accompanied by source SHA256 hashes and same-seed regression tests.

In the TempCompen development dataset, r1/r2 identify two experiments spanning the time series. This version does not estimate temporal autocorrelation, search over periods, test differential rhythmicity between strains or temperatures, or jointly model gene expression. Conditioning on the observed event depth N reflects the loss of information at low read counts, but cannot recover missing phases or resolve PSI-dependent nonrandom missingness.

Coverage eligibility preserves the V2 implementation byte for byte and uses only N, time, and experiment labels. Additional diagnostics describe depth by absolute time and phase, phases with repeated zero coverage across cycles, effective sample size, phase gaps, and the harmonic design condition number. These additional diagnostics do not change the filtering rules.

## Four outcome states

| Status | Meaning |
|---|---|
| `rhythm_candidate` | Coverage passes and the test is numerically reliable, with raw p ≤ alpha; an exploratory candidate. |
| `not_detected` | Coverage passes and the test is reliable, but the threshold is not reached; this does not establish absence of rhythm. |
| `coverage_insufficient` | The original coverage rules are not satisfied; the primary p-value is missing. |
| `model_unresolved` | Coverage passes, but fitting or resampling cannot produce a reliable test; the primary p-value is missing. |

Unidentifiable cases, including support for only one splice form, also receive `model_unresolved`. Adequate total coverage does not guarantee that changes in splicing proportions can be estimated reliably. The software does not assign a “confirmed nonrhythmic” label.

With 199 final bootstrap draws, p-value resolution is 0.005. Increasing the budget reduces Monte Carlo error but does not repair model misspecification. The Monte Carlo interval in each event JSON describes only final-simulation tail probability conditional on the fitted null. It excludes pilot, parameter, and model uncertainty and is not a confidence interval for the complete rhythmicity test. Events near the threshold are flagged with `threshold_ambiguous`; this diagnostic does not automatically change their classification.

## Traceable outputs

- `results.tsv`: one row per input event, with status, reason, available p-value, effect size, phase, and event JSON path.
- `coverage.tsv`: per-event metrics and failure reasons from the frozen coverage rules.
- `events/*.json`: input counts, times, experiments, effective lengths, original annotations, additional coverage diagnostics, fits, bootstrap records, Monte Carlo diagnostics, and predicted curves by experiment.
- `manifest.json`: SHA256 hashes of inputs and accompanying provenance files, parameters, software/dependency versions, code and output hashes, status counts, and run completion state.

Event seeds are derived from the root seed and an unambiguous `record_key` using SHA256, independently of event input order and process count. Invalid inputs are rejected before execution. An unexpected interruption marks the run as failed rather than complete. Historical p=1 sentinel values from numerical failures are retained only in internal diagnostic objects and are never used as primary results.

## Validation and next steps

See the [validation notes](validation/VALIDATION.md), [input specification](docs/input-format.md), [engine interface](docs/engine.md), and [development results](validation/RESULTS.md). The validation notes and development results are currently in Chinese. Both datasets are development datasets without established rhythmic-AS positive controls. Full input/coverage scans and model fitting on the fixed panel are reported separately.

A future R interface could initially use `system2()` with result-import and plotting functions, allowing R and the CLI to share the same engine and schema. Continuous BLN should provide independently calibrated sensitivity results without altering the primary BB results. Currently, `--engine bln_continuous` explicitly reports that the engine is not implemented.

The repository includes source code, tests, simulated regression fixtures, validation scripts, and summaries. Original rMATS/BAM data, the complete 332-event panel, per-event analysis outputs, and local audit files are not distributed in this public repository. References to those paths in the validation documents describe the original development environment. A machine-readable summary is available in [published_summary.json](validation/published_summary.json).

Run the development tests with:

```bash
python -m unittest discover -s tests -v
```
