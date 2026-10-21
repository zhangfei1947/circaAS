# BB engine implementation and limits

The public engine is Python, preserving the numerical implementation used in the
V2/V3–V6 benchmarks. It does not depend on any sibling research directory.
An R implementation has not been substituted for the validated algorithm.

## Model

For each event and biological condition, inclusion counts conditional on total
informative reads follow a beta-binomial distribution. The logit mean has the
rMATS effective-length offset `log(LI / LS)`, experiment intercepts, and a fixed
24-hour cosine and sine pair. `reps` identifies complete time-series experiments,
not independent per-time-point replicate labels. The event's dispersion is shared
across observations. The engine does not jointly model gene abundance and does
not impute zero-depth observations.

The observed statistic is the original V2 likelihood-ratio statistic. A pilot of
99 null draws estimates bias in the null dispersion, then the bootstrap generator
uses `clip(2 * rho_hat - mean(rho_pilot), 0, .98)` and a refitted null mean. The final
bootstrap refits the original V2 statistic. `SeedSequence(seed).spawn(2)` splits
pilot and final streams, unchanged from V3. The threshold for unsuccessful pilot
or final fits remains 5%. Failed final fits count as exceedances in the upper-bound
p-value. This is approximate plug-in calibration, not an exact or double-bootstrap
test. No new shrinkage prior, filtration, or empirical p-value remapping is added.

## Public contract

```python
from circadian_as.engines import get_engine
result = get_engine('bb_bias_corrected').analyze(record, n_boot=199, seed=173)
```

`record` requires `I`, `S`, `N`, `times`, `reps`, `LI`, and `LS`. Counts and times
are nonempty equal-length vectors. Counts are finite nonnegative integers and
`N == I + S` for every sample. `LI`/`LS` may be positive scalars or per-sample
vectors. Invalid inputs raise `ValueError`; they are never silently discarded.
A supplied `period` or `period_hours` must be 24. A budget must be a positive
integer; a seed must be a nonnegative integer.

Results contain:

- `converged`: observed model fit valid, independently of resampling reliability.
- `reliable`: observed fit and pilot/final bootstrap valid under the frozen rule.
- `p_value`: calibrated p-value only when reliable; otherwise `None`.
- `statistic`, `peak_time_hours`, `psi_peak_to_trough_blockmean`: descriptive
  observed-fit estimates when the fit converges. They do not assert rhythmicity.
- `fit`, `bootstrap`: full original diagnostic dictionaries. Their historical
  sentinel `p_bootstrap=1` or `p_asymptotic=1` on failure **must not be interpreted
  as successful tests**. Consumers must use public `p_value` and `reliable`.
- `predictions`: reference PSI excluding the effective-length offset on a
  0–24-hour, 0.25-hour grid, by fitted experiment. Null predictions are also
  included. A block with only zero-depth observations is listed under
  `unavailable_blocks`; no intercept or PSI is invented for it.
- `failure_reason`: numerical/fit/resampling explanation, or `None` on success.

NumPy values become JSON scalars/lists. Nonfinite diagnostic values become JSON
null. No JSON NaN or Infinity is emitted. The engine itself makes no coverage gate
or final biological classification; the caller must retain coverage-insufficient,
not-detected, and numerical-unresolved as separate states.

`Engine` is a runtime-checkable protocol. `get_engine('bln_continuous')` reserves
the same interface and explicitly raises `NotImplementedError` in this release.
There is no silent BB fallback and no fake BLN result.

## Statistical interpretation

The default budget 199 has minimum attainable calibrated p-value 0.005. A larger
budget can improve Monte Carlo resolution but does not itself correct model
misspecification. High-dispersion, repeated missing-phase patterns have not
obtained uniformly valid 5% error control in development benchmarks. This release
is for coverage-aware candidate analysis, not a validated genome-wide discovery
claim. Bootstrap failure is an unresolved test, never evidence for no rhythm.

## Vendored provenance and parity

`src/circadian_as/engines/_legacy/provenance.json` records SHA-256 hashes of the
three historical sources and distributed copies. `bb_model.py` and
`small_sample.py` are byte-identical copies. In `bb_bias_corrected.py` only import
plumbing changes: remove `Path`/`importlib.util`, replace the dynamic sibling-path
module loader with `from . import bb_model as base`, and make `small_sample`
relative. Model formulas, optimizer options, boundaries, RNG streams, pilot size,
and failure rules are untouched.

`_legacy/parity_reference.json` contains two reproducible regression fixtures,
one per actual dataset coverage pattern, with newly simulated counts from the
frozen V6 heldout set. They are **simulation fixtures, not real biological AS
findings**. Reference values were produced by importing the old V3 engine, fitting
and running pilot 99 + final 199 with the same record/seed. All fit and bootstrap
fields were exactly equal to the vendored engine in that validation environment.
The packaged regression test allows small optimizer floating-point differences
across SciPy/BLAS builds while requiring exact discrete diagnostics and p-values.
Zero-depth, single-splice, low-depth, invalid-input, and failure-semantics tests
cover public behavior independently of this numerical regression test.

Parity compares identical sample order in the old and new engine. The engine does
not reorder its input. The software importer may canonically sort samples for
stable runs; changed ordering also changes which random draws attach to which
sample, so finite-bootstrap p-values need not reproduce a historical run performed
in another order, even with the same seed. Record ordering, seed, budget, software
version, and input hashes must all be retained by the workflow.
