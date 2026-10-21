"""Frozen V2 eligibility with additional descriptive phase diagnostics."""
import numpy as np
from ._coverage_legacy import RULE, assess_coverage
from .provenance import safe

LEGACY_SHA256 = '0f7cf431ea0967fd345d07402fef6ad8693525a469ed75d9d0a48cd1cde95490'


def diagnose(record):
    n = np.asarray(record['N'], dtype=float)
    times = np.asarray(record['times'], dtype=float)
    result = assess_coverage(n, times, record['reps'])
    by_time = [{'time': float(t), 'total_count': int(n[times == t].sum()),
                'positive_samples': int((n[times == t] > 0).sum())}
               for t in np.unique(times)]
    phase = np.round(times % 24, 8)
    by_phase = []
    for p in np.unique(phase):
        cells = [row for row in by_time if abs(row['time'] % 24 - p) < 1e-7]
        by_phase.append({'phase': float(p), 'observed_absolute_times': len(cells),
                         'zero_absolute_times': sum(c['total_count'] == 0 for c in cells),
                         'total_count': int(n[phase == p].sum())})
    result.update(time_depth=by_time, phase_depth=by_phase,
                  repeated_zero_phases=[p['phase'] for p in by_phase
                      if p['observed_absolute_times'] >= 2 and p['zero_absolute_times'] >= 2],
                  diagnostic_scope='Event informative-read coverage, not measured gene expression',
                  gate_version='V2_frozen', gate_source_sha256=LEGACY_SHA256)
    return safe(result)
