"""Coverage-first execution with explicit non-negative failure states."""
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
import csv
from datetime import datetime, timezone
import multiprocessing
import math
from pathlib import Path
import platform
import sys

from . import __version__
from .coverage import diagnose, RULE, LEGACY_SHA256
from .inputs import iter_records
from .provenance import event_digest, event_seed, sha256_file, write_json

STATES = ('rhythm_candidate', 'not_detected', 'coverage_insufficient', 'model_unresolved')
FIELDS = ['record_key', 'dataset', 'condition', 'mode', 'event_type', 'event_id',
          'gene_id', 'status', 'reason', 'p_value', 'p_resolution', 'peak_time_hours',
          'psi_peak_to_trough_blockmean', 'event_seed', 'artifact']
COVERAGE_FIELDS = ['record_key', 'sufficient', 'reason', 'n_samples', 'n_positive_samples',
                   'n_positive_times', 'n_zero_times', 'n_low_times', 'n_distinct_phases',
                   'max_phase_gap', 'total_count', 'min_positive_each_repeat',
                   'effective_samples', 'harmonic_condition', 'design_rank']


def classify(coverage, analysis, alpha):
    if not coverage['sufficient']:
        return 'coverage_insufficient', coverage['reason']
    p = (analysis or {}).get('p_value')
    if (not analysis or not analysis.get('reliable') or not isinstance(p, (int, float))
            or isinstance(p, bool) or not math.isfinite(p) or not 0 <= p <= 1):
        return 'model_unresolved', (analysis or {}).get('failure_reason') or 'test_unresolved'
    if analysis['p_value'] <= alpha:
        return 'rhythm_candidate', 'exploratory_raw_p_at_or_below_alpha'
    return 'not_detected', 'rhythm_not_detected_not_evidence_of_absence'


def analyze_one(task):
    record, config = task
    coverage = diagnose(record)
    seed = event_seed(config['seed'], record['record_key'])
    analysis = None
    if coverage['sufficient']:
        from .engines import get_engine
        analysis = get_engine(config['engine']).analyze(record, config['bootstrap'], seed)
        analysis['monte_carlo'] = monte_carlo_diagnostic(analysis, config['alpha'])
    state, reason = classify(coverage, analysis, config['alpha'])
    return {'record': record, 'coverage': coverage, 'analysis': analysis,
            'status': state, 'reason': reason, 'event_seed': seed}


def monte_carlo_diagnostic(analysis, alpha):
    """Final-tail simulation error only; excludes pilot and parameter uncertainty."""
    boot = analysis.get('bootstrap', {})
    b, k, f = (boot.get(key) for key in ('bootstrap_attempted', 'bootstrap_exceed', 'bootstrap_failures'))
    if not b or k is None or f is None:
        return None
    from scipy.stats import beta
    lower = 0. if k == 0 else float(beta.ppf(.025, k, b - k + 1))
    upper = 1. if k + f == b else float(beta.ppf(.975, k + f + 1, b - k - f))
    threshold = (alpha * (b + 1) - 1) / b
    return dict(tail_probability_interval95=[lower, upper],
                threshold_ambiguous=lower <= threshold <= upper,
                p_resolution=1 / (b + 1),
                interpretation='Conditional final-bootstrap tail uncertainty; failures bounded both ways; '
                               'excludes pilot, null-parameter and model uncertainty; does not change status')


def _now():
    return datetime.now(timezone.utc).isoformat()


def execute(input_path, out_dir, *, coverage_only=False, engine='bb_bias_corrected',
            bootstrap=199, seed=20261020, workers=1, alpha=.05):
    if bootstrap < 1 or seed < 0 or workers < 1 or not 0 < alpha < 1:
        raise ValueError('Require bootstrap>=1, seed>=0, workers>=1, 0<alpha<1')
    if not coverage_only:
        from .engines import get_engine
        get_engine(engine)  # Reject a reserved/unknown engine before creating output.
    input_path, out_dir = Path(input_path).resolve(), Path(out_dir).resolve()
    seen = set()
    for record in iter_records(input_path):
        key = record['record_key']
        if key in seen:
            raise ValueError(f'Duplicate record_key: {key}')
        seen.add(key)
    if not seen:
        raise ValueError('Input contains no records')
    out_dir.mkdir(parents=True, exist_ok=False)
    config = dict(engine=engine, bootstrap=bootstrap, seed=seed, workers=workers,
                  alpha=alpha, period_hours=24, coverage_only=coverage_only)
    import numpy, scipy
    source = Path(__file__).parent
    manifest = dict(schema_version='1.0', package_version=__version__, status='running',
        started_at=_now(), input_path=str(input_path), input_sha256=sha256_file(input_path),
        n_records=len(seen), config=config, coverage_rule=RULE,
        coverage_source_sha256=LEGACY_SHA256,
        code_sha256={str(p.relative_to(source)): sha256_file(p) for p in sorted(source.rglob('*.py'))},
        dependencies={'python': sys.version, 'numpy': numpy.__version__, 'scipy': scipy.__version__,
                      'platform': platform.platform()},
        seed_rule='SHA256(UTF8(compact JSON [root_seed,record_key])) first 8 bytes big endian',
        limitations=['Exploratory raw p; no genome-wide FDR certification',
                     'Fixed 24-hour period; experiment-specific fixed intercepts',
                     'Not detected does not establish absence of rhythm',
                     'Coverage gate is operational and does not guarantee statistical power',
                     'BB calibration remains imperfect for high-dispersion heterogeneous BLN',
                     'Conditional on observed event depth; PSI-dependent missingness is not modeled'])
    for name in ('prepare_manifest.json', 'panel_manifest.json', 'selection_audit.json', 'source_snapshot_audit.json'):
        companion = input_path.parent / name
        if companion.is_file():
            import json
            manifest[name] = {'sha256': sha256_file(companion), 'content': json.loads(companion.read_text())}
    write_json(out_dir / 'manifest.json', manifest)
    counts = Counter()
    try:
        with (out_dir / 'coverage.tsv').open('w', newline='') as cover_file:
            cw = csv.DictWriter(cover_file, COVERAGE_FIELDS, delimiter='\t', extrasaction='ignore')
            cw.writeheader()
            if coverage_only:
                for record in iter_records(input_path):
                    coverage = diagnose(record)
                    cw.writerow(dict(record_key=record['record_key'], **coverage))
                    counts['sufficient' if coverage['sufficient'] else 'coverage_insufficient'] += 1
            else:
                (out_dir / 'events').mkdir()
                tasks = ((record, config) for record in iter_records(input_path))
                pool = None
                if workers > 1:
                    pool = ProcessPoolExecutor(workers, mp_context=multiprocessing.get_context('spawn'))
                try:
                    # Bounded batches avoid materializing a whole-genome future list.
                    import itertools
                    with (out_dir / 'results.tsv').open('w', newline='') as result_file:
                        rw = csv.DictWriter(result_file, FIELDS, delimiter='\t', extrasaction='ignore')
                        rw.writeheader()
                        while batch := list(itertools.islice(tasks, max(1, workers * 2))):
                            results = pool.map(analyze_one, batch) if pool else map(analyze_one, batch)
                            for result in results:
                                record, analysis = result['record'], result['analysis'] or {}
                                artifact = 'events/' + event_digest(record['record_key']) + '.json'
                                write_json(out_dir / artifact, result)
                                cw.writerow(dict(record_key=record['record_key'], **result['coverage']))
                                row = dict(record, **{k: result[k] for k in ('status', 'reason', 'event_seed')},
                                    p_value=analysis.get('p_value'),
                                    p_resolution=1 / (bootstrap + 1) if analysis.get('reliable') else None,
                                    peak_time_hours=analysis.get('peak_time_hours'),
                                    psi_peak_to_trough_blockmean=analysis.get('psi_peak_to_trough_blockmean'),
                                    artifact=artifact)
                                rw.writerow(row)
                                counts[result['status']] += 1
                                if sum(counts.values()) % 10 == 0:
                                    print(f"Completed {sum(counts.values())}/{len(seen)} events", file=sys.stderr, flush=True)
                finally:
                    if pool:
                        pool.shutdown(wait=True, cancel_futures=True)
        manifest.update(status='complete', finished_at=_now(), status_counts=dict(counts))
        manifest['output_sha256'] = {str(p.relative_to(out_dir)): sha256_file(p)
                                     for p in sorted(out_dir.rglob('*'))
                                     if p.is_file() and p.name != 'manifest.json'}
    except BaseException as exc:
        manifest.update(status='failed', finished_at=_now(), status_counts=dict(counts),
                        error=f'{type(exc).__name__}: {exc}')
        raise
    finally:
        write_json(out_dir / 'manifest.json', manifest)
    return manifest
