"""Stable engine interface. Coverage decisions belong to the calling workflow."""
from __future__ import annotations

from collections.abc import Mapping
from typing import Any, Protocol, runtime_checkable

import numpy as np

from ._legacy import bb_bias_corrected as _bb


@runtime_checkable
class Engine(Protocol):
    name: str

    def analyze(self, record: Mapping[str, Any], n_boot: int = 199,
                seed: int = 173) -> dict[str, Any]: ...


def _json_safe(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return [_json_safe(item) for item in value.tolist()]
    if isinstance(value, np.generic):
        return _json_safe(value.item())
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def _integer(value: Any, name: str, minimum: int = 0) -> int:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)):
        raise ValueError(f'{name} must be an integer >= {minimum}')
    if value < minimum:
        raise ValueError(f'{name} must be an integer >= {minimum}')
    return int(value)


def _numeric(value: Any, name: str) -> np.ndarray:
    if any(isinstance(item, (bool, np.bool_)) for item in np.asarray(value, dtype=object).flat):
        raise ValueError(f'{name} must not contain booleans')
    raw = np.asarray(value)
    if raw.dtype.kind not in 'iuf' or raw.dtype.kind == 'b':
        raise ValueError(f'{name} must contain numeric values (not strings or booleans)')
    result = raw.astype(float)
    if not np.all(np.isfinite(result)):
        raise ValueError(f'{name} must contain only finite values')
    return result


def validate_record(record: Mapping[str, Any]) -> dict[str, np.ndarray]:
    """Validate without dropping rows; zero-depth rows remain in the bootstrap."""
    missing = {'I', 'S', 'N', 'times', 'reps', 'LI', 'LS'} - set(record)
    if missing:
        raise ValueError(f'Missing engine fields: {", ".join(sorted(missing))}')
    for period_key in ('period', 'period_hours'):
        if period_key in record:
            value = _numeric(record[period_key], period_key)
            if value.ndim != 0 or float(value) != 24:
                raise ValueError('This engine supports only the fixed 24-hour period')
    arrays = {name: _numeric(record[name], name) for name in ('I', 'S', 'N', 'times')}
    length = arrays['I'].size
    if not length or any(arr.ndim != 1 or len(arr) != length for arr in arrays.values()):
        raise ValueError('I, S, N and times must be nonempty one-dimensional arrays of equal length')
    for name in ('I', 'S', 'N'):
        arr = arrays[name]
        if np.any(arr < 0) or np.any(arr != np.floor(arr)) or np.any(arr > 2**53 - 1):
            raise ValueError(f'{name} must contain nonnegative exactly representable integer counts')
    if not np.array_equal(arrays['I'] + arrays['S'], arrays['N']):
        raise ValueError('N must equal I + S in every sample')
    for name in ('LI', 'LS'):
        arr = _numeric(record[name], name)
        if arr.ndim == 0:
            arr = np.full(length, float(arr))
        if arr.ndim != 1 or len(arr) != length or np.any(arr <= 0):
            raise ValueError(f'{name} must be positive and scalar or one value per sample')
        arrays[name] = arr
    raw_reps = np.asarray(record['reps'], dtype=object)
    if raw_reps.ndim != 1 or len(raw_reps) != length:
        raise ValueError('reps must have one experiment identifier per sample')
    for value in raw_reps:
        if value is None or isinstance(value, (bool, np.bool_)):
            raise ValueError('reps must not contain missing or boolean identifiers')
        if not isinstance(value, (str, int, float, np.integer, np.floating)):
            raise ValueError('reps must contain scalar string or numeric identifiers')
        if isinstance(value, (float, np.floating)) and not np.isfinite(value):
            raise ValueError('reps must contain finite identifiers')
        if not str(value).strip():
            raise ValueError('reps must not contain empty identifiers')
    arrays['reps'] = raw_reps.astype(str)
    # Compute the same expression used in the historical workflow. Reject an
    # unrepresentable ratio instead of silently creating nonfinite offsets.
    with np.errstate(over='ignore', divide='ignore', invalid='ignore'):
        arrays['offset'] = np.log(arrays['LI'] / arrays['LS'])
    if not np.all(np.isfinite(arrays['offset'])):
        raise ValueError('LI / LS must have a finite representable log ratio')
    return arrays


class BBEngine:
    """Original BB LRT with pilot-99 bias-corrected-null calibration.

    The public p_value is absent when the fit or bootstrap is unreliable.
    Historical sentinel values are preserved only inside diagnostic subobjects.
    """
    name = 'bb_bias_corrected'

    def analyze(self, record: Mapping[str, Any], n_boot: int = 199,
                seed: int = 173) -> dict[str, Any]:
        n_boot = _integer(n_boot, 'n_boot', 1)
        seed = _integer(seed, 'seed')
        arr = validate_record(record)
        result: dict[str, Any] = {
            'engine': self.name, 'period_hours': 24.0,
            'converged': False, 'reliable': False, 'p_value': None,
            'statistic': None, 'peak_time_hours': None,
            'psi_peak_to_trough_blockmean': None,
            'fit': {}, 'bootstrap': {}, 'predictions': None,
            'failure_reason': None,
        }
        try:
            fit = _bb.fit(arr['I'], arr['N'], arr['times'], arr['offset'], arr['reps'])
            result['fit'] = fit
            result['converged'] = bool(fit['converged'])
            boot = _bb.bootstrap(arr['I'], arr['N'], arr['times'], arr['offset'],
                                 arr['reps'], n_boot=n_boot, seed=seed, observed=fit)
            result['bootstrap'] = boot
            result['reliable'] = bool(fit['converged'] and boot['bootstrap_reliable'])
            if result['reliable']:
                result['p_value'] = boot['p_bootstrap']
            if not fit['converged']:
                result['failure_reason'] = fit.get('failure_reason', 'observed_fit_invalid')
            elif not result['reliable']:
                result['failure_reason'] = boot.get('bootstrap_skipped_reason', 'bootstrap_failure_rate_above_5pct')
            if fit['converged']:
                result['statistic'] = fit['statistic']
                result['peak_time_hours'] = fit['peak_time']
                result['psi_peak_to_trough_blockmean'] = fit['psi_peak_to_trough_blockmean']
                grid = np.linspace(0, 24, 97)
                result['predictions'] = {
                    'times': grid,
                    'interpretation': 'Reference PSI, excluding effective-length offset; fitted separately for each experiment intercept',
                    'unavailable_blocks': sorted(set(arr['reps']) - set(fit['block_levels'])),
                    'blocks': [
                        {'block': block,
                         'psi': _bb.predict_psi(grid, fit, np.repeat(block, len(grid))),
                         'null_psi': _bb.predict_psi(grid, fit, np.repeat(block, len(grid)), harmonic=False)}
                        for block in fit['block_levels']
                    ],
                }
        except (ValueError, RuntimeError, FloatingPointError, np.linalg.LinAlgError) as exc:
            # Input validation precedes this block. A numerical failure is an
            # unresolved test, never evidence of absence or a public p of one.
            result['reliable'] = False
            result['p_value'] = None
            result['failure_reason'] = f'numerical_exception:{type(exc).__name__}:{exc}'
        return _json_safe(result)


def get_engine(name: str = 'bb_bias_corrected') -> Engine:
    if name == 'bb_bias_corrected':
        return BBEngine()
    if name == 'bln_continuous':
        raise NotImplementedError('Continuous BLN sensitivity analysis is a reserved interface; it is not implemented in this release')
    raise ValueError(f'Unknown engine: {name!r}')


__all__ = ['Engine', 'BBEngine', 'get_engine', 'validate_record']
