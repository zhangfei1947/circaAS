"""Numerical regression against frozen V3, plus public missing-result semantics."""
import inspect
import json
from pathlib import Path
import unittest
from unittest.mock import patch

import numpy as np

from circadian_as.engines import BBEngine, Engine, get_engine, validate_record
from circadian_as.engines import _legacy


FIXTURES = Path(_legacy.__file__).parent / 'parity_reference.json'


def simple_record():
    return {'I': [10, 15, 20, 15, 10, 5] * 2,
            'S': [20, 15, 10, 15, 20, 25] * 2,
            'N': [30] * 12, 'times': list(range(0, 48, 4)),
            'reps': ['r1'] * 6 + ['r2'] * 6, 'LI': 2., 'LS': 1.}


class EngineTests(unittest.TestCase):
    def assert_nested_equal(self, actual, expected):
        if isinstance(expected, dict):
            self.assertEqual(set(actual), set(expected))
            for key, value in expected.items():
                self.assert_nested_equal(actual[key], value)
        elif isinstance(expected, list):
            self.assertEqual(len(actual), len(expected))
            for left, right in zip(actual, expected):
                self.assert_nested_equal(left, right)
        elif isinstance(expected, float):
            # Exact on the fixture-generation platform; allow optimizer roundoff
            # across supported BLAS/SciPy versions, without relaxing counts/p.
            np.testing.assert_allclose(actual, expected, rtol=2e-7, atol=2e-8)
        else:
            self.assertEqual(actual, expected)

    def test_frozen_legacy_parity_two_actual_coverage_patterns(self):
        fixtures = json.loads(FIXTURES.read_text())
        self.assertEqual({x['dataset'] for x in fixtures}, {'neurospora', 'tempcompen'})
        for fixture in fixtures:
            with self.subTest(dataset=fixture['dataset']):
                result = BBEngine().analyze(fixture['record'], n_boot=199, seed=fixture['seed'])
                self.assert_nested_equal(result['fit'], fixture['expected_fit'])
                self.assert_nested_equal(result['bootstrap'], fixture['expected_bootstrap'])
                self.assertEqual(result['bootstrap']['pilot_attempted'], 99)
                self.assertEqual(result['bootstrap']['bootstrap_attempted'], 199)
                self.assertEqual(result['p_value'], fixture['expected_bootstrap']['p_bootstrap'])
                self.assertTrue(result['converged'])
                self.assertTrue(result['reliable'])
                json.dumps(result, allow_nan=False)
                for block in result['predictions']['blocks']:
                    self.assertEqual(len(block['psi']), 97)
                    self.assertTrue(all(0 <= p <= 1 for p in block['psi']))

    def test_zero_depth_and_single_splice_are_not_negative_results(self):
        for inclusion, skipping in ((0, 0), (0, 30), (30, 0)):
            record = simple_record()
            record.update(I=[inclusion]*12, S=[skipping]*12, N=[inclusion+skipping]*12)
            result = BBEngine().analyze(record, n_boot=1)
            self.assertFalse(result['converged'])
            self.assertFalse(result['reliable'])
            self.assertIsNone(result['p_value'])
            self.assertIsNone(result['predictions'])
            self.assertEqual(result['bootstrap']['bootstrap_attempted'], 0)
            json.dumps(result, allow_nan=False)

    def test_low_positive_depth_is_retained_without_a_hidden_gate(self):
        record = simple_record()
        record.update(I=[0, 1]*6, S=[1, 0]*6, N=[1]*12)
        validated = validate_record(record)
        self.assertEqual(len(validated['N']), 12)
        self.assertTrue(np.all(validated['N'] == 1))
        result = BBEngine().analyze(record, n_boot=1)
        self.assertEqual(result['fit']['n_samples'], 12)
        if not result['reliable']:
            self.assertIsNone(result['p_value'])

    def test_zero_rows_retained_and_unobserved_blocks_not_predicted(self):
        record = simple_record()
        record['I'] += [0]; record['S'] += [0]; record['N'] += [0]
        record['times'] += [0]; record['reps'] += ['unobserved']
        result = BBEngine().analyze(record, n_boot=1)
        self.assertTrue(result['converged'])
        self.assertEqual(result['fit']['n_zero_depth'], 1)
        self.assertEqual(result['predictions']['unavailable_blocks'], ['unobserved'])
        self.assertEqual({b['block'] for b in result['predictions']['blocks']}, {'r1', 'r2'})

    def test_invalid_inputs_raise_instead_of_dropping_samples(self):
        cases = [('I', [np.nan]*12), ('I', [-1]*12), ('I', [1.5]*12),
                 ('I', [True]+[10]*11), ('N', [10]*12), ('times', [0]*11),
                 ('LI', 0), ('LS', float('inf')), ('reps', [None]*12),
                 ('period', 22), ('period_hours', 12)]
        for name, value in cases:
            with self.subTest(field=name, value=str(value)):
                record = simple_record(); record[name] = value
                with self.assertRaises(ValueError):
                    BBEngine().analyze(record, n_boot=1)
        for budget in (0, -1, True, 1.5):
            with self.assertRaises(ValueError):
                BBEngine().analyze(simple_record(), n_boot=budget)

    def test_unreliable_bootstrap_does_not_expose_historical_p_one(self):
        with patch('circadian_as.engines._bb.bootstrap', return_value={
                'bootstrap_reliable': False, 'p_bootstrap': 1.,
                'bootstrap_skipped_reason': 'pilot_failure_rate_above_5pct'}):
            result = BBEngine().analyze(simple_record(), n_boot=199)
        self.assertTrue(result['converged'])
        self.assertFalse(result['reliable'])
        self.assertIsNone(result['p_value'])
        self.assertEqual(result['bootstrap']['p_bootstrap'], 1.)
        self.assertEqual(result['failure_reason'], 'pilot_failure_rate_above_5pct')

    def test_numerical_exception_is_unresolved(self):
        with patch('circadian_as.engines._bb.fit', side_effect=np.linalg.LinAlgError('singular')):
            result = BBEngine().analyze(simple_record(), n_boot=1)
        self.assertFalse(result['reliable'])
        self.assertIsNone(result['p_value'])
        self.assertIn('numerical_exception', result['failure_reason'])

    def test_registry_and_no_runtime_history_dependency(self):
        self.assertIsInstance(get_engine(), Engine)
        with self.assertRaises(NotImplementedError):
            get_engine('bln_continuous')
        with self.assertRaises(ValueError):
            get_engine('unknown')
        for path in Path(_legacy.__file__).parent.glob('*.py'):
            source = path.read_text()
            self.assertNotIn('spec_from_file_location', source)
            self.assertNotIn('sys.path', source)
            self.assertNotIn('parents[', source)


if __name__ == '__main__':
    unittest.main()
