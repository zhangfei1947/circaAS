import hashlib
from pathlib import Path
import unittest

from circadian_as.coverage import LEGACY_SHA256
from circadian_as.provenance import event_seed
from circadian_as.runner import classify, monte_carlo_diagnostic


class WorkflowTests(unittest.TestCase):
    def test_four_states(self):
        sufficient = {'sufficient': True, 'reason': 'ok'}
        self.assertEqual(classify(sufficient, {'reliable': True, 'p_value': .01}, .05)[0], 'rhythm_candidate')
        self.assertEqual(classify(sufficient, {'reliable': True, 'p_value': .5}, .05)[0], 'not_detected')
        self.assertEqual(classify(sufficient, {'reliable': False, 'p_value': 1}, .05)[0], 'model_unresolved')
        self.assertEqual(classify({'sufficient': False, 'reason': 'missing'}, None, .05)[0], 'coverage_insufficient')
        for p in (float('nan'), float('inf'), -1, 2, True):
            self.assertEqual(classify(sufficient, {'reliable': True, 'p_value': p}, .05)[0], 'model_unresolved')

    def test_frozen_coverage_source(self):
        import circadian_as
        path = Path(circadian_as.__file__).parent / '_coverage_legacy.py'
        self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), LEGACY_SHA256)

    def test_seed_key_separation(self):
        self.assertEqual(event_seed(17, 'a'), event_seed(17, 'a'))
        self.assertNotEqual(event_seed(17, 'a'), event_seed(17, 'b'))
        self.assertNotEqual(event_seed(17, 'a'), event_seed(18, 'a'))

    def test_mc_tail_is_diagnostic_and_bounds_failed_draws(self):
        def diagnostic(k, f):
            return monte_carlo_diagnostic({'bootstrap': {'bootstrap_attempted': 199,
                'bootstrap_exceed': k, 'bootstrap_failures': f}}, .05)
        self.assertTrue(diagnostic(9, 0)['threshold_ambiguous'])
        self.assertEqual(diagnostic(0, 0)['tail_probability_interval95'][0], 0)
        self.assertEqual(diagnostic(199, 0)['tail_probability_interval95'][1], 1)
        self.assertGreater(diagnostic(9, 5)['tail_probability_interval95'][1],
                           diagnostic(9, 0)['tail_probability_interval95'][1])
        self.assertIsNone(monte_carlo_diagnostic({'bootstrap': {}}, .05))


if __name__ == '__main__':
    unittest.main()
