"""Input checks exercise alignment and rejection of plausible silent corruption."""
import copy
import json
from pathlib import Path
import tempfile
import unittest

from circadian_as.inputs import iter_records, prepare_rmats, validate_record


class InputTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.counts = self.root / 'rmats'
        self.counts.mkdir()
        self.metadata = self.root / 'samples.tsv'
        self.bam1 = self.root / 'b1.txt'
        self.bam2 = self.root / 'b2.txt'
        self.bam1.write_text('/x/b.bam,/x/a.bam\n')
        self.bam2.write_text('/x/c.bam\n')
        self.meta_header = 'dataset\tsample_id\tcondition\ttime\treplicate\tbam\n'
        self.meta_rows = [
            'test\ta\tWT\t0\tr1\t/x/a.bam\n',
            'test\tb\tWT\t4\tr1\t/x/b.bam\n',
            'other\ta\tWT\tNaN\tr1\t/ignored/a.bam\n',
        ]
        self.metadata.write_text(self.meta_header + ''.join(self.meta_rows))
        self.raw = self.counts / 'JC.raw.input.RI.txt'
        self.ann = self.counts / 'fromGTF.RI.txt'
        self.raw_header = 'ID\tIJC_SAMPLE_1\tSJC_SAMPLE_1\tIJC_SAMPLE_2\tSJC_SAMPLE_2\tIncFormLen\tSkipFormLen\n'
        self.raw_rows = ['007\t4,1\t6,2\t\t\t100\t50\n', '8\t0,0\t0,0\t\t\t100\t50\n']
        self.raw.write_text(self.raw_header + ''.join(self.raw_rows))
        self.ann.write_text('ID\tGeneID\tgeneSymbol\tchr\n007\t"gene\t"symbol"\tchr1\n8\t"g2"\t"s2"\tchr1\n')

    def tearDown(self):
        self.temp.cleanup()

    def prepare(self, name='prepared', bam2=None, **kwargs):
        return prepare_rmats(self.counts, self.metadata, self.bam1, bam2,
                             self.root / name, 'test', event_types=['RI'], **kwargs)

    def records(self, name='prepared'):
        return list(iter_records(self.root / name / 'records.jsonl'))

    def test_literal_quotes_ids_and_bam_order(self):
        manifest = self.prepare()
        records = self.records()
        self.assertEqual(manifest['record_count'], 2)
        self.assertEqual(records[0]['event_id'], 'RI:007')
        self.assertEqual(records[0]['I'], [1, 4])
        self.assertEqual(records[0]['S'], [2, 6])
        self.assertEqual(records[0]['annotation']['GeneID'], '"gene')
        self.assertEqual(records[0]['gene_id'], '"gene')
        self.assertEqual(records[1]['N'], [0, 0])
        self.assertEqual(records[0]['source_positions'][0]['vector_index'], 1)
        self.assertEqual(json.loads(records[0]['record_key']), ['test', 'WT', 'JC', 'RI', '007'])

    def test_metadata_reordering_is_invariant(self):
        self.prepare('one')
        self.metadata.write_text(self.meta_header + ''.join(reversed(self.meta_rows)))
        self.prepare('two')
        for name in ('records.jsonl', 'samples.jsonl', 'annotations.jsonl'):
            self.assertEqual((self.root/'one'/name).read_bytes(), (self.root/'two'/name).read_bytes())

    def test_group2_explicit_mapping(self):
        self.metadata.write_text(self.meta_header + ''.join(self.meta_rows) + 'test\tc\tWT\t8\tr1\t/x/c.bam\n')
        self.raw.write_text(self.raw_header + '007\t4,1\t6,2\t9\t10\t100\t50\n')
        self.prepare(bam2=self.bam2)
        self.assertEqual(self.records()[0]['I'], [1, 4, 9])
        self.assertEqual(self.records()[0]['source_positions'][2]['group'], 2)

    def test_group2_never_ignored(self):
        self.raw.write_text(self.raw_header + '007\t4,1\t6,2\t9\t10\t100\t50\n')
        with self.assertRaisesRegex(ValueError, 'no BAM group mapping'):
            self.prepare()

    def test_no_basename_guessing(self):
        self.bam1.write_text('/different/b.bam,/x/a.bam')
        with self.assertRaisesRegex(ValueError, 'map exactly'):
            self.prepare()

    def test_invalid_counts_and_lengths(self):
        for changes in [('4,-1','6,2','100'), ('4,1.5','6,2','100'), ('4','6,2','100'), ('4,1','6,2','0')]:
            with self.subTest(changes=changes):
                self.raw.write_text(self.raw_header + f'007\t{changes[0]}\t{changes[1]}\t\t\t{changes[2]}\t50\n')
                with self.assertRaises(ValueError):
                    self.prepare()
                self.assertFalse((self.root/'prepared').exists())

    def test_missing_annotation_and_duplicate_ids(self):
        for body in ['missing\t4,1\t6,2\t\t\t100\t50\n', self.raw_rows[0]*2]:
            self.raw.write_text(self.raw_header + body)
            with self.assertRaises(ValueError):
                self.prepare()
        self.raw.write_text(self.raw_header + self.raw_rows[0])
        self.ann.write_text(self.ann.read_text() + '007\tg\ts\tchr1\n')
        with self.assertRaisesRegex(ValueError, 'duplicate annotation'):
            self.prepare()

    def test_bad_metadata(self):
        modifications = [
            self.meta_rows[0] + self.meta_rows[0] + self.meta_rows[1],
            self.meta_rows[0].replace('\t0\t', '\tNaN\t') + self.meta_rows[1],
            self.meta_rows[0] + self.meta_rows[1].replace('\t4\t', '\t0\t'),
            self.meta_rows[0],
        ]
        for body in modifications:
            with self.subTest(body=body):
                self.metadata.write_text(self.meta_header + body)
                with self.assertRaises(ValueError):
                    self.prepare()

    def test_output_never_overwritten_and_jcec_explicit(self):
        self.prepare()
        with self.assertRaises(FileExistsError):
            self.prepare()
        (self.counts/'JCEC.raw.input.RI.txt').write_text(self.raw_header + '007\t8,3\t6,2\t\t\t100\t50\n')
        self.prepare('jcec', mode='JCEC')
        self.assertEqual(self.records('jcec')[0]['I'], [3, 8])
        self.assertEqual(self.records('jcec')[0]['mode'], 'JCEC')
        with self.assertRaises(ValueError):
            self.prepare('invalid', mode='BOTH')

    def test_external_record_validation(self):
        self.prepare()
        original = self.records()[0]
        changes = [('N', [100, 10]), ('I', [True, 4]), ('times', [float('inf'), 4]), ('LI', -1), ('reps', ['r1']), ('times', [0, 0]), ('sample_ids', [[], 'b'])]
        for key, value in changes:
            record = copy.deepcopy(original)
            record[key] = value
            with self.subTest(key=key), self.assertRaises(ValueError):
                validate_record(record)
        malformed = copy.deepcopy(original)
        malformed['source_positions'][0]['group'] = True
        with self.assertRaises(ValueError):
            validate_record(malformed)
        path = self.root/'duplicates.jsonl'
        path.write_text(json.dumps(original)+'\n'+json.dumps(original)+'\n')
        with self.assertRaisesRegex(ValueError, 'duplicate record identity'):
            list(iter_records(path))


if __name__ == '__main__':
    unittest.main()
