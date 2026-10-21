#!/usr/bin/env python3
"""Independently audit full imports against raw snapshots, then select frozen V2 IDs.

This script deliberately does not import the new package or old preparation code.
It validates all raw JC event-condition records, not only the development panel.
"""
import argparse
import collections
import csv
import hashlib
import json
from pathlib import Path

ETYPES = ('SE', 'A3SS', 'A5SS', 'MXE', 'RI')
FIELDS = ('dataset', 'condition', 'event_type', 'event_id', 'gene_id',
          'gene_symbol', 'I', 'S', 'N', 'times', 'reps', 'LI', 'LS', 'sample_ids')


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()


def table(path):
    lines = Path(path).read_text().splitlines()
    head = lines[0].split('\t')
    assert len(head) == len(set(head)), (path, 'duplicate columns')
    rows = []
    for ln, line in enumerate(lines[1:], 2):
        values = line.split('\t')
        assert len(values) == len(head), (path, ln, 'ragged row')
        rows.append(dict(zip(head, values)))
    return rows


def indexed(path):
    rows = table(path)
    d = {r['ID']: r for r in rows}
    assert len(d) == len(rows), (path, 'duplicate IDs')
    return d


def audit_samples(root):
    samples = table(root / 'inputs/samples.tsv')
    assert len({(r['dataset'], r['sample_id']) for r in samples}) == len(samples)
    assert len({(r['dataset'], r['condition'], r['time'], r['replicate']) for r in samples}) == len(samples)
    plan = json.loads((root / 'inputs/neurospora/inputs/plan.json').read_text())['samples']
    nmap = {r['library_id']: r for r in plan}
    with (root / 'inputs/tempcompen/Table_S4_sample_metadata.tsv').open() as f:
        tm = {r['sample_id']: r for r in csv.DictReader(f, delimiter='\t')}
    used = [root / 'inputs/samples.tsv', root / 'inputs/neurospora/inputs/plan.json',
            root / 'inputs/tempcompen/Table_S4_sample_metadata.tsv']
    bamlists = {'neurospora': root / 'inputs/neurospora/inputs/all_rna_bams.txt',
                'tempcompen': root / 'inputs/tempcompen/bamfile.all.txt'}
    for ds, path in bamlists.items():
        ss = [r for r in samples if r['dataset'] == ds]
        bams = path.read_text().strip().strip(',').split(',')
        assert [r['bam'] for r in ss] == bams, (ds, 'BAM order mismatch')
        used.append(path)
        for r in ss:
            if ds == 'neurospora':
                original = nmap[r['sample_id']]
                assert r['condition'] == original['genotype']
                assert float(r['time']) == original['time_hours']
                assert r['replicate'] == str(original['replicate'])
                assert r['bam'] == original['bam_path']
            else:
                original = tm[r['sample_id']]
                assert r['condition'] == original['condition']
                assert float(r['time']) == (int(original['day'][1:])-4)*24+int(original['phase_hour'])
                assert r['replicate'] == original['replicate']
                assert r['bam'] == original['bam_path']
    assert set(r['dataset'] for r in samples) == set(bamlists)
    return samples, used


def expected_records(root, samples, used):
    for ds in ('neurospora', 'tempcompen'):
        ss = [r for r in samples if r['dataset'] == ds]
        folder = root / 'inputs' / ('neurospora/output' if ds == 'neurospora' else ds)
        for et in ETYPES:
            rawpath = folder / f'JC.raw.input.{et}.txt'
            annpath = folder / f'fromGTF.{et}.txt'
            used.extend([rawpath, annpath])
            raw, ann = indexed(rawpath), indexed(annpath)
            assert raw.keys() == ann.keys(), (ds, et, 'annotation/count IDs differ')
            for eid, r in raw.items():
                inc = list(map(int, r['IJC_SAMPLE_1'].split(',')))
                skp = list(map(int, r['SJC_SAMPLE_1'].split(',')))
                assert len(inc) == len(skp) == len(ss)
                assert all(x >= 0 for x in inc + skp)
                assert r['IJC_SAMPLE_2'] in ('', 'NA') and r['SJC_SAMPLE_2'] in ('', 'NA')
                li, ls = int(r['IncFormLen']), int(r['SkipFormLen'])
                assert li > 0 and ls > 0
                for cond in dict.fromkeys(r['condition'] for r in ss):
                    indices = [i for i, s in enumerate(ss) if s['condition'] == cond]
                    inc2 = [inc[i] for i in indices]; skp2 = [skp[i] for i in indices]
                    yield dict(profile_id=f'{ds}|{cond}|{et}:{eid}', dataset=ds, condition=cond,
                               event_type=et, event_id=f'{et}:{eid}', gene_id=ann[eid]['GeneID'].strip('"'),
                               gene_symbol=ann[eid]['geneSymbol'].strip('"'), annotation=ann[eid], source_positions=[dict(group=1,vector_index=i,bam=ss[i]['bam']) for i in indices], I=inc2, S=skp2,
                               N=[i+s for i,s in zip(inc2, skp2)], times=[float(ss[i]['time']) for i in indices],
                               reps=[ss[i]['replicate'] for i in indices], LI=li, LS=ls,
                               sample_ids=[ss[i]['sample_id'] for i in indices])


def legacy_id(r):
    return f"{r['dataset']}|{r['condition']}|{r['event_id']}"


def compare(got, expected, context):
    indices = [got['sample_ids'].index(sid) for sid in expected['sample_ids']]
    assert len(set(indices)) == len(got['sample_ids']) == len(expected['sample_ids'])
    arrays = {'I', 'S', 'N', 'times', 'reps', 'sample_ids'}
    for key in FIELDS:
        value = [got[key][i] for i in indices] if key in arrays else got.get(key)
        # The new importer preserves unmatched quotes in display gene labels;
        # the historical panel stripped them. Full raw annotation is audited below.
        if key in {'gene_id', 'gene_symbol'}:
            value = value.strip('"')
        assert value == expected[key], (context, expected['profile_id'], key, value, expected[key])
    if context == 'raw':
        assert got['annotation'] == expected['annotation'], (expected['profile_id'], 'raw annotations')
        assert [got['source_positions'][i] for i in indices] == expected['source_positions']
        assert got['mode'] == 'JC'
        assert json.loads(got['record_key']) == [got['dataset'], got['condition'], 'JC', got['event_type'], got['event_id'].split(':', 1)[1]]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--records', nargs='+', type=Path, required=True, help='Full prepared records.jsonl, one or both datasets')
    p.add_argument('--source-root', type=Path, required=True)
    p.add_argument('--reference-panel', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True, help='New validation directory; will not overwrite')
    a = p.parse_args()
    assert not a.output.exists(), 'Refusing to replace validation directory'
    samples, used = audit_samples(a.source_root)
    refs = json.loads(a.reference_panel.read_text())
    frozen = {r['profile_id']: r for r in refs}
    assert len(frozen) == len(refs) == 332
    actual = {}
    for path in a.records:
        with path.open() as f:
            for line in f:
                r = json.loads(line)
                assert legacy_id(r) not in actual, ('duplicate profile', legacy_id(r))
                actual[legacy_id(r)] = r
    counts = collections.Counter(); seen = set(); selected = {}
    for expected in expected_records(a.source_root, samples, used):
        key = expected['profile_id']; seen.add(key)
        assert key in actual, ('missing imported event', key)
        compare(actual[key], expected, 'raw')
        counts[expected['dataset']] += 1
        if key in frozen:
            compare(actual[key], frozen[key], 'V2 fixed panel')
            selected[key] = actual[key]
    assert seen == actual.keys(), 'Unexpected imported event IDs'
    assert len(seen) == 111920 and len(selected) == 332
    source_manifest = a.source_root / 'qc/source_manifest.json'
    frozen_sources = {x['path']: x['sha256'] for x in json.loads(source_manifest.read_text())}
    snapshot_checks = []
    for path in dict.fromkeys(used):
        rel = str(path.relative_to(a.source_root))
        if rel in frozen_sources:
            assert sha(path) == frozen_sources[rel], ('source snapshot drift', rel)
            snapshot_checks.append(rel)
    paths = list(dict.fromkeys(used + a.records + [a.reference_panel, source_manifest]))
    a.output.mkdir(parents=True)
    panelpath = a.output / 'records.jsonl'
    with panelpath.open('w') as f:
        for reference in refs:
            f.write(json.dumps(selected[reference['profile_id']], allow_nan=False, separators=(',', ':')) + '\n')
    summary = dict(source_snapshot_hashes_verified=snapshot_checks, full_raw_import_count=len(seen), full_by_dataset=dict(counts),
                   selected_count=len(selected), selected_by_dataset=dict(collections.Counter(r['dataset'] for r in selected.values())),
                   comparison_fields=list(FIELDS) + ['annotation', 'source_positions', 'mode', 'record_key'],
                   array_alignment='Compared by sample_id; new canonical sorting need not equal raw BAM order.',
                   annotation_note='Full literal annotation preserved. Historical V2 display gene fields stripped quotes; matched after strip only for legacy comparison.', sample_counts=dict(collections.Counter(r['dataset'] for r in samples)),
                   sample_order_verified_against='original BAM lists; time/experiment labels checked against plan.json and Table_S4',
                   selection='Exact historical V2 332 profile IDs, fixed before new software results; no new p/PSI selection.',
                   limits='Full raw import audit is separate from selected-panel fit test. Both datasets are development data; no biological ground truth or error-rate estimate.',
                   sources=[dict(path=str(path.resolve()), sha256=sha(path), bytes=path.stat().st_size) for path in paths],
                   output=dict(path=str(panelpath.resolve()), sha256=sha(panelpath)),
                   selected_profile_ids=[r['profile_id'] for r in refs])
    (a.output / 'selection_audit.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps({k: summary[k] for k in ['full_raw_import_count', 'full_by_dataset', 'selected_count', 'selected_by_dataset']}))

if __name__ == '__main__':
    main()
