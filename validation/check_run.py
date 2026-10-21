#!/usr/bin/env python3
"""Independent output/semantics audit; no fitting or package imports.

Compare two run directories with --compare to check deterministic event output
across serial/parallel execution or input record ordering.
"""
import argparse
import collections
import csv
import hashlib
import json
import math
from pathlib import Path

STATES = {'rhythm_candidate', 'not_detected', 'coverage_insufficient', 'model_unresolved'}


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1024*1024), b''):
            h.update(b)
    return h.hexdigest()


def load_table(path):
    with path.open() as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    mapped = {r['record_key']: r for r in rows}
    assert len(mapped) == len(rows), (path, 'duplicate records')
    return mapped


def read_json(path):
    def bad(value):
        raise ValueError(f'Non-standard JSON numeric constant {value}: {path}')
    return json.loads(path.read_text(), parse_constant=bad)


def audit(directory):
    manifest = read_json(directory/'manifest.json')
    assert manifest['status'] == 'complete', 'Run incomplete/failed'
    assert manifest['config']['coverage_only'] is False
    source=Path(manifest['input_path'])
    assert sha(source)==manifest['input_sha256'], 'input SHA mismatch'
    records = {}
    for line in source.read_text().splitlines():
        r=json.loads(line)
        assert r['record_key'] not in records
        records[r['record_key']]=r
    results=load_table(directory/'results.tsv'); coverage=load_table(directory/'coverage.tsv')
    assert records.keys()==results.keys()==coverage.keys()
    assert manifest['n_records']==len(records)
    for name, checksum in manifest['output_sha256'].items():
        assert sha(directory/name)==checksum, (name, 'output SHA mismatch')
    assert manifest['code_sha256'], 'Missing code provenance'
    for name in ('prepare_manifest.json', 'panel_manifest.json', 'selection_audit.json', 'source_snapshot_audit.json'):
        if name in manifest:
            companion=source.parent/name
            assert sha(companion)==manifest[name]['sha256']
            assert read_json(companion)==manifest[name]['content']
    code_root=Path(__file__).resolve().parents[1]/'src/circadian_as'
    code_match={name:sha(code_root/name)==checksum for name,checksum in manifest['code_sha256'].items()}
    counts=collections.Counter(); by_dataset=collections.defaultdict(collections.Counter); events={}; checked=0
    for key, record in records.items():
        row=results[key]; c=coverage[key]
        digest=hashlib.sha256(key.encode()).hexdigest()
        assert row['artifact']==f'events/{digest}.json'
        event=read_json(directory/row['artifact']); events[key]=event
        assert event['record']==record
        seed_payload=json.dumps([manifest['config']['seed'],key],ensure_ascii=False,separators=(',',':')).encode()
        seed=int.from_bytes(hashlib.sha256(seed_payload).digest()[:8],'big')
        assert event['event_seed']==int(row['event_seed'])==seed
        status=event['status']; assert row['status']==status and status in STATES
        assert row['reason']==event['reason'] and row['reason']
        cov=event['coverage']; analysis=event['analysis']
        assert c['sufficient']==str(cov['sufficient'])
        assert c['reason']==cov['reason']
        for name,value in cov.items():
            if name not in c or name in ('reason','sufficient'): continue
            if value is None:
                assert c[name] in ('', 'inf', 'nan')
            elif isinstance(value,(int,float)):
                assert float(c[name])==value
        if not cov['sufficient']:
            assert status=='coverage_insufficient' and analysis is None
        else:
            assert status!='coverage_insufficient' and analysis is not None
            p=analysis['p_value']
            reliable=analysis['reliable'] and analysis['converged'] and p is not None
            if not reliable:
                assert status=='model_unresolved' and p is None
            else:
                assert isinstance(p,(int,float)) and math.isfinite(p) and 0<p<=1
                assert status==('rhythm_candidate' if p<=manifest['config']['alpha'] else 'not_detected')
                assert float(row['p_value'])==p
                assert float(row['p_resolution'])==1/(manifest['config']['bootstrap']+1)
                assert abs(p*(manifest['config']['bootstrap']+1)-round(p*(manifest['config']['bootstrap']+1)))<1e-8
                if status=='not_detected':
                    assert 'not_evidence_of_absence' in row['reason']
                checked+=1
        if status in {'model_unresolved','coverage_insufficient'}:
            assert row['p_value']=='' and row['p_resolution']==''
        if analysis and analysis.get('monte_carlo') is not None:
            mc=analysis['monte_carlo']; lo,hi=mc['tail_probability_interval95']
            b=analysis['bootstrap']['bootstrap_attempted']
            threshold=(manifest['config']['alpha']*(b+1)-1)/b
            assert 0<=lo<=hi<=1
            assert mc['p_resolution']==1/(b+1)
            assert mc['threshold_ambiguous']==(lo<=threshold<=hi)
            assert 'does not change status' in mc['interpretation']
        for identity in ('dataset','condition','mode','event_type','event_id','gene_id'):
            assert row[identity]==record[identity]
        counts[status]+=1; by_dataset[record['dataset']][status]+=1
    assert dict(counts)==manifest['status_counts']
    assert len(list((directory/'events').glob('*.json')))==len(records)
    return dict(run=str(directory.resolve()),n_records=len(records),status_counts=dict(counts),
                status_by_dataset={k:dict(v) for k,v in by_dataset.items()},
                successful_p_values_checked=checked,input_sha256=sha(source),
                manifest_sha256=sha(directory/'manifest.json'),
                all_record_seed_and_artifact_checks=True,
                current_code_hashes_match=all(code_match.values()),
                code_files_different_from_recorded_version=[k for k,v in code_match.items() if not v]),events


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('run',type=Path)
    p.add_argument('--compare',type=Path)
    p.add_argument('--out',type=Path)
    a=p.parse_args(); report,events=audit(a.run)
    if a.compare:
        other,events2=audit(a.compare)
        assert events.keys()==events2.keys(), 'Comparison event IDs differ'
        for key in events:
            assert events[key]==events2[key], ('Non-deterministic result', key)
        report['comparison']={'run':other['run'],'all_event_json_identical':True,'events':len(events)}
    if a.out:
        a.out.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__':
    main()
