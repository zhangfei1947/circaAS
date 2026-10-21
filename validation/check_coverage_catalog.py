#!/usr/bin/env python3
"""Compare complete new coverage-only output with frozen V2 coverage catalog."""
import argparse
import collections
import csv
import gzip
import json
import math
from pathlib import Path
from make_validation_panel import sha


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--runs',type=Path,nargs='+',required=True)
    p.add_argument('--reference',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    with gzip.open(a.reference,'rt') as f:
        reference={r['profile_id']:r for r in csv.DictReader(f,delimiter='\t')}
    got={}; hashes={}; counts=collections.defaultdict(collections.Counter); numeric=0
    for run in a.runs:
        manifest=json.loads((run/'manifest.json').read_text())
        assert manifest['status']=='complete' and manifest['config']['coverage_only'] is True
        assert sha(Path(manifest['input_path']))==manifest['input_sha256']
        assert sha(run/'coverage.tsv')==manifest['output_sha256']['coverage.tsv']
        hashes[str(run.resolve())]=sha(run/'manifest.json')
        with (run/'coverage.tsv').open() as f:
            for row in csv.DictReader(f,delimiter='\t'):
                ds,cond,mode,et,eid=json.loads(row['record_key'])
                assert mode=='JC'
                key=f'{ds}|{cond}|{et}:{eid}'
                assert key not in got
                got[key]=row
                old=reference[key]
                for col,value in row.items():
                    if col=='record_key':continue
                    if col in {'sufficient','reason'}:
                        assert value==old[col], (key,col,value,old[col])
                    else:
                        x=float(value) if value else math.inf; y=float(old[col])
                        assert x==y or math.isclose(x,y,rel_tol=1e-10,abs_tol=1e-8), (key,col,x,y)
                        numeric+=1
                counts[ds]['sufficient' if row['sufficient']=='True' else 'coverage_insufficient']+=1
    assert got.keys()==reference.keys()
    report=dict(n_events=len(got),status_by_dataset={k:dict(v) for k,v in counts.items()},
                numeric_comparisons=numeric,all_gate_decisions_reasons_and_diagnostics_match=True,
                frozen_catalog_sha256=sha(a.reference),run_manifest_hashes=hashes,
                scope='All raw JC event-condition records; no PSI or outcome selection; frozen V2 gate parity only, not power guarantee')
    a.out.write_text(json.dumps(report,indent=2)+'\n'); print(json.dumps(report,indent=2))

if __name__=='__main__':main()
