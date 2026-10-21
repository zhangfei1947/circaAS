"""Write four synthetic examples; these are not a statistical validation dataset."""
import argparse
import json
import math
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', required=True)
    args = parser.parse_args()
    path = Path(args.out)
    times = list(range(0, 48, 4)) * 2
    reps = ['r1'] * 12 + ['r2'] * 12
    with path.open('x') as stream:
        for index, label in enumerate(('harmonic_example', 'constant_example', 'no_reads', 'one_splice_form'), 1):
            total = [0 if label == 'no_reads' else 100] * 24
            inc = [round(n * (.5 + .25 * math.cos(2 * math.pi * t / 24)))
                   if label == 'harmonic_example' else (n if label == 'one_splice_form' else n // 2)
                   for n, t in zip(total, times)]
            record = dict(dataset='synthetic_demo', condition='example', mode='JC', event_type='SE',
                          event_id=f'SE:{index}', profile_id=f'synthetic_demo/example/JC/SE:{index}',
                          record_key=json.dumps(['synthetic_demo', 'example', 'JC', 'SE', str(index)], separators=(',', ':')),
                          sample_ids=[f'{r}_t{t:02d}' for r, t in zip(reps, times)], times=times, reps=reps,
                          I=inc, S=[n-i for n, i in zip(total, inc)], N=total, LI=100, LS=100,
                          gene_id=f'demo_gene_{index}', gene_symbol=label, annotation={'ID': str(index)},
                          note='Synthetic software demo; names describe inputs, not guaranteed inferential outcomes')
            stream.write(json.dumps(record, allow_nan=False) + '\n')


if __name__ == '__main__':
    main()
