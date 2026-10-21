"""Command-line entry point; heavy scientific imports follow argument parsing."""
import argparse
import json
import os
from pathlib import Path


def main(argv=None):
    parser = argparse.ArgumentParser(prog='circadian-as', description='Exploratory fixed-24h count-based AS analysis')
    from . import __version__
    parser.add_argument('--version', action='version', version=__version__)
    subs = parser.add_subparsers(dest='command', required=True)
    prep = subs.add_parser('prepare-rmats', help='Import raw rMATS counts with explicit sample/BAM ordering')
    for flag in ('counts-dir', 'samples', 'bam-list-1', 'dataset', 'out'):
        prep.add_argument('--' + flag, required=True)
    prep.add_argument('--bam-list-2')
    prep.add_argument('--mode', choices=['JC', 'JCEC'], default='JC')
    prep.add_argument('--event-types', nargs='+', choices=['SE', 'A3SS', 'A5SS', 'MXE', 'RI'])
    for name in ('coverage', 'run'):
        p = subs.add_parser(name)
        p.add_argument('--input', required=True)
        p.add_argument('--out', required=True)
        if name == 'run':
            p.add_argument('--engine', choices=['bb_bias_corrected', 'bln_continuous'], default='bb_bias_corrected')
            p.add_argument('--bootstrap', type=int, default=199)
            p.add_argument('--seed', type=int, default=20261020)
            p.add_argument('--workers', type=int, default=1)
            p.add_argument('--alpha', type=float, default=.05)
    subs.add_parser('engines')
    plot = subs.add_parser('plot', help='Plot an event JSON artifact with observed PSI and read depth')
    plot.add_argument('--event', required=True)
    plot.add_argument('--out', required=True)
    args = parser.parse_args(argv)
    # One BLAS thread per worker. Set before importing NumPy/SciPy.
    for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
        os.environ[key] = '1'
    try:
        if args.command == 'prepare-rmats':
            from .inputs import prepare_rmats
            result = prepare_rmats(args.counts_dir, args.samples, args.bam_list_1,
                                   args.bam_list_2, args.out, args.dataset, args.mode, args.event_types)
        elif args.command in ('run', 'coverage'):
            from .runner import execute
            kwargs = {k: getattr(args, k) for k in ('engine', 'bootstrap', 'seed', 'workers', 'alpha')} if args.command == 'run' else {}
            result = execute(args.input, args.out, coverage_only=args.command == 'coverage', **kwargs)
        elif args.command == 'engines':
            result = {'bb_bias_corrected': 'available: exploratory fixed-24h BB with pilot bias correction',
                      'bln_continuous': 'reserved: not implemented'}
        else:
            from .plotting import plot_event
            plot_event(args.event, args.out)
            result = {'plot': str(Path(args.out).resolve())}
        print(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False))
    except (ValueError, FileExistsError, FileNotFoundError, NotImplementedError) as exc:
        parser.exit(2, f'Error: {exc}\n')
