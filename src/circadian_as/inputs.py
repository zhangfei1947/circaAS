"""Strict, literal-tab rMATS count ingestion with explicit BAM ordering.

BAM paths are identifiers: no basename fallback or filesystem resolution is used.
"""
from __future__ import annotations

from collections import defaultdict
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import tempfile

EVENT_TYPES = ('SE', 'A5SS', 'A3SS', 'MXE', 'RI')
_COUNT_COLUMNS = ('ID', 'IJC_SAMPLE_1', 'SJC_SAMPLE_1', 'IJC_SAMPLE_2', 'SJC_SAMPLE_2', 'IncFormLen', 'SkipFormLen')
_INTEGER = re.compile(r'^[0-9]+$')


def _sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def _table(path, required):
    """Read TSV literally, preserving rMATS quotes as data, including unmatched quotes."""
    path = Path(path)
    with path.open(encoding='utf-8-sig', newline='') as stream:
        header = stream.readline().rstrip('\r\n').split('\t')
        if len(set(header)) != len(header) or any(not x for x in header):
            raise ValueError(f'{path}: empty or duplicate column name')
        missing = set(required) - set(header)
        if missing:
            raise ValueError(f'{path}: missing columns {sorted(missing)}')
        for line_no, line in enumerate(stream, 2):
            values = line.rstrip('\r\n').split('\t')
            if len(values) != len(header):
                raise ValueError(f'{path}:{line_no}: expected {len(header)} tab-separated fields; found {len(values)}')
            yield line_no, dict(zip(header, values))


def _nonempty(value, context):
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f'{context}: required nonempty string')
    return value


def _integer(value, context):
    if not _INTEGER.fullmatch(value):
        raise ValueError(f'{context}: expected nonnegative integer, got {value!r}')
    out = int(value)
    if out > 2**53 - 1:
        raise ValueError(f'{context}: count exceeds exact floating-point integer range')
    return out


def _bam_order(path):
    if path is None:
        return []
    text = Path(path).read_text(encoding='utf-8-sig').strip()
    if not text:
        raise ValueError(f'{path}: empty BAM order file')
    # Standard rMATS comma list, optionally wrapped at line boundaries.
    values = [v.strip() for line in text.splitlines() for v in line.split(',')]
    if any(not x for x in values) or len(set(values)) != len(values):
        raise ValueError(f'{path}: empty or duplicate BAM identifier')
    return values


def _samples(path, dataset, groups):
    rows = []
    for line_no, row in _table(path, ('sample_id', 'condition', 'time', 'replicate', 'bam')):
        if 'dataset' in row and row['dataset'] != dataset:
            continue
        for name in ('sample_id', 'condition', 'replicate', 'bam'):
            _nonempty(row[name], f'{path}:{line_no}:{name}')
        try:
            time = float(row['time'])
        except ValueError as exc:
            raise ValueError(f'{path}:{line_no}: invalid time') from exc
        if not math.isfinite(time):
            raise ValueError(f'{path}:{line_no}: time must be finite')
        row = dict(row, dataset=dataset, time=time)
        rows.append(row)
    if not rows:
        raise ValueError(f'{path}: no samples for dataset {dataset!r}')
    for field in ('sample_id', 'bam'):
        if len({r[field] for r in rows}) != len(rows):
            raise ValueError(f'{path}: duplicate {field} within dataset {dataset!r}')
    cells = [(r['condition'], r['time'], r['replicate']) for r in rows]
    if len(set(cells)) != len(cells):
        raise ValueError(f'{path}: duplicate condition/time/replicate; give independent experiments distinct replicate IDs')
    all_bams = groups[0] + groups[1]
    if len(set(all_bams)) != len(all_bams):
        raise ValueError('BAM identifiers duplicated across SAMPLE_1 and SAMPLE_2')
    metadata_bams = {r['bam'] for r in rows}
    if metadata_bams != set(all_bams):
        raise ValueError('metadata and BAM order must map exactly; '
                         f'missing metadata={sorted(set(all_bams)-metadata_bams)[:5]}, '
                         f'unmapped metadata={sorted(metadata_bams-set(all_bams))[:5]}; no basename matching is performed')
    locations = {bam: {'group': gi + 1, 'vector_index': i, 'bam': bam}
                 for gi, group in enumerate(groups) for i, bam in enumerate(group)}
    for row in rows:
        row['source_position'] = locations[row['bam']]
    return sorted(rows, key=lambda r: r['sample_id'])


def validate_record(record):
    """Validate one canonical record. Raise ValueError; return None on success.

    Extra fields are permitted. Identities and aligned arrays are mandatory;
    counts must be actual JSON integers (not rounded floats or booleans).
    """
    if not isinstance(record, dict):
        raise ValueError('record must be a JSON object')
    for key in ('dataset', 'condition', 'mode', 'event_type', 'event_id', 'profile_id', 'record_key'):
        _nonempty(record.get(key), key)
    if record['mode'] not in ('JC', 'JCEC') or record['event_type'] not in EVENT_TYPES:
        raise ValueError('invalid count mode or event type')
    prefix = record['event_type'] + ':'
    if not record['event_id'].startswith(prefix) or not record['event_id'][len(prefix):]:
        raise ValueError('event_id must have TYPE:raw_ID format')
    expected_key = [record['dataset'], record['condition'], record['mode'], record['event_type'], record['event_id'][len(prefix):]]
    try:
        parsed_key = json.loads(record['record_key'])
    except (TypeError, ValueError) as exc:
        raise ValueError('record_key must encode the identity tuple as JSON') from exc
    if parsed_key != expected_key:
        raise ValueError('record_key does not match record identity')
    names = ('sample_ids', 'times', 'reps', 'I', 'S', 'N')
    if any(not isinstance(record.get(k), list) for k in names):
        raise ValueError('aligned record arrays must be lists')
    n = len(record['sample_ids'])
    if n == 0 or any(len(record[k]) != n for k in names):
        raise ValueError('record arrays have empty or inconsistent lengths')
    for key in ('sample_ids', 'reps'):
        for value in record[key]:
            _nonempty(value, key)
    if len(set(record['sample_ids'])) != n:
        raise ValueError('duplicate sample_id in record')
    for value in record['times']:
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
            raise ValueError('times must be finite numbers')
    if len(set(zip(record['times'], record['reps']))) != n:
        raise ValueError('duplicate time/replicate cell within record')
    for key in ('I', 'S', 'N'):
        if any(type(v) is not int or v < 0 or v > 2**53 - 1 for v in record[key]):
            raise ValueError(f'{key} must contain nonnegative exact-range integers')
    if any(i + s != total for i, s, total in zip(record['I'], record['S'], record['N'])):
        raise ValueError('N must equal I + S')
    for key in ('LI', 'LS'):
        value = record.get(key)
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
            raise ValueError(f'{key} must be a finite positive effective length')
    if not isinstance(record.get('annotation'), dict) or record['annotation'].get('ID') != expected_key[-1]:
        raise ValueError('annotation ID does not match raw event ID')
    for key in ('gene_id', 'gene_symbol'):
        if not isinstance(record.get(key), str):
            raise ValueError(f'{key} must be a string (empty allowed)')
    if 'source_positions' in record:
        positions = record['source_positions']
        if not isinstance(positions, list) or len(positions) != n:
            raise ValueError('source_positions must align with samples')
        for position in positions:
            if not isinstance(position, dict) or type(position.get('group')) is not int or position.get('group') not in (1, 2) or type(position.get('vector_index')) is not int or position['vector_index'] < 0:
                raise ValueError('invalid source vector position')
            _nonempty(position.get('bam'), 'source_positions.bam')
        if len({(p['group'], p['vector_index']) for p in positions}) != n or len({p['bam'] for p in positions}) != n:
            raise ValueError('duplicate source vector position or BAM within record')


def iter_records(path):
    """Yield validated records from a JSONL file; reject duplicate identities."""
    seen = set()
    with Path(path).open(encoding='utf-8') as stream:
        for line_no, line in enumerate(stream, 1):
            try:
                record = json.loads(line)
                validate_record(record)
                # Canonicalize equivalent whitespace in an externally supplied key.
                identity = tuple(json.loads(record['record_key']))
                if identity in seen:
                    raise ValueError('duplicate record identity')
                seen.add(identity)
            except (ValueError, TypeError) as exc:
                raise ValueError(f'{path}:{line_no}: {exc}') from exc
            yield record


def _gene(value):
    return value[1:-1] if len(value) >= 2 and value[0] == value[-1] == '"' else value


def prepare_rmats(counts_dir, samples_path, bam_list1, bam_list2, out_dir,
                  dataset, mode='JC', event_types=None):
    """Atomically prepare all event/condition records from counts-only rMATS.

    SAMPLE_2 needs an explicit bam_list2 when nonempty; an absent group is
    represented by empty strings, never by implicit sample order assumptions.
    Existing output directories are never overwritten.
    """
    _nonempty(dataset, 'dataset')
    if mode not in ('JC', 'JCEC'):
        raise ValueError('mode must be JC or JCEC; modes cannot be pooled')
    types = list(EVENT_TYPES if event_types is None else event_types)
    if not types or len(set(types)) != len(types) or set(types) - set(EVENT_TYPES):
        raise ValueError('event_types must be a nonempty distinct subset of SE,A5SS,A3SS,MXE,RI')
    types = sorted(types)
    out = Path(out_dir).absolute()
    if out.exists():
        raise FileExistsError(f'output already exists: {out}')
    counts_dir = Path(counts_dir)
    groups = [_bam_order(bam_list1), _bam_order(bam_list2)]
    if not groups[0]:
        raise ValueError('bam_list1 is required')
    samples = _samples(samples_path, dataset, groups)
    conditions = defaultdict(list)
    for sample in samples:
        conditions[sample['condition']].append(sample)
    sources = [Path(samples_path), Path(bam_list1)]
    if bam_list2 is not None:
        sources.append(Path(bam_list2))
    records = []
    annotations = []
    counts_by_type = {}
    for event_type in types:
        ann_path = counts_dir / f'fromGTF.{event_type}.txt'
        raw_path = counts_dir / f'{mode}.raw.input.{event_type}.txt'
        sources.extend([ann_path, raw_path])
        ann = {}
        for line_no, row in _table(ann_path, ('ID', 'GeneID', 'geneSymbol')):
            event = _nonempty(row['ID'], f'{ann_path}:{line_no}:ID')
            if event in ann:
                raise ValueError(f'{ann_path}: duplicate annotation ID {event!r}')
            ann[event] = row
        seen = set()
        for line_no, row in _table(raw_path, _COUNT_COLUMNS):
            event = _nonempty(row['ID'], f'{raw_path}:{line_no}:ID')
            if event in seen:
                raise ValueError(f'{raw_path}: duplicate count ID {event!r}')
            if event not in ann:
                raise ValueError(f'{raw_path}: annotation missing for count ID {event!r}')
            seen.add(event)
            lengths = [_integer(row[k], f'{raw_path}:{line_no}:{k}') for k in ('IncFormLen', 'SkipFormLen')]
            if min(lengths) <= 0:
                raise ValueError(f'{raw_path}:{line_no}: effective lengths must be positive')
            vectors = []
            for gi, group in enumerate(groups, 1):
                pair = []
                for form in ('IJC', 'SJC'):
                    field = f'{form}_SAMPLE_{gi}'
                    value = row[field]
                    if value and not group:
                        raise ValueError(f'{raw_path}:{line_no}: {field} nonempty but no BAM group mapping provided')
                    vector = [] if value == '' else [_integer(v, f'{raw_path}:{line_no}:{field}') for v in value.split(',')]
                    if len(vector) != len(group):
                        raise ValueError(f'{raw_path}:{line_no}: {field} length {len(vector)} != BAM group length {len(group)}')
                    pair.append(vector)
                vectors.append(pair)
            for condition in sorted(conditions):
                subset = conditions[condition]
                positions = [s['source_position'] for s in subset]
                inclusion = [vectors[p['group']-1][0][p['vector_index']] for p in positions]
                skipping = [vectors[p['group']-1][1][p['vector_index']] for p in positions]
                record = {
                    'dataset': dataset, 'condition': condition, 'mode': mode,
                    'event_type': event_type, 'event_id': event_type + ':' + event,
                    'profile_id': f'{dataset}/{condition}/{mode}/{event_type}:{event}',
                    'record_key': json.dumps([dataset, condition, mode, event_type, event], separators=(',', ':'), ensure_ascii=False),
                    'sample_ids': [s['sample_id'] for s in subset],
                    'times': [s['time'] for s in subset], 'reps': [s['replicate'] for s in subset],
                    'I': inclusion, 'S': skipping, 'N': [i+s for i, s in zip(inclusion, skipping)],
                    'LI': lengths[0], 'LS': lengths[1],
                    'gene_id': _gene(ann[event]['GeneID']), 'gene_symbol': _gene(ann[event]['geneSymbol']),
                    'annotation': ann[event], 'source_positions': positions,
                }
                validate_record(record)
                records.append(record)
            annotations.append({'event_type': event_type, 'event_id': event_type + ':' + event, 'annotation': ann[event]})
        # Extra annotation rows are retained in the source manifest count but do
        # not generate fabricated count records.
        counts_by_type[event_type] = {'count_events': len(seen), 'annotation_events': len(ann), 'annotation_without_counts': len(set(ann)-seen)}
    records.sort(key=lambda r: (r['condition'], r['event_type'], r['annotation']['ID']))
    annotations.sort(key=lambda a: (a['event_type'], a['annotation']['ID']))
    manifest = {
        'schema_version': '1.0', 'created_utc': datetime.now(timezone.utc).isoformat(),
        'dataset': dataset, 'mode': mode, 'event_types': types,
        'sample_count': len(samples), 'condition_count': len(conditions),
        'record_count': len(records), 'events': counts_by_type,
        'bam_group_sizes': [len(g) for g in groups],
        'sources': [{'path': str(p.absolute()), 'sha256': _sha(p), 'bytes': p.stat().st_size} for p in sources],
        'mapping_policy': 'exact literal BAM identifier; no basename fallback',
        'sort_policy': 'records condition/event_type/raw_ID; samples sample_id lexicographic',
        'filtering': 'none; all count events x conditions retained',
        'outputs': {},
    }
    out.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix='.' + out.name + '.tmp-', dir=out.parent))
    try:
        for name, rows in (('records.jsonl', records), ('samples.jsonl', samples), ('annotations.jsonl', annotations)):
            destination = temporary / name
            with destination.open('w', encoding='utf-8', newline='\n') as stream:
                for row in rows:
                    stream.write(json.dumps(row, ensure_ascii=False, allow_nan=False, separators=(',', ':')) + '\n')
            manifest['outputs'][name] = {'sha256': _sha(destination), 'bytes': destination.stat().st_size}
        (temporary / 'prepare_manifest.json').write_text(json.dumps(manifest, ensure_ascii=False, allow_nan=False, indent=2) + '\n', encoding='utf-8')
        if out.exists():
            raise FileExistsError(f'output appeared during preparation: {out}')
        temporary.rename(out)
    except BaseException:
        shutil.rmtree(temporary, ignore_errors=True)
        raise
    return manifest
