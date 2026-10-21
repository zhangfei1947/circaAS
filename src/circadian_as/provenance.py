"""Portable hashes, deterministic event seeds and strict JSON serialization."""
import hashlib
import json
import math
from pathlib import Path


def sha256_file(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def event_digest(record_key):
    return hashlib.sha256(record_key.encode('utf-8')).hexdigest()


def event_seed(seed, record_key):
    payload = json.dumps([seed, record_key], ensure_ascii=False, separators=(',', ':'))
    return int.from_bytes(hashlib.sha256(payload.encode('utf-8')).digest()[:8], 'big')


def safe(value):
    if isinstance(value, dict):
        return {str(k): safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [safe(v) for v in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if hasattr(value, 'item'):
        return safe(value.item())
    return value


def write_json(path, obj):
    path = Path(path)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(safe(obj), ensure_ascii=False, indent=2, allow_nan=False) + '\n')
    temporary.replace(path)
