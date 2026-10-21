# rMATS input and canonical event records

`prepare_rmats(counts_dir, samples_path, bam_list1, bam_list2, out_dir, dataset,
mode='JC', event_types=None)` imports counts-only rMATS files. The output directory
must not exist. Preparation validates all inputs before publishing a complete
output directory; it never replaces an existing preparation. `event_types=None`
requires all five types: SE, A5SS, A3SS, MXE, RI. An explicit subset is permitted.

For each selected type, supply `{mode}.raw.input.{type}.txt` and
`fromGTF.{type}.txt`. Only one mode is imported per call. JC is the default; JCEC
must be explicit. These two evidence modes are never summed or pooled.

Tables are **literal tab-separated text**, not quote-aware CSV. Quotes in rMATS
annotation values are preserved even if unmatched. Annotation IDs remain strings;
`007` and `7` are different IDs. Count rows require ID, IJC_SAMPLE_1,
SJC_SAMPLE_1, IJC_SAMPLE_2, SJC_SAMPLE_2, IncFormLen, and SkipFormLen. Counts must
be nonnegative integers and effective lengths must be positive integers. Empty
vectors are accepted only for an absent group. Every count ID must have a unique
matching annotation row. Annotation rows without counts are counted in the
manifest but do not produce fabricated observations.

## Samples and BAM ordering

Samples TSV requires these columns:

```text
sample_id	condition	time	replicate	bam
WT_t0_r1	WT	0	r1	/data/WT_t0_r1.bam
```

Time is a finite numeric value in hours on one shared origin for each condition.
Replicate labels identify experiments spanning the time series (e.g. r1 and r2).
The importer does not infer the experimental structure from names. A `dataset`
column, when present, is filtered by exact equality to the requested dataset
before validation of that dataset's rows.

The actual rMATS BAM order files are mandatory. Supply the original SAMPLE_1
comma-separated BAM list as `bam_list1`. Supply `bam_list2` for a nonempty
SAMPLE_2, or Python `None` if SAMPLE_2 is empty. Comma-separated lists may also
have line breaks between entries. Each metadata BAM must match exactly one
literal identifier in these lists, including relative path components. The
importer does not resolve paths or use basename matching: BAM files need not be
accessible locally. Each ordered BAM must have exactly one metadata row.
Metadata row ordering has no effect on prepared records.

Duplicate sample IDs, duplicate BAM identifiers (including across groups), and
duplicate condition/time/replicate cells are rejected. Independent experiments
must have distinct replicate labels. Repeated identifiers are never silently
collapsed. Missing or extra metadata rows and count-vector length mismatches
cause preparation to fail.

## Prepared outputs

- `records.jsonl`: one record per count event and condition, including all-zero
  and otherwise insufficient-coverage records. No coverage filtering occurs.
- `samples.jsonl`: selected metadata, normalized finite time values and exact
  source group/vector positions, sorted by sample_id.
- `annotations.jsonl`: literal annotations for count events and their event types.
- `prepare_manifest.json`: schema version, dataset, mode, dimensions, BAM group
  sizes, per-type counts, source paths and SHA-256 hashes, output hashes, sorting
  and mapping policies. This is the preparation provenance; fitting adds its own
  run provenance.

Public `iter_records(path)` yields validated JSON objects from a JSONL file and
rejects duplicate identities. `validate_record(record)` returns `None` or raises
`ValueError`. Extra fields are permitted for future extensions.

Each record contains:

| Field | Meaning |
|---|---|
| `dataset`, `condition`, `mode`, `event_type` | Explicit analysis identity |
| `event_id` | Type prefix plus literal source ID, e.g. `RI:007` |
| `profile_id` | Human-readable `dataset/condition/mode/type:ID`; display only |
| `record_key` | JSON-encoded tuple `[dataset,condition,mode,event_type,raw_ID]`; unambiguous machine identity |
| `sample_ids`, `times`, `reps` | Aligned lists, sorted lexicographically by sample_id |
| `I`, `S`, `N` | Aligned integer count lists, with `N=I+S` |
| `LI`, `LS` | Positive scalar effective inclusion/skipping lengths |
| `gene_id`, `gene_symbol` | Annotation strings with only matched enclosing double quotes removed |
| `annotation` | Complete literal source annotation row, including raw `ID` |
| `source_positions` | Aligned objects with `group` (1 or 2), zero-based `vector_index`, and exact `bam` |

The model's count-probability offset is `log(LI/LS)`. PSI is based on effective
length adjusted inclusion/skipping abundance; it is not simply `I/N` when
lengths differ. Zero informative coverage supplies no PSI observation and is
never replaced with PSI zero by preparation.

Records are sorted by condition, event type and literal raw ID. Output content
is invariant to metadata row reordering, although the source metadata file hash
correctly changes when its bytes change. Preparation records timestamps, so
manifests are not intended to be byte-identical across runs.
