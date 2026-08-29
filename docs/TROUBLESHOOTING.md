# ArchLink Troubleshooting

This document records runtime failures encountered during ArchLink validation
and the recovery commands for the released repository.

## Required databases

ArchLink uses two different CheckM installations:

- **CheckM1** is used internally by UniItem to score candidate Leiden
  partitions. Its database must contain `pfam/Pfam-A.hmm.dat`.
- **CheckM2** is used after secondary binning to assess the final bins.

Configure both paths in the YAML file. The CheckM1 path must be a local path
on the machine running the analysis:

```yaml
common:
  path:
    checkm1_data_path: "/absolute/path/to/CheckM1"
    checkm2_path: "/absolute/path/to/checkm2_env"
```

The same CheckM1 path can be supplied without editing YAML:

```bash
export CHECKM_DATA_PATH=/absolute/path/to/CheckM1
```

Before UniItem starts, ArchLink verifies `Pfam-A.hmm.dat` and synchronizes
the CheckM1 root with `checkm data setRoot`. An incomplete database now stops
the run with an actionable error instead of producing empty quality tables.

## FragGeneScan

The bundled C source must be compiled once on each Linux installation:

```bash
bash scripts/build_fraggenescan.sh
```

The main entry point also checks for the executable and attempts the same
build. The marker extraction script is called explicitly through Perl, so the
Unix executable bit on `test_getmarker_2quarter.pl` is not required.

## Marker seed and clustering compatibility

ArchLink calls the marker script as `perl <script> ...` and checks that the
expected `.seed` file exists and is non-empty. The released environment pins
NumPy `1.23.5`, SciPy `1.10.1`, scikit-learn `1.1.2`, and joblib `1.2.0` to
match the bundled random-forest models. KMeans uses the supported
`algorithm="lloyd"` value.

## Resume after an interrupted run

After contrastive learning and primary clustering have completed, resume the
secondary binning and linking stages with:

```bash
bash scripts/resume_archlink_linking.sh \
  --config /absolute/path/to/configuration.yaml \
  --start-at auto
```

Available starting points are `binning`, `cosine`, and `connect`. `auto`
reuses complete, non-empty outputs. The resume entry point is a standalone
Python file because multiprocessing with the `spawn` method cannot reliably
re-import a Python program supplied through a shell here-document.

The script requires:

- `cluster_res/estimate_res.txt` from primary clustering;
- `bam.graph`;
- a valid CheckM1 database;
- CheckM2 output and secondary bin FASTA files before cosine/connect stages.

An interrupted process with exit code `143` was terminated by `SIGTERM`; it
does not indicate a model or algorithm failure. Re-run the resume command
after the scheduler or operating-system interruption has ended.

## Missing `extracted_edges.npz`

The secondary RF-enhanced Leiden input must contain:

```text
extracted_edges.npz
namelist.txt
length_weight.txt
```

ArchLink reuses these files when they are complete. If they are absent or
empty, it regenerates the initial edges and RF-enhanced edges. The program
now raises an explicit error if an expected edge file cannot be loaded.

## MAG export semantics

The paired comparison uses the same bin IDs before and after linking:

| Set | Source | Filter |
| --- | --- | --- |
| `before_link` | `{output}/binning/bins/` | CheckM2 contamination < 10% |
| `after_link` | `{output}/linking/connect/` | Same bin IDs as `before_link` |

Export the paired sets for downstream QC with:

```bash
python scripts/export_archlink_mags.py \
  --output-path /absolute/path/to/output
```

The command writes `mags/export/EXPORT.txt`, which records the source,
contamination threshold, counts, and any missing after-link IDs.

Validate that a run produced actual MAG files:

```bash
python scripts/output_validation.py \
  --output-path /absolute/path/to/output
```

The full pipeline exits non-zero when either before-link or after-link output
contains zero non-empty FASTA files. A zero-MAG run must not be marked as a
successful downstream stage.
