# ArchLink

ArchLink is a graph-based metagenomic binning and context-aware scaffolding workflow for improving microbial genome reconstruction from short-read assemblies.

This repository contains the source code, configuration templates, bundled helper scripts, pretrained model artifacts, and helper binaries required to run the released ArchLink pipeline on Linux.

## Repository layout

- `archlink.py`: top-level pipeline entry point
- `configuration.yaml`: full configuration template
- `environment.yml`: conda environment specification
- `example/`: bundled synthetic FASTA/BAM toy example, expected outputs, and validators
- `figures/`: deterministic main-figure generation entry points
- `source_data/`: versioned Source Data workbook and panel provenance map
- `benchmarks/`: public benchmark/evaluation entry points
- `contrastive_learning/`: representation learning and initial binning
- `generate01/`, `binning02/`: graph construction and bin refinement
- `generate_cos03/`, `connect04/`: context-aware linking and scaffolding
- `save_models/`: pretrained Transformer checkpoint, random-forest models, and helper binaries
- `FragGeneScan-master/`: bundled FragGeneScan source, Perl wrapper, and training files
- `scripts/repository_audit.py`: repository completeness check for public release
- `scripts/resume_archlink_linking.py`: resume secondary binning and linking
- `scripts/export_archlink_mags.py`: export paired before/after MAG sets
- `docs/TROUBLESHOOTING.md`: runtime failures and recovery commands

## Included release artifacts

The following release artifacts are bundled in this repository:

- pretrained Transformer checkpoint:
  - `save_models/bacteria_transformer2.pth`
- random-forest model files and feature definitions:
  - `save_models/best_random_forest_model_focus0_D_B2.pkl`
  - `save_models/feature_columns_focus0_D_B2.pkl`
  - `save_models/best_random_forest_model_gas_connect_COMB_A_weight1_A_weight23.pkl`
  - `save_models/feature_columns_gas_connect_COMB_A_weight1_A_weight23.pkl`
  - `save_models/best_random_forest_model_gas_connect_COMB_C1_cosine_C2_cosine3.pkl`
  - `save_models/feature_columns_gas_connect_COMB_C1_cosine_C2_cosine3.pkl`
- helper binaries used by the linking stages:
  - `save_models/generateG13`
  - `save_models/matching`

## Software requirements

ArchLink is developed for Linux with Python 3.10.

Create the conda environment:

```bash
conda env create -f environment.yml
conda activate ly_archlink
```

The environment file includes the core packaged dependencies used directly by the repository, including:

- NumPy `1.23.5`
- SciPy `1.10.1`
- scikit-learn `1.1.2`
- joblib `1.2.0`
- pandas `1.5.3`
- matplotlib `3.7.5`
- openpyxl `3.1.5`
- PyTorch
- `make` and a C compiler for building the bundled FragGeneScan executable
- HMMER
- bedtools
- samtools
- prodigal
- CheckM2
- Perl runtime

The pinned NumPy/SciPy/scikit-learn/joblib versions are required for loading the bundled
random-forest pickle files, which were trained with scikit-learn `1.1.2`. Do not upgrade
scikit-learn in this environment unless the bundled models are retrained and revalidated.

CheckM2 may be run from a separate environment, such as `ly_checkm2`, when that installation
requires a different scikit-learn version. Configure its executable with
`common.path.checkm2_bin` or `common.path.checkm2_path` in `configuration.yaml`.
ArchLink's UniItem stage separately requires a CheckM1 database containing
`pfam/Pfam-A.hmm.dat`; configure it with `common.path.checkm1_data_path` or
the `CHECKM_DATA_PATH` environment variable.

## Quick repository check

Reviewers can confirm that the public repository contains the expected source files and bundled artifacts with:

```bash
python scripts/repository_audit.py
python archlink.py --help
python example/create_toy_data.py --verify-only
bash example/run_toy.sh
```

## External runtime expectations

ArchLink calls several external executables during the pipeline:

- `hmmsearch`
- `bedtools`
- `samtools`
- `prodigal`
- `checkm2`
- Perl for `FragGeneScan-master/run_FragGeneScan.pl`

FragGeneScan is bundled in this repository and is invoked from `FragGeneScan-master/`.

### Build FragGeneScan

The repository contains the FragGeneScan source code, Perl wrapper, and training files. The compiled
`FragGeneScan-master/FragGeneScan` executable is generated locally because it is platform-specific.
The build respects the `CC` environment variable; otherwise it uses the system C compiler `cc`.

After creating the conda environment, build it with:

```bash
bash scripts/build_fraggenescan.sh
```

Equivalent manual commands are:

```bash
make -C FragGeneScan-master clean
make -C FragGeneScan-master fgs
```

The main ArchLink entry point also checks for this executable and attempts to build it automatically
before the binning stage. If compilation fails, the error reports the required directory and command.

The marker extraction helper is a Perl script. ArchLink invokes it as
`perl contrastive_learning/auxiliary/test_getmarker_2quarter.pl`, so its Unix executable bit is not
required.

CheckM2 can be configured in either of two ways:

- set `common.path.checkm2_bin` to an executable name or full path
- or set `common.path.checkm2_path` to the root of an environment containing `bin/checkm2`

## Configuration

Edit `configuration.yaml` before running full analyses. At minimum, update:

- `common.path.contig_file`
- `common.path.bam_file`
- `common.path.base_path`
- `common.path.ID`
- `common.path.checkm2_bin` or `common.path.checkm2_path`
- `common.path.checkm1_data_path` or `CHECKM_DATA_PATH`
- `common.path.LD_LIBRARY_PATH`

The template uses repository-relative defaults so that the file can be versioned safely.

## Running ArchLink

Run the full pipeline with:

```bash
python archlink.py --config configuration.yaml
```

### Clustering modes

The default `full` mode preserves the manuscript-style parameter search. For large datasets or
development tests, use the faster `fast` mode:

```bash
python archlink.py --config configuration.yaml --clustering-mode fast
```

The initial clustering stage searches 540 Leiden parameter combinations in `full` mode and 8
representative combinations in `fast` mode. The post-random-forest clustering stage searches 50
combinations in `full` mode and 4 combinations in `fast` mode. Fast-mode results are intended for
pipeline testing and debugging; use `full` mode for final analyses and reported results.

The mode can also be stored in the YAML file under
`contrastive_learning.share_params.clustering_mode`. The release configuration uses `full`, while
`example/config.minimal.yaml` uses `fast`.

When rerunning with a different mode, use a new output `ID` or remove the previous clustering
outputs first. ArchLink reuses existing result files when their filenames already exist.

## Reviewer toy example

The repository bundles a synthetic eight-contig FASTA, coordinate-sorted BAM,
BAI index, two-bin mapping and deterministic expected outputs. The portable
smoke test demonstrates one accepted link, a tied competition that is explicitly
abstained, unconnected contig ends and the released scaffold FASTA builder:

```bash
bash example/run_toy.sh
```

See [example/README.md](./example/README.md) for the PowerShell command and the
optional full end-to-end tier requiring CheckM databases and Linux executables.

## Figure generation

The combined Source Data workbook is versioned at
`source_data/Source_Data.xlsx`. Regenerate all data-driven main-figure panels
and a machine-readable build manifest with:

```bash
python figures/make_all_figures.py \
  --source-data source_data/Source_Data.xlsx \
  --output figures/output
```

The exact sheet-to-panel map, fixed seed, sample/tool order, color map,
statistics and non-scripted panel limitations are documented in
[figures/README.md](./figures/README.md).

## Benchmark and evaluation entry points

Public benchmark-facing materials currently include:

- [benchmarks/README.md](./benchmarks/README.md)
- `benchmarks/linking_precision_recall.py`

This benchmark script evaluates predicted contig links against per-bin `ground_truth.txt` files and reports precision/recall summaries for the linking stage.

## Inputs

ArchLink requires:

- an assembled contig FASTA file
- a directory containing coordinate-sorted BAM files for the same sample or sample set

## Outputs

Main outputs are written under the configured `output_path` and include:

- contrastive embeddings
- initial and refined bins
- CheckM2 quality reports
- linking graphs
- bin-level scaffolding outputs

The pipeline validates that both the before-link and after-link MAG
collections contain non-empty FASTA files before reporting successful
completion. See [docs/TROUBLESHOOTING.md](./docs/TROUBLESHOOTING.md) for
resume and MAG export commands.

## Reproducibility notes

- The algorithmic manuscript baseline is Git commit `af2231f`.
- The fixed submission release is `v1.0.0-nbt-submission`; reproducibility-only
  additions after `af2231f` do not change the core ArchLink algorithm.
- The pretrained Transformer checkpoint, random-forest model files, and helper binaries required by the released linking stages are bundled in `save_models/`.
- The main configuration template, a minimal example configuration, and the environment specification are versioned in the repository.
- `contrastive_learning/train_CLmodel.py` sets explicit PyTorch random seeds for model training code paths.
- Model and Source Data SHA-256 checksums are versioned in
  `save_models/SHA256SUMS.txt` and `source_data/SHA256SUMS.txt`.

The complete submission snapshot and known limitations are listed in
[RELEASE_NOTES_v1.0.0-nbt-submission.md](./RELEASE_NOTES_v1.0.0-nbt-submission.md).

## Citation

If you use ArchLink, please cite the accompanying manuscript.
