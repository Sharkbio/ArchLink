# ArchLink

ArchLink is a graph-based metagenomic binning and context-aware scaffolding workflow for improving microbial genome reconstruction from short-read assemblies.

This repository contains the source code, configuration templates, bundled helper scripts, pretrained model artifacts, and helper binaries required to run the released ArchLink pipeline on Linux.

## Repository layout

- `archlink.py`: top-level pipeline entry point
- `configuration.yaml`: full configuration template
- `environment.yml`: conda environment specification
- `example/`: minimal reviewer-facing run layout and configuration template
- `benchmarks/`: public benchmark/evaluation entry points
- `contrastive_learning/`: representation learning and initial binning
- `generate01/`, `binning02/`: graph construction and bin refinement
- `generate_cos03/`, `connect04/`: context-aware linking and scaffolding
- `save_models/`: pretrained Transformer checkpoint, random-forest models, and helper binaries
- `FragGeneScan-master/`: bundled FragGeneScan source, Perl wrapper, and training files
- `scripts/repository_audit.py`: repository completeness check for public release

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
conda activate archlink
```

The environment file includes the core packaged dependencies used directly by the repository, including:

- PyTorch
- `make` and a C compiler for building the bundled FragGeneScan executable
- HMMER
- bedtools
- samtools
- prodigal
- CheckM2
- Perl runtime

## Quick repository check

Reviewers can confirm that the public repository contains the expected source files and bundled artifacts with:

```bash
python scripts/repository_audit.py
python archlink.py --help
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
- `common.path.LD_LIBRARY_PATH`

The template uses repository-relative defaults so that the file can be versioned safely.

## Running ArchLink

Run the full pipeline with:

```bash
python archlink.py --config configuration.yaml
```

## Minimal example status

A reviewer-facing minimal example layout is documented in [example/README.md](./example/README.md).

The repository already includes:

- the executable pipeline
- configuration templates
- pretrained model weights
- random-forest model files
- benchmark/evaluation entry points

The repository does not yet bundle a redistributable toy input dataset consisting of a small contig FASTA and matching sorted BAM files. Until such a toy dataset is added, the minimal example section documents the expected layout and launch command, but cannot serve as a fully self-contained end-to-end demo.

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

## Reproducibility notes

- The repository contains the released source code used by the ArchLink workflow.
- The pretrained Transformer checkpoint, random-forest model files, and helper binaries required by the released linking stages are bundled in `save_models/`.
- The main configuration template, a minimal example configuration, and the environment specification are versioned in the repository.
- `contrastive_learning/train_CLmodel.py` sets explicit PyTorch random seeds for model training code paths.

## Current gaps before manuscript submission

This repository is substantially closer to reviewer-ready than an incomplete code drop, but two items still merit explicit completion before submission:

1. create a fixed GitHub release/tag for the submitted software version
2. add a small redistributable toy dataset so the minimal example becomes fully executable from the public repository alone

If manuscript figures and source-data generation scripts are intended to be part of the code release, they should also be added as explicit directories such as `figures/` and `source_data/` with one documented entry point per figure panel or table.

## Citation

If you use ArchLink, please cite the accompanying manuscript.
