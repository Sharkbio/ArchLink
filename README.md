# ArchLink

ArchLink is a graph-based metagenomic binning and context-aware scaffolding workflow for improving microbial genome reconstruction from short-read assemblies.

This repository contains the source code, configuration template, bundled helper scripts, and the model/binary layout required to run the published ArchLink pipeline. The default execution target is Linux.

## Repository contents

- `archlink.py`: top-level pipeline entry point
- `configuration.yaml`: editable configuration template
- `environment.yml`: conda environment specification
- `contrastive_learning/`: contrastive representation learning and initial binning
- `generate01/`, `binning02/`: graph construction and bin refinement
- `generate_cos03/`, `connect04/`: context-aware linking and scaffolding
- `save_models/`: pretrained model artifacts and helper binaries used by the linking stages
- `FragGeneScan-master/`: bundled FragGeneScan source and training files

## Requirements

ArchLink is developed for Linux with Python 3.10. Core software dependencies are listed in `environment.yml`.

Create the conda environment:

```bash
conda env create -f environment.yml
conda activate archlink
```

External executables used by the pipeline:

- `hmmsearch` from HMMER
- `prodigal`
- Perl runtime for `run_FragGeneScan.pl`
- `samtools`
- `checkm2`

The current pipeline expects the `checkm2` executable under `${checkm2_path}/bin/checkm2`, as configured in `configuration.yaml`.
FragGeneScan is bundled in this repository under `FragGeneScan-master/` and is invoked from that location by the pipeline.

## Model artifacts

The `save_models/` directory is used for pretrained random-forest models, helper binaries, and the Transformer checkpoint used in the linking stage.

Artifacts currently bundled in this repository layout:

- `generateG13`
- `matching`
- `bacteria_transformer2.pth`
- `best_random_forest_model_focus0_D_B2.pkl`
- `feature_columns_focus0_D_B2.pkl`
- `best_random_forest_model_gas_connect_COMB_A_weight1_A_weight23.pkl`
- `feature_columns_gas_connect_COMB_A_weight1_A_weight23.pkl`
- `best_random_forest_model_gas_connect_COMB_C1_cosine_C2_cosine3.pkl`
- `feature_columns_gas_connect_COMB_C1_cosine_C2_cosine3.pkl`

## Configuration

Edit `configuration.yaml` before use. At minimum, update these fields:

- `common.path.contig_file`
- `common.path.bam_file`
- `common.path.base_path`
- `common.path.ID`
- `common.path.checkm2_path`
- `common.path.LD_LIBRARY_PATH`

The template uses repository-relative defaults so it can be versioned safely. Runtime outputs are written under `{base_path}/{ID}`.

## Running ArchLink

```bash
python archlink.py --config configuration.yaml
```

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
- scaffolding outputs for bin-level connection results

## Reproducibility notes

- The repository includes the source code used for the ArchLink workflow and a versioned environment file.
- Random-forest model artifacts, helper binaries, and the Transformer checkpoint currently used by the released pipeline are included in `save_models/`.
- The configuration template and environment file are versioned with the repository to support reproducible setup.

## Citation

If you use ArchLink, please cite the accompanying manuscript.
