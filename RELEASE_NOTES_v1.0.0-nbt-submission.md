# ArchLink v1.0.0-nbt-submission

Fixed software snapshot for the manuscript:

> ArchLink: multi-evidence reconstruction of genomes and genomic context from short-read metagenomes

## Release identity

- Algorithmic manuscript baseline: `af2231f` (`Improve ArchLink recovery and runtime validation`)
- Release tag: `v1.0.0-nbt-submission`
- Manuscript software snapshot date: 2026-08-31
- Reproducibility packaging completed: 2026-09-01
- Repository: `Sharkbio/ArchLink`
- Processed-results archive: Zenodo DOI `10.5281/zenodo.20487052`

The release commit adds documentation, synthetic toy inputs, expected outputs,
figure-generation entry points and the Source Data workbook. Core algorithmic
files remain those of the `af2231f` baseline.

## Bundled model artifacts

The released Transformer checkpoint, random-forest models and feature schemas
are in `save_models/`. Exact SHA-256 values are recorded in
`save_models/SHA256SUMS.txt`.

## Configuration and runtime

- Full configuration template: `configuration.yaml`
- Reviewer full-workflow toy configuration: `example/config.toy.yaml`
- Portable smoke test: `bash example/run_toy.sh`
- Conda environment: `environment.yml`
- Python: 3.10
- NumPy: 1.23.5
- pandas: 1.5.3
- SciPy: 1.10.1
- scikit-learn: 1.1.2
- joblib: 1.2.0
- matplotlib: 3.7.5
- openpyxl: 3.1.5

The pinned scikit-learn version is required for compatibility with the bundled
random-forest pickle files. PyTorch and external bioinformatics executables are
resolved through `environment.yml`; CheckM1 and CheckM2 databases are not
redistributed.

## Source Data and figures

- Source Data workbook: `source_data/Source_Data.xlsx`
- Workbook SHA-256: `6ae85c320a615795cae0209eb089e50c5330a52fa4ad26e967ea7f73f217036d`
- Central provenance map: `source_data/Figure_Panel_Data_Source_Map_v3.csv`
- One-command quantitative build: `python figures/make_all_figures.py`
- Fixed figure seed: `20260831`
- Default resolution: 300 dpi

## Known limitations

- The full ArchLink pipeline targets Linux and calls external executables such
  as HMMER, bedtools, samtools, prodigal, CheckM and CheckM2.
- The portable toy smoke test validates paired-end evidence, abstention and
  released FASTA joining without retraining models or running CheckM. The
  separate end-to-end tier requires the complete external runtime.
- Fig. 1 is a manually assembled conceptual schematic.
- The Source Data for Fig. 4a contains rank-level counts, not the MAG-level
  transition table needed to recreate the original alluvial geometry; the
  provided script renders the quantitative rank-count summary.
- Fig. 6b and Fig. 6c depend on antiSMASH/GenBank context and external
  AlphaFold/PyMOL rendering; final multi-panel assembly remains manual.
- Only the Waste Water dataset has independent Oxford Nanopore long-read
  junction validation. Deep Marine `connect_00055` is a computational candidate
  architecture and is not experimental functional validation.

## Verification commands

```bash
python scripts/repository_audit.py
bash example/run_toy.sh
python figures/make_all_figures.py \
  --source-data source_data/Source_Data.xlsx \
  --output figures/output
git diff --check
```
