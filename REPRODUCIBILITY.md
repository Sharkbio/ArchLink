# Reproducibility Notes

## Fixed manuscript release

The ArchLink algorithm used for the manuscript is anchored at Git commit
`af2231f` (`Improve ArchLink recovery and runtime validation`). The public
submission archive is tagged `v1.0.0-nbt-submission`. The release commit adds
reviewer-facing toy data, figure entry points, Source Data and documentation;
it does not alter the algorithmic files present at `af2231f`.

The processed-results archive is available at Zenodo DOI
`10.5281/zenodo.20487052`. The versioned Source Data workbook in this repository
has SHA-256
`6ae85c320a615795cae0209eb089e50c5330a52fa4ad26e967ea7f73f217036d`.

This repository includes the configuration template (`configuration.yaml`), conda environment specification (`environment.yml`), and bundled model artifacts in `save_models/` that are required to run the released ArchLink pipeline.

The ArchLink runtime environment is named `ly_archlink` and pins NumPy `1.23.5`, SciPy `1.10.1`,
scikit-learn `1.1.2`, and joblib `1.2.0`. These versions are intentional: the bundled random-forest
pickle files were trained with scikit-learn `1.1.2` and may fail to load with newer scikit-learn
tree representations. CheckM2 can be invoked from a separate environment, for example
`ly_checkm2`, through the configured `checkm2_bin` or `checkm2_path`.
ArchLink's UniItem stage separately requires a CheckM1 database containing
`pfam/Pfam-A.hmm.dat`; configure it with `common.path.checkm1_data_path` or
the `CHECKM_DATA_PATH` environment variable.

Randomness control is partially fixed in the released training and clustering code paths. Examples include `torch.manual_seed(1)` and deterministic cuDNN settings in `contrastive_learning/train_CLmodel.py`, `torch.manual_seed(0)` in `contrastive_learning/simclr.py`, and `random_state=7` in the seed-kmeans implementation in `contrastive_learning/cluster.py`.

Because the project was developed as a multi-stage workflow over time, the repository does not currently provide a single centralized seed manifest covering every downstream benchmark and figure-generation run used in the manuscript.

## Reviewer smoke test

`bash example/run_toy.sh` validates the bundled FASTA/BAM/BAI fixture and checks
the observed two-bin mapping, candidate links, accepted join, ambiguous
unconnected end and scaffold checksum against `example/expected/`. This tier is
independent of CheckM databases and GPU training. The optional
`example/run_toy_end_to_end.sh` uses `example/config.toy.yaml` and requires the
full external runtime described in the main README.

## Figure entry points

`python figures/make_all_figures.py` reads explicit sheets from
`source_data/Source_Data.xlsx`, uses seed `20260831`, fixed dataset/tool order and
fixed colors, and writes `figure_build_manifest.json` with software versions and
the Source Data checksum. Fig. 1 and externally rendered/assembled Fig. 6
components are documented rather than falsely claimed as one-command plots.
