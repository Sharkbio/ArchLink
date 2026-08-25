# Reproducibility Notes

This repository includes the configuration template (`configuration.yaml`), conda environment specification (`environment.yml`), and bundled model artifacts in `save_models/` that are required to run the released ArchLink pipeline.

The ArchLink runtime environment is named `ly_archlink` and pins NumPy `1.23.5`, SciPy `1.10.1`,
scikit-learn `1.1.2`, and joblib `1.2.0`. These versions are intentional: the bundled random-forest
pickle files were trained with scikit-learn `1.1.2` and may fail to load with newer scikit-learn
tree representations. CheckM2 can be invoked from a separate environment, for example
`ly_checkm2`, through the configured `checkm2_bin` or `checkm2_path`.

Randomness control is partially fixed in the released training and clustering code paths. Examples include `torch.manual_seed(1)` and deterministic cuDNN settings in `contrastive_learning/train_CLmodel.py`, `torch.manual_seed(0)` in `contrastive_learning/simclr.py`, and `random_state=7` in the seed-kmeans implementation in `contrastive_learning/cluster.py`.

Because the project was developed as a multi-stage workflow over time, the repository does not currently provide a single centralized seed manifest covering every downstream benchmark and figure-generation run used in the manuscript.
