# Reproducibility Notes

This repository includes the configuration template (`configuration.yaml`), conda environment specification (`environment.yml`), and bundled model artifacts in `save_models/` that are required to run the released ArchLink pipeline.

Randomness control is partially fixed in the released training and clustering code paths. Examples include `torch.manual_seed(1)` and deterministic cuDNN settings in `contrastive_learning/train_CLmodel.py`, `torch.manual_seed(0)` in `contrastive_learning/simclr.py`, and `random_state=7` in the seed-kmeans implementation in `contrastive_learning/cluster.py`.

Because the project was developed as a multi-stage workflow over time, the repository does not currently provide a single centralized seed manifest covering every downstream benchmark and figure-generation run used in the manuscript.
