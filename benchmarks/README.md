# Benchmark and evaluation notes

This directory collects public-facing evaluation entry points for ArchLink.

## Current contents

- `linking_precision_recall.py`: evaluates predicted contig links against per-bin `ground_truth.txt` files

## Expected benchmark input layout

The evaluator expects a directory containing one subdirectory per bin. Each bin directory should contain:

- `ground_truth.txt`
- `connect_dir.r`

Example:

```text
benchmark_case/
  bin_001/
    ground_truth.txt
    connect_dir.r
  bin_002/
    ground_truth.txt
    connect_dir.r
```

## Run the evaluator

```bash
python benchmarks/linking_precision_recall.py benchmark_case
```

## Scope

This repository currently exposes the link-level evaluation script used for benchmarking scaffolding predictions. If the manuscript reports additional benchmark tables or figure-generation workflows, those should be added here as explicit scripts and documented in the main README before submission.
