# Minimal example workflow

This directory documents the smallest reviewer-facing run layout for ArchLink.

## What is included

- `config.minimal.yaml`: a minimal configuration template for a single test run

## What is not included yet

The public repository does not currently bundle a redistributable toy dataset with:

- a small contig FASTA file
- a matching directory of coordinate-sorted BAM files

Without those two inputs, reviewers can inspect the full pipeline, validate bundled models, and run repository checks, but they cannot execute an end-to-end example from this directory alone.

## Expected layout for a minimal example

```text
example/
  config.minimal.yaml
  data/
    mini_run/
      contigs.fasta
      bam/
        sample1.sorted.bam
        sample1.sorted.bam.bai
```

## Reviewer smoke checks

From the repository root:

```bash
python scripts/repository_audit.py
python archlink.py --help
```

## End-to-end command once toy input data are added

```bash
python archlink.py --config example/config.minimal.yaml
```
