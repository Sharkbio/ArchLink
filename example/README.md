# Reviewer-facing toy example

The repository includes a deterministic, redistributable FASTA/BAM fixture that
demonstrates ArchLink's paired-end evidence handling, explicit abstention on an
ambiguous competition, and released FASTA joining code without requiring the
large CheckM databases or model-training stages.

## Bundled layout

```text
example/
  data/
    toy_contigs.fasta
    toy_reads.sorted.bam
    toy_reads.sorted.bam.bai
  precomputed/
    contig_bins.tsv
  expected/
    expected_bins.tsv
    expected_candidates.tsv
    expected_links.tsv
    expected_unconnected_ends.tsv
    expected_scaffolds.fasta
  config.toy.yaml
  create_toy_data.py
  run_toy.py
  run_toy.sh
  run_toy.ps1
  run_toy_end_to_end.sh
```

The eight contigs represent two simplified genomes.  Twelve paired reads support
the accepted `contig01:R -- contig02:L` connection.  Two five-read candidates
compete at `contig03:R`, so neither is accepted.  Contigs 7 and 8 have no
cross-contig evidence and remain unjoined.

## Tier 1: portable smoke test

This test runs in seconds and does not invoke CheckM, CheckM2, model training or
GPU code.  It extracts candidate evidence from the bundled BAM, applies the
documented reciprocal-best rule with a minimum of eight spanning read pairs,
and uses `connect04.make_fa2` from the released ArchLink code to create the
scaffold FASTA.

Linux/macOS:

```bash
bash example/run_toy.sh
```

PowerShell:

```powershell
./example/run_toy.ps1
```

Equivalent explicit commands:

```bash
python example/create_toy_data.py --verify-only
python example/run_toy.py
python scripts/check_toy_output.py \
  --observed example/output/toy_run \
  --expected example/expected \
  --data example/data
```

The validator checks input presence and BAM/BAI structure, two-bin membership,
candidate and final-join tables, the ambiguous unconnected contig end, and the
SHA-256 checksum of the complete scaffold FASTA.

## Tier 2: full end-to-end workflow

`config.toy.yaml` supplies the same FASTA/BAM paths to the complete ArchLink
entry point.  Edit the CheckM1/CheckM2 and `LD_LIBRARY_PATH` placeholders first,
activate the pinned conda environment, then run:

```bash
bash example/run_toy_end_to_end.sh
```

This tier exercises the full binning, model and linking workflow and therefore
requires the external databases and Linux executables documented in the main
README.  The portable Tier 1 smoke test is the reviewer default.
