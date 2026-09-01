# Synthetic toy input data

These files are fully synthetic and contain no biological or personal data.

- `toy_contigs.fasta`: eight deterministic 1-kb contigs representing two simplified genomes.
- `toy_reads.sorted.bam`: 25 deterministic paired-read alignments, coordinate sorted.
- `toy_reads.sorted.bam.bai`: BAI index for the BAM file.

Regenerate and validate them with:

```bash
python example/create_toy_data.py
python example/create_toy_data.py --verify-only
```

The fixture contains one decisive contig-end link, one tied two-way competition
that must be abstained, one within-contig read-pair control, and contigs with no
cross-contig evidence.
