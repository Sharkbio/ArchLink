#!/usr/bin/env python3
"""Validate toy smoke-test inputs and outputs against versioned expectations."""

from __future__ import annotations

import argparse
import csv
import hashlib
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from example.toy_bam import iter_bam_records, read_bam_header, validate_bai  # noqa: E402


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def compare_tsv(observed: Path, expected: Path) -> list[dict[str, str]]:
    observed_rows = read_tsv(observed)
    expected_rows = read_tsv(expected)
    if observed_rows != expected_rows:
        raise AssertionError(
            f"TSV differs: {observed.name}\n"
            f"observed={observed_rows!r}\nexpected={expected_rows!r}"
        )
    return observed_rows


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--observed", type=Path, required=True)
    parser.add_argument("--expected", type=Path, default=ROOT / "example" / "expected")
    parser.add_argument("--data", type=Path, default=ROOT / "example" / "data")
    args = parser.parse_args()

    fasta = args.data / "toy_contigs.fasta"
    bam = args.data / "toy_reads.sorted.bam"
    bai = bam.with_suffix(bam.suffix + ".bai")
    for path in (fasta, bam, bai):
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(f"Missing or empty toy input: {path}")
    references, lengths = read_bam_header(bam)
    if len(references) != 8 or lengths != [1000] * 8:
        raise AssertionError("Toy BAM must contain eight 1-kb contig references")
    validate_bai(bai, len(references))
    bam_records = list(iter_bam_records(bam))
    if len(bam_records) != 50 or len({record.query_name for record in bam_records}) != 25:
        raise AssertionError("Toy BAM read-pair count is incorrect")

    mapping = {
        "bins.tsv": "expected_bins.tsv",
        "candidate_links.tsv": "expected_candidates.tsv",
        "final_joins.tsv": "expected_links.tsv",
        "unconnected_ends.tsv": "expected_unconnected_ends.tsv",
    }
    compared = {
        name: compare_tsv(args.observed / name, args.expected / expected_name)
        for name, expected_name in mapping.items()
    }
    bin_count = len({row["bin_id"] for row in compared["bins.tsv"]})
    accepted = compared["final_joins.tsv"]
    abstained = [
        row for row in compared["candidate_links.tsv"] if row["status"] == "abstained"
    ]
    unconnected = compared["unconnected_ends.tsv"]
    if bin_count != 2:
        raise AssertionError(f"Expected 2 toy bins, found {bin_count}")
    if len(accepted) != 1:
        raise AssertionError(f"Expected 1 accepted join, found {len(accepted)}")
    if len(abstained) != 2 or not all("tied at contig03:R" in row["reason"] for row in abstained):
        raise AssertionError("Ambiguous contig03 competition was not abstained")
    if not any(row["contig_id"] == "contig03" and row["end"] == "R" for row in unconnected):
        raise AssertionError("The ambiguous contig03 right end must remain unconnected")

    observed_scaffolds = args.observed / "scaffolds.fasta"
    expected_scaffolds = args.expected / "expected_scaffolds.fasta"
    if sha256(observed_scaffolds) != sha256(expected_scaffolds):
        raise AssertionError("Scaffold FASTA checksum differs from expected output")

    print("Toy input/output validation passed.")
    print(f"Bins: {bin_count}")
    print(f"Candidate links: {len(compared['candidate_links.tsv'])}")
    print(f"Accepted joins: {len(accepted)}")
    print(f"Unconnected ends: {len(unconnected)}")
    print(f"Scaffold SHA-256: {sha256(observed_scaffolds)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
