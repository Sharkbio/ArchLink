#!/usr/bin/env python3
"""Generate the deterministic, redistributable ArchLink toy FASTA/BAM fixture."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

from toy_bam import AlignmentSpec, iter_bam_records, read_bam_header, validate_bai, write_bam_and_bai


CONTIG_LENGTH = 1000
CONTIG_NAMES = [f"contig{i:02d}" for i in range(1, 9)]


def deterministic_sequence(label: str, length: int, gc_percent: int) -> str:
    """Create stable non-repetitive sequence without relying on random module details."""

    bases: list[str] = []
    counter = 0
    while len(bases) < length:
        digest = hashlib.sha256(f"{label}:{counter}".encode("ascii")).digest()
        counter += 1
        for value in digest:
            if value % 100 < gc_percent:
                bases.append("G" if value & 1 else "C")
            else:
                bases.append("A" if value & 1 else "T")
            if len(bases) == length:
                break
    return "".join(bases)


def build_references() -> list[tuple[str, str]]:
    references = []
    for index, name in enumerate(CONTIG_NAMES):
        genome = "genome_A" if index < 4 else "genome_B"
        gc_percent = 42 if genome == "genome_A" else 61
        references.append((name, deterministic_sequence(name, CONTIG_LENGTH, gc_percent)))
    return references


def write_fasta(path: Path, references: list[tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="ascii", newline="\n") as handle:
        for name, sequence in references:
            handle.write(f">{name}\n")
            for offset in range(0, len(sequence), 80):
                handle.write(sequence[offset : offset + 80] + "\n")


def build_alignments(references: list[tuple[str, str]]) -> list[AlignmentSpec]:
    sequences = [sequence for _, sequence in references]
    pairs: list[tuple[str, int, int, int, int]] = []
    for index in range(12):
        pairs.append((f"accepted_{index:02d}", 0, 930, 1, 20))
    for index in range(5):
        pairs.append((f"ambiguous_a_{index:02d}", 2, 930, 3, 20))
        pairs.append((f"ambiguous_b_{index:02d}", 2, 930, 4, 20))
    for index in range(3):
        pairs.append((f"within_{index:02d}", 5, 100, 5, 400))

    alignments: list[AlignmentSpec] = []
    for name, ref_a, pos_a, ref_b, pos_b in pairs:
        alignments.append(
            AlignmentSpec(
                query_name=name,
                reference_id=ref_a,
                position=pos_a,
                mate_reference_id=ref_b,
                mate_position=pos_b,
                read1=True,
                sequence=sequences[ref_a][pos_a : pos_a + 50],
            )
        )
        alignments.append(
            AlignmentSpec(
                query_name=name,
                reference_id=ref_b,
                position=pos_b,
                mate_reference_id=ref_a,
                mate_position=pos_a,
                read1=False,
                sequence=sequences[ref_b][pos_b : pos_b + 50],
            )
        )
    return alignments


def verify_fixture(fasta: Path, bam: Path) -> None:
    if not fasta.is_file() or not bam.is_file():
        raise FileNotFoundError("Toy FASTA/BAM fixture is incomplete")
    names, lengths = read_bam_header(bam)
    if names != CONTIG_NAMES or lengths != [CONTIG_LENGTH] * len(CONTIG_NAMES):
        raise AssertionError("BAM header does not match the toy FASTA contigs")
    records = list(iter_bam_records(bam))
    if len(records) != 50 or len({record.query_name for record in records}) != 25:
        raise AssertionError("Unexpected number of toy BAM alignment records")
    validate_bai(bam.with_suffix(bam.suffix + ".bai"), len(CONTIG_NAMES))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parent / "data",
    )
    parser.add_argument(
        "--verify-only",
        action="store_true",
        help="Validate the bundled files without regenerating them.",
    )
    args = parser.parse_args()
    output_dir = args.output_dir.resolve()
    fasta = output_dir / "toy_contigs.fasta"
    bam = output_dir / "toy_reads.sorted.bam"

    if not args.verify_only:
        references = build_references()
        write_fasta(fasta, references)
        write_bam_and_bai(bam, references, build_alignments(references))
    verify_fixture(fasta, bam)
    print(f"FASTA: {fasta}")
    print(f"BAM:   {bam}")
    print(f"BAI:   {bam.with_suffix(bam.suffix + '.bai')}")
    print("Toy fixture validation passed (8 contigs, 25 read pairs).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
