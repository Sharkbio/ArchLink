#!/usr/bin/env python3
"""Run the portable ArchLink reviewer smoke test on the bundled toy fixture.

This test extracts paired-end support from the BAM, applies a documented
reciprocal-best/abstention rule to the tiny candidate set, and delegates FASTA
joining to ArchLink's released ``connect04.make_fa2`` implementation.  It does
not retrain models or invoke CheckM; those remain part of the optional full
end-to-end workflow in ``config.toy.yaml``.
"""

from __future__ import annotations

import argparse
import csv
import shutil
import sys
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from connect04 import make_fa2  # noqa: E402
from example.toy_bam import BamRecord, iter_bam_records  # noqa: E402


END_WINDOW_BP = 100
MIN_SPANNING_PAIRS = 8
LINK_FIELDS = [
    "contig_a",
    "end_a",
    "contig_b",
    "end_b",
    "spanning_pairs",
    "status",
    "reason",
]


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    name: str | None = None
    chunks: list[str] = []
    with path.open("r", encoding="ascii") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(chunks)))
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if name is not None:
        records.append((name, "".join(chunks)))
    return records


def contig_end(record: BamRecord, lengths: dict[str, int]) -> str:
    if record.position < END_WINDOW_BP:
        return "L"
    if record.position + 50 > lengths[record.reference_name] - END_WINDOW_BP:
        return "R"
    return "internal"


def count_cross_contig_pairs(
    bam_path: Path, lengths: dict[str, int]
) -> dict[tuple[tuple[str, str], tuple[str, str]], set[str]]:
    evidence: dict[tuple[tuple[str, str], tuple[str, str]], set[str]] = defaultdict(set)
    seen: dict[str, BamRecord] = {}
    for record in iter_bam_records(bam_path):
        mate = seen.pop(record.query_name, None)
        if mate is None:
            seen[record.query_name] = record
            continue
        if mate.reference_name == record.reference_name:
            continue
        first = (mate.reference_name, contig_end(mate, lengths))
        second = (record.reference_name, contig_end(record, lengths))
        if "internal" in (first[1], second[1]):
            continue
        edge = tuple(sorted((first, second)))
        evidence[edge].add(record.query_name)
    if seen:
        raise AssertionError(f"Unpaired BAM records remain: {sorted(seen)[:5]}")
    return evidence


def classify_candidates(
    evidence: dict[tuple[tuple[str, str], tuple[str, str]], set[str]]
) -> list[dict[str, object]]:
    supports = {edge: len(names) for edge, names in evidence.items()}
    incident: dict[tuple[str, str], list[tuple[tuple[str, str], tuple[str, str]]]] = defaultdict(list)
    for edge in supports:
        incident[edge[0]].append(edge)
        incident[edge[1]].append(edge)

    rows: list[dict[str, object]] = []
    for edge in sorted(supports):
        support = supports[edge]
        failures: list[str] = []
        for endpoint in edge:
            endpoint_edges = incident[endpoint]
            best = max(supports[item] for item in endpoint_edges)
            best_count = sum(supports[item] == best for item in endpoint_edges)
            if support < best:
                failures.append(f"not best at {endpoint[0]}:{endpoint[1]}")
            elif best_count > 1:
                failures.append(f"tied at {endpoint[0]}:{endpoint[1]}")
        if support < MIN_SPANNING_PAIRS:
            failures.append(f"support < {MIN_SPANNING_PAIRS}")
        accepted = not failures
        rows.append(
            {
                "contig_a": edge[0][0],
                "end_a": edge[0][1],
                "contig_b": edge[1][0],
                "end_b": edge[1][1],
                "spanning_pairs": support,
                "status": "accepted" if accepted else "abstained",
                "reason": "unique reciprocal best" if accepted else "; ".join(failures),
            }
        )
    return rows


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_unconnected_ends(
    path: Path,
    contigs: list[str],
    candidates: list[dict[str, object]],
) -> None:
    connected: set[tuple[str, str]] = set()
    incident: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in candidates:
        endpoints = [
            (str(row["contig_a"]), str(row["end_a"])),
            (str(row["contig_b"]), str(row["end_b"])),
        ]
        for endpoint in endpoints:
            incident[endpoint].append(row)
        if row["status"] == "accepted":
            connected.update(endpoints)

    rows = []
    for contig in contigs:
        for end in ("L", "R"):
            endpoint = (contig, end)
            if endpoint in connected:
                continue
            endpoint_rows = incident.get(endpoint, [])
            reasons = " | ".join(sorted({str(row["reason"]) for row in endpoint_rows}))
            rows.append(
                {
                    "contig_id": contig,
                    "end": end,
                    "reason": reasons or "no candidate evidence",
                }
            )
    write_tsv(path, ["contig_id", "end", "reason"], rows)


def write_matching_result(path: Path, accepted: list[dict[str, object]]) -> None:
    with path.open("w", encoding="ascii", newline="\n") as handle:
        for row in accepted:
            if row["end_a"] != "R" or row["end_b"] != "L":
                raise ValueError("Toy fixture currently supports only R-to-L accepted joins")
            handle.write(f"{row['contig_a']}+\t{row['contig_b']}+\n")


def write_expected(output: Path, expected: Path) -> None:
    expected.mkdir(parents=True, exist_ok=True)
    mapping = {
        "bins.tsv": "expected_bins.tsv",
        "candidate_links.tsv": "expected_candidates.tsv",
        "final_joins.tsv": "expected_links.tsv",
        "unconnected_ends.tsv": "expected_unconnected_ends.tsv",
        "scaffolds.fasta": "expected_scaffolds.fasta",
    }
    for source_name, target_name in mapping.items():
        shutil.copyfile(output / source_name, expected / target_name)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", type=Path, default=ROOT / "example" / "data" / "toy_contigs.fasta")
    parser.add_argument("--bam", type=Path, default=ROOT / "example" / "data" / "toy_reads.sorted.bam")
    parser.add_argument("--bins", type=Path, default=ROOT / "example" / "precomputed" / "contig_bins.tsv")
    parser.add_argument("--output", type=Path, default=ROOT / "example" / "output" / "toy_run")
    parser.add_argument("--write-expected", type=Path, help="Refresh expected files from the observed output.")
    args = parser.parse_args()

    if not args.fasta.is_file() or not args.bam.is_file() or not args.bins.is_file():
        raise FileNotFoundError("Toy FASTA, BAM and precomputed bin table are required")
    contigs = read_fasta(args.fasta)
    lengths = {name: len(sequence) for name, sequence in contigs}
    candidates = classify_candidates(count_cross_contig_pairs(args.bam, lengths))
    accepted = [row for row in candidates if row["status"] == "accepted"]

    args.output.mkdir(parents=True, exist_ok=True)
    for name in (
        "bins.tsv",
        "candidate_links.tsv",
        "final_joins.tsv",
        "unconnected_ends.tsv",
        "matching_result.tsv",
        "scaffolds.fasta",
    ):
        target = args.output / name
        if target.exists():
            target.unlink()
    shutil.copyfile(args.bins, args.output / "bins.tsv")
    write_tsv(args.output / "candidate_links.tsv", LINK_FIELDS, candidates)
    write_tsv(args.output / "final_joins.tsv", LINK_FIELDS, accepted)
    write_unconnected_ends(args.output / "unconnected_ends.tsv", [name for name, _ in contigs], candidates)
    matching_result = args.output / "matching_result.tsv"
    write_matching_result(matching_result, accepted)
    make_fa2.main(str(args.fasta), str(matching_result), str(args.output / "scaffolds.fasta"))

    if args.write_expected:
        write_expected(args.output, args.write_expected.resolve())
    print(f"Candidate links: {len(candidates)}")
    print(f"Accepted joins: {len(accepted)}")
    print(f"Abstained candidates: {sum(row['status'] == 'abstained' for row in candidates)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
