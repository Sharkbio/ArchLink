#!/usr/bin/env python3
"""Minimal standards-compliant BAM/BAI helpers for the synthetic toy fixture.

The writer intentionally supports only the small subset needed by this fixture:
single-block BGZF, one 50M CIGAR operation, coordinate-sorted paired reads, and
one linear-index interval per short contig.  The reader is sufficient for the
portable smoke test and does not replace pysam/samtools in the full pipeline.
"""

from __future__ import annotations

import gzip
import struct
import zlib
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, Sequence


BGZF_EOF = bytes.fromhex(
    "1f8b08040000000000ff0600424302001b0003000000000000000000"
)
BASE_CODES = {
    "=": 0,
    "A": 1,
    "C": 2,
    "M": 3,
    "G": 4,
    "R": 5,
    "S": 6,
    "V": 7,
    "T": 8,
    "W": 9,
    "Y": 10,
    "H": 11,
    "K": 12,
    "D": 13,
    "B": 14,
    "N": 15,
}


@dataclass(frozen=True)
class AlignmentSpec:
    query_name: str
    reference_id: int
    position: int
    mate_reference_id: int
    mate_position: int
    read1: bool
    sequence: str


@dataclass(frozen=True)
class BamRecord:
    query_name: str
    reference_name: str
    position: int
    mate_reference_name: str
    mate_position: int
    flag: int
    mapping_quality: int


def reg2bin(begin: int, end: int) -> int:
    """Return the BAI bin for a zero-based half-open interval."""

    end -= 1
    if begin >> 14 == end >> 14:
        return 4681 + (begin >> 14)
    if begin >> 17 == end >> 17:
        return 585 + (begin >> 17)
    if begin >> 20 == end >> 20:
        return 73 + (begin >> 20)
    if begin >> 23 == end >> 23:
        return 9 + (begin >> 23)
    if begin >> 26 == end >> 26:
        return 1 + (begin >> 26)
    return 0


def encode_sequence(sequence: str) -> bytes:
    codes = [BASE_CODES.get(base.upper(), 15) for base in sequence]
    if len(codes) % 2:
        codes.append(0)
    return bytes((codes[i] << 4) | codes[i + 1] for i in range(0, len(codes), 2))


def make_bgzf_block(payload: bytes) -> bytes:
    compressor = zlib.compressobj(level=6, wbits=-15)
    compressed = compressor.compress(payload) + compressor.flush()
    block_size = 18 + len(compressed) + 8
    if block_size > 65536:
        raise ValueError("Toy BAM payload exceeds the single-BGZF-block limit")
    header = (
        b"\x1f\x8b\x08\x04"
        + struct.pack("<I", 0)
        + b"\x00\xff"
        + struct.pack("<H", 6)
        + b"BC"
        + struct.pack("<H", 2)
        + struct.pack("<H", block_size - 1)
    )
    footer = struct.pack("<II", zlib.crc32(payload) & 0xFFFFFFFF, len(payload))
    return header + compressed + footer


def build_alignment_record(spec: AlignmentSpec) -> bytes:
    read_name = spec.query_name.encode("ascii") + b"\x00"
    read_length = len(spec.sequence)
    flag = 0x1 | (0x40 if spec.read1 else 0x80)
    bin_mq_nl = (reg2bin(spec.position, spec.position + read_length) << 16) | (
        60 << 8
    ) | len(read_name)
    flag_nc = (flag << 16) | 1
    core = struct.pack(
        "<iiIIiiii",
        spec.reference_id,
        spec.position,
        bin_mq_nl,
        flag_nc,
        read_length,
        spec.mate_reference_id,
        spec.mate_position,
        0,
    )
    cigar = struct.pack("<I", read_length << 4)  # M operation is code 0.
    qualities = bytes([40]) * read_length
    auxiliary = b"RGZtoy\x00"
    body = core + read_name + cigar + encode_sequence(spec.sequence) + qualities + auxiliary
    return struct.pack("<i", len(body)) + body


def write_bam_and_bai(
    bam_path: Path,
    references: Sequence[tuple[str, str]],
    alignments: Iterable[AlignmentSpec],
) -> None:
    """Write a coordinate-sorted BAM and matching BAI index."""

    bam_path.parent.mkdir(parents=True, exist_ok=True)
    header_text = "@HD\tVN:1.6\tSO:coordinate\n"
    header_text += "".join(
        f"@SQ\tSN:{name}\tLN:{len(sequence)}\n" for name, sequence in references
    )
    header_text += "@RG\tID:toy\tSM:ArchLink_toy\n"
    encoded_header = header_text.encode("ascii")
    payload = bytearray(b"BAM\x01")
    payload.extend(struct.pack("<i", len(encoded_header)))
    payload.extend(encoded_header)
    payload.extend(struct.pack("<i", len(references)))
    for name, sequence in references:
        encoded_name = name.encode("ascii") + b"\x00"
        payload.extend(struct.pack("<i", len(encoded_name)))
        payload.extend(encoded_name)
        payload.extend(struct.pack("<i", len(sequence)))

    sorted_alignments = sorted(
        alignments,
        key=lambda item: (item.reference_id, item.position, item.query_name, not item.read1),
    )
    record_offsets: dict[int, list[tuple[int, int, int]]] = {
        index: [] for index in range(len(references))
    }
    for spec in sorted_alignments:
        start = len(payload)
        record = build_alignment_record(spec)
        payload.extend(record)
        end = len(payload)
        record_offsets[spec.reference_id].append((start, end, spec.position))

    if len(payload) >= 60000:
        raise ValueError("Toy BAM must remain small enough for one BGZF block")
    bam_path.write_bytes(make_bgzf_block(bytes(payload)) + BGZF_EOF)

    bai_path = bam_path.with_suffix(bam_path.suffix + ".bai")
    with bai_path.open("wb") as handle:
        handle.write(b"BAI\x01")
        handle.write(struct.pack("<i", len(references)))
        for ref_id in range(len(references)):
            offsets = record_offsets[ref_id]
            if not offsets:
                handle.write(struct.pack("<i", 0))
                # An explicit zero linear-index entry is accepted by htslib and
                # avoids edge-case failures in lightweight pure-Python readers.
                handle.write(struct.pack("<i", 1))
                handle.write(struct.pack("<Q", 0))
                continue
            chunk_begin = offsets[0][0]
            chunk_end = offsets[-1][1]
            bin_id = reg2bin(offsets[0][2], offsets[-1][2] + 50)
            handle.write(struct.pack("<i", 1))
            handle.write(struct.pack("<I", bin_id))
            handle.write(struct.pack("<i", 1))
            handle.write(struct.pack("<QQ", chunk_begin, chunk_end))
            handle.write(struct.pack("<i", 1))
            handle.write(struct.pack("<Q", chunk_begin))
        handle.write(struct.pack("<Q", 0))


def read_bam_header(path: Path) -> tuple[list[str], list[int]]:
    with gzip.open(path, "rb") as handle:
        if handle.read(4) != b"BAM\x01":
            raise ValueError(f"Not a BAM file: {path}")
        (header_length,) = struct.unpack("<i", handle.read(4))
        handle.read(header_length)
        (reference_count,) = struct.unpack("<i", handle.read(4))
        names: list[str] = []
        lengths: list[int] = []
        for _ in range(reference_count):
            (name_length,) = struct.unpack("<i", handle.read(4))
            names.append(handle.read(name_length)[:-1].decode("ascii"))
            (length,) = struct.unpack("<i", handle.read(4))
            lengths.append(length)
    return names, lengths


def iter_bam_records(path: Path) -> Iterator[BamRecord]:
    with gzip.open(path, "rb") as handle:
        if handle.read(4) != b"BAM\x01":
            raise ValueError(f"Not a BAM file: {path}")
        (header_length,) = struct.unpack("<i", handle.read(4))
        handle.read(header_length)
        (reference_count,) = struct.unpack("<i", handle.read(4))
        references: list[str] = []
        for _ in range(reference_count):
            (name_length,) = struct.unpack("<i", handle.read(4))
            references.append(handle.read(name_length)[:-1].decode("ascii"))
            handle.read(4)

        while True:
            size_bytes = handle.read(4)
            if not size_bytes:
                break
            if len(size_bytes) != 4:
                raise ValueError("Truncated BAM block size")
            (block_size,) = struct.unpack("<i", size_bytes)
            core = handle.read(32)
            if len(core) != 32:
                raise ValueError("Truncated BAM alignment core")
            (
                reference_id,
                position,
                bin_mq_nl,
                flag_nc,
                _read_length,
                mate_reference_id,
                mate_position,
                _template_length,
            ) = struct.unpack("<iiIIiiii", core)
            variable = handle.read(block_size - 32)
            read_name_length = bin_mq_nl & 0xFF
            query_name = variable[: read_name_length - 1].decode("ascii")
            flag = flag_nc >> 16
            mapping_quality = (bin_mq_nl >> 8) & 0xFF
            yield BamRecord(
                query_name=query_name,
                reference_name=references[reference_id],
                position=position,
                mate_reference_name=references[mate_reference_id],
                mate_position=mate_position,
                flag=flag,
                mapping_quality=mapping_quality,
            )


def validate_bai(path: Path, expected_reference_count: int) -> None:
    with path.open("rb") as handle:
        if handle.read(4) != b"BAI\x01":
            raise ValueError(f"Not a BAI index: {path}")
        (reference_count,) = struct.unpack("<i", handle.read(4))
    if reference_count != expected_reference_count:
        raise ValueError(
            f"BAI reference count {reference_count} != expected {expected_reference_count}"
        )
