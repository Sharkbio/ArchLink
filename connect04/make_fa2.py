#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Make FASTA from matching result.

Usage:
    python make_fa.py <input_fasta> <matching_result> <output_fasta>

Arguments:
    input_fasta      : Original assembly in FASTA format.
    matching_result  : Output file from the 'matching' program (tab-separated, with orientation suffixes like '+', '-').
    output_fasta     : Output FASTA file with contigs joined according to the matching result.
"""

import sys
import re
from typing import Dict, List, Tuple


def print_help() -> None:
    """Print usage help."""
    print(__doc__.strip())


def complement(seq: str) -> str:
    """
    Return the complement of a DNA sequence.
    
    Handles both uppercase and lowercase bases. Non-ATCG characters are left unchanged.
    """
    comp_table = str.maketrans('ATCGatcg', 'TAGCtagc')
    return seq.translate(comp_table)


def parse_fasta(file_path: str) -> Dict[str, str]:
    """
    Parse a FASTA file into a dictionary mapping sequence IDs to sequences.
    
    Only the first word of the header (after '>') is used as the ID.
    """
    sequences = {}
    current_id = None
    current_seq_lines = []

    try:
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    if current_id is not None:
                        sequences[current_id] = ''.join(current_seq_lines)
                    # Extract ID: take first token after '>'
                    current_id = line[1:].split()[0]
                    current_seq_lines = []
                else:
                    current_seq_lines.append(line)
            # Handle last sequence
            if current_id is not None:
                sequences[current_id] = ''.join(current_seq_lines)
    except FileNotFoundError:
        raise FileNotFoundError(f"FASTA file not found: {file_path}")
    
    return sequences


def parse_matching_result(file_path: str) -> List[List[str]]:
    """
    Parse the matching result file.
    
    Each line is split by tab. Lines starting with 'iter' or 'self' are skipped.
    Duplicate entries (based on the first element without orientation) are removed.
    """
    seen_ids = set()
    filtered_lines = []

    try:
        with open(file_path, 'r') as f:
            for line_num, line in enumerate(f, 1):
                stripped = line.strip()
                if not stripped:
                    continue
                if stripped.startswith(('iter', 'self')):
                    continue

                parts = stripped.split('\t')
                if not parts:
                    continue

                # Use the ID part (without last char, which is '+' or '-')
                first_id = parts[0][:-1] if len(parts[0]) > 1 else parts[0]

                if first_id in seen_ids:
                    print(f"Warning: Skipping duplicate entry for '{first_id}' at line {line_num}", file=sys.stderr)
                    continue

                seen_ids.add(first_id)
                filtered_lines.append(parts)
    except FileNotFoundError:
        raise FileNotFoundError(f"Matching result file not found: {file_path}")

    return filtered_lines


def main(input_fasta,matching_result,output_fasta):

    # Load FASTA sequences
    try:
        record_dict = parse_fasta(input_fasta)
    except Exception as e:
        print(f"Error reading FASTA file: {e}", file=sys.stderr)
        sys.exit(1)

    # Parse matching result
    try:
        ordered_groups = parse_matching_result(matching_result)
    except Exception as e:
        print(f"Error reading matching result: {e}", file=sys.stderr)
        sys.exit(1)

    used_ids = set()

    try:
        with open(output_fasta, 'w') as out_f:
            # Process each group from matching result
            for group in ordered_groups:
                full_seq = ""
                group_ids = []

                for item in group:
                    if not item:
                        continue
                    # Last character is orientation ('+' or '-')
                    if item[-1] in '+-':
                        seq_id = item[:-1]
                        orientation = item[-1]
                    else:
                        seq_id = item
                        orientation = '+'  # default

                    if seq_id not in record_dict:
                        print(f"Warning: Sequence '{seq_id}' not found in FASTA. Skipping.", file=sys.stderr)
                        continue

                    seq = record_dict[seq_id]
                    if orientation == '-':
                        seq = complement(seq)[::-1]  # reverse complement

                    full_seq += seq
                    group_ids.append(item)
                    used_ids.add(seq_id)

                if full_seq:
                    header = ">{}".format("_".join(group_ids))
                    out_f.write(f"{header}\n")
                    out_f.write(f"{full_seq}\n")

            # Write unused sequences (not in any group)
            for seq_id, seq in record_dict.items():
                if seq_id not in used_ids:
                    out_f.write(f">{seq_id}\n")
                    out_f.write(f"{seq}\n")

    except IOError as e:
        print(f"Error writing output file: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"Successfully wrote output to: {output_fasta}")


if __name__ == "__main__":
    main()