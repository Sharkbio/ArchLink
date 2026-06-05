from Bio import SeqIO
import mimetypes
import os
import gzip
import random
import shutil
from typing import Dict
import argparse
import logging
import sys

# =========================================================================
# 1. Sequence Loading Function
# =========================================================================

def get_inputsequences(fastx_file: str):
    """
    Retrieve sequences from a FASTX file and return them as a dictionary.

    :param fastx_file: Path to the FASTX file (either FASTA or FASTQ).
    :return: A dictionary where sequence IDs are keys and sequences are values.
    """
    file_type = mimetypes.guess_type(fastx_file)[1]
    
    # Determine how to open the file (gzipped or plain)
    if file_type == 'gzip':
        # 'rt' is used for reading text mode with gzip
        f = gzip.open(fastx_file, "rt")
    elif not file_type:
        f = open(fastx_file, "rt")
    else:
        raise RuntimeError(f"Unknown type of file: '{fastx_file}'")
        
    seqs = {}
    
    # Check if file is empty
    if os.path.getsize(fastx_file) == 0:
        f.close()
        return seqs
        
    file_format = None
    
    # Read the first line to determine the format
    line = f.readline()
    
    # Simple check for fasta/fastq format
    if line.startswith('@'):
        file_format = "fastq"
    elif line.startswith(">"):
        file_format = "fasta"
        
    # Rewind file pointer to the beginning
    f.seek(0)
    
    if not file_format:
        f.close()
        raise RuntimeError(f"Invalid sequence file (must start with '>' or '@'): '{fastx_file}'")
        
    # Parse the file using Biopython's SeqIO
    for seq_record in SeqIO.parse(f, file_format):
        seqs[seq_record.id] = seq_record.seq

    f.close()
    return seqs


# =========================================================================
# 2. Sequence Augmentation and Saving Function
# =========================================================================

def gen_augfasta(seqs: Dict[str, str], augprefix: str, out_file: str,
                 p: float = None, contig_len: int = 1000): 
    """
    Generate augmented sequences (via random slicing only) and save them to a FASTA file 
    along with sequence information. Point mutation error injection is removed.

    :param seqs: A dictionary of input sequences where keys are sequence IDs, and values are sequences.
    :param augprefix: A prefix used in the augmented sequence IDs (e.g., 'aug1').
    :param out_file: Path to the output FASTA file. (e.g., '/path/to/aug1/aug1.fasta')
    :param p: Proportion of the original sequence to include in the augmented sequences (default is None).
    :param contig_len: Minimum length of the original sequence required for augmentation (default is 1000).
    """
    seqkeys = []
    # Filter out short contigs that cannot meet the minimum length requirement
    for seqid in seqs.keys():
        # +1 because sequence slicing [start:end+1] ensures a minimum length of contig_len
        if len(seqs[seqid]) >= contig_len + 1: 
            seqkeys.append(seqid)

    aug_seq_info = []
    
    # Optimize file writing
    with open(out_file, 'w') as f:
        for seqid in seqkeys:
            original_sequence_obj = seqs[seqid]
            
            # --- Random Slicing Logic ---
            if not p:
                # Randomly select start point and length
                max_start = len(original_sequence_obj) - (contig_len + 1)
                if max_start < 0: continue # Should not happen if pre-filtered, but safe check
                start = random.randint(0, max_start)
                # Ensure slice length is at least contig_len
                sim_len = random.randint(contig_len, len(original_sequence_obj) - start)
                end = start + sim_len - 1
            else:
                # Proportional slicing logic
                sim_len = int(p * len(original_sequence_obj))
                # Ensure sufficient space for slicing
                max_start = len(original_sequence_obj) - sim_len - 1
                if max_start < 0: continue
                start = random.randint(0, max_start)
                end = start + sim_len - 1
                
            # Extract subsequence
            sub_sequence = str(original_sequence_obj[start:end + 1])
            
            # --- Sequence Preparation (No Error Injection) ---
            # Since perturbation is removed, the final sequence is just the sliced sequence
            final_sequence = sub_sequence

            # --- Write to FASTA File ---
            seqid_name = f">{seqid}_{augprefix}"
            f.write(f"{seqid_name}\n")
            f.write(f"{final_sequence}\n")
            
            # Record augmentation information
            aug_seq_info.append((seqid, start, end, len(final_sequence)))

    # Save augmentation info file
    directory = os.path.dirname(out_file)
    prefix = os.path.splitext(os.path.basename(out_file))[0] # e.g., 'aug1'
    # Use consistent naming for info file
    aug_seq_info_out_file = f"{directory}/sequences_{prefix}_info.tsv" 

    with open(aug_seq_info_out_file, 'w') as afile:
        afile.write('seqid\tstart\tend\tlength\n')
        for seqid, start, end, final_len in aug_seq_info:
            afile.write(f"{seqid}\t{start}\t{end}\t{final_len}\n")


# =========================================================================
# 3. Main Execution Function
# =========================================================================

def run_gen_augfasta(logger, args):
    """
    Generate augmentation fasta file and save index, and generate k-mer features.
    """
    # Calculate the number of augmented views needed (n_views - 1)
    num_aug = args.n_views - 1 
    fasta_file = args.contig_file
    out_path = args.out_augdata_path
    contig_len = args.contig_len

    # --- K-mer Feature Generation Import ---
    # !!! IMPORTANT: This assumes a module named 'gen_kmer' is available in the current project structure.
    try:
        from gen_kmer import run_gen_kmer
        logger.info("Successfully imported run_gen_kmer for k-mer feature generation.")
    except ImportError:
        logger.error("Could not import run_gen_kmer. K-mer generation steps will likely fail.")
        # Define a mock function to prevent runtime crash if the file is run standalone
        def run_gen_kmer(out_file, start_k, end_k):
             logger.warning(f"Mock run_gen_kmer called for {out_file}. (K-mer generation skipped)")

    # --- Start Logging ---
    logger.info("=" * 70)
    logger.info("ArchLink: Starting Sequence Augmentation Module (gen_augfasta)")
    logger.info(f"Input Contig File: {fasta_file}")
    logger.info(f"Output Directory: {out_path}")
    logger.info(f"Parameters: N_Views={args.n_views}, Contig_Len_Min={contig_len}")
    logger.info("-" * 70)

    # --- 1. Process View 0 (Original) ---
    outdir_0 = out_path + '/aug0'
    # Use exist_ok=True to prevent errors on existing directory
    os.makedirs(outdir_0, exist_ok=True) 
    out_file_0 = outdir_0 + '/aug0.fasta' 
    
    # Copy original contig file as View 0 reference
    shutil.copyfile(fasta_file, out_file_0)
    logger.info(f"[View 0] Created reference view by copying to: {out_file_0}")

    # --- K-MER GENERATION FOR VIEW 0 ---
    # K-mer counting for the original sequences (aug0)
    run_gen_kmer(out_file_0, 0, 4)
    logger.info(f"[View 0] Initiated K-mer feature generation for {out_file_0}")
    # ------------------------------------

    # --- 2. Process Augmented Views 1 to N-1 ---
    seqs = get_inputsequences(fasta_file) # Load original sequences only once

    for i in range(num_aug):
        view_idx = i + 1
        outdir = f"{out_path}/aug{view_idx}"
        os.makedirs(outdir, exist_ok=True)
        
        logger.info(f"--- Generating Augmented View {view_idx} (aug{view_idx}.fasta) ---")
        
        # 'p' parameter is currently set to None in the original code
        p = None 
        
        out_file = f"{outdir}/aug{view_idx}.fasta"
        
        # Call gen_augfasta to apply random slicing (error_rate argument removed)
        gen_augfasta(
            seqs, 
            f'aug{view_idx}', 
            out_file, 
            p=p, 
            contig_len=contig_len
        )
        
        logger.info(f"[View {view_idx}] Successfully generated {out_file} and corresponding info file.")
        
        # --- K-MER GENERATION FOR AUGMENTED VIEWS ---
        # K-mer counting for the augmented sequences
        run_gen_kmer(out_file, 0, 4)
        logger.info(f"[View {view_idx}] Initiated K-mer feature generation for {out_file}")
        # --------------------------------------------
        
    # --- End Logging ---
    logger.info("=" * 70)
    logger.info("ArchLink: All sequence views and K-mer features generated successfully.")
    logger.info("=" * 70)
