import gzip
import os
import sys
import glob
from collections import defaultdict
import logging

# 配置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def gen_bins(fastafile: str, resultfile: str, outputdir: str) -> None:
    """
    根据 binning 结果将 contigs 分组并生成每个 bin 的 FASTA 文件
    """
    logger.info("Processing FASTA file: %s", fastafile)
    sequences = {}

    opener, mode = (gzip.open, 'rt') if fastafile.endswith("gz") else (open, 'r')

    with opener(fastafile, mode) as f:
        current_seq_name = None
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current_seq_name = line.split(" ")[0][1:] if " " in line else line[1:]
                sequences[current_seq_name] = []
            elif current_seq_name:
                sequences[current_seq_name].append(line)

    logger.info("Loaded %d contigs.", len(sequences))

    logger.info("Reading result file: %s", resultfile)
    contig_to_bin = {}
    try:
        with open(resultfile, 'r') as f:
            next(f)  # skip header
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    contig_to_bin[parts[0]] = parts[1]
    except FileNotFoundError:
        logger.error("Result file '%s' not found.", resultfile)
        sys.exit(1)

    bins = defaultdict(list)
    for contig_id, bin_id in contig_to_bin.items():
        if contig_id in sequences:
            bins[bin_id].append(contig_id)

    os.makedirs(outputdir, exist_ok=True)

    logger.info("Generating %d bin FASTA files to '%s'.", len(bins), outputdir)
    for bin_id, contig_ids in bins.items():
        out_file = os.path.join(outputdir, f"{bin_id}.fa")
        with open(out_file, 'w') as f:
            for contig_id in contig_ids:
                f.write(f">{contig_id}\n{''.join(sequences[contig_id])}\n")

    logger.info("All FASTA files generated in '%s'.", outputdir)


def main(args,logger):
    CONTIG_FILE = args.contig_file
    RESULT_DIR = os.path.join(args.output_path, 'binning/s_cluster/cluster_res')
    OUTPUT_DIR = os.path.join(args.output_path, 'binning/bins')

    if not os.path.isdir(RESULT_DIR):
        logger.error("'%s' is not a valid directory.", RESULT_DIR)
        sys.exit(1)

    result_files = glob.glob(os.path.join(RESULT_DIR, '*_200000.tsv'))

    if not result_files:
        logger.error("No '_200000.tsv' files found in '%s'.", RESULT_DIR)
        sys.exit(1)
    elif len(result_files) > 1:
        logger.error("Multiple '_200000.tsv' files found in '%s'. Only one expected.", RESULT_DIR)
        for f in result_files:
            logger.error("  - %s", f)
        sys.exit(1)

    RESULT_FILE = result_files[0]
    logger.info("Using result file: %s", RESULT_FILE)
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    gen_bins(fastafile=CONTIG_FILE, resultfile=RESULT_FILE, outputdir=OUTPUT_DIR)
    logger.info("Processing completed!")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Generate bin FASTA files from clustering results")
    parser.add_argument('--contig_file', type=str, required=True, help='Input contigs FASTA file')
    parser.add_argument('--output_path', type=str, required=True, help='Base output directory')
    args = parser.parse_args()
    main(args)
