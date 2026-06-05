import leidenalg
import numpy as np
import pandas as pd
import time
import os
import logging
import argparse
from igraph import Graph
import multiprocessing
from typing import List, Optional, Union

# -------- Logging --------
logger = logging.getLogger('ArchLink_Clustering')
logger.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(message)s')
console_hdr = logging.StreamHandler()
console_hdr.setFormatter(formatter)
logger.addHandler(console_hdr)

def save_result(labels: np.ndarray, output_file: str, namelist: List[str]):
    """Save clustering results to TSV."""
    try:
        with open(output_file, 'w') as f:
            for i, label in enumerate(labels):
                f.write(f"group{label}\t{namelist[i]}\n")
        logger.info(f"Results saved to {output_file}")
    except Exception as e:
        logger.error(f"Failed to save results to {output_file}: {e}")

def gen_bins_from_tsv(contig_file: str, tsv_file: str, output_prefix: str):
    """Generate per-bin files from clustering results."""
    try:
        df = pd.read_csv(tsv_file, sep='\t', header=None, names=['bin_id', 'contig_name'])
        for bin_id in df['bin_id'].unique():
            out = f"{output_prefix}_{bin_id}.txt"
            os.makedirs(os.path.dirname(output_prefix), exist_ok=True)
            df[df['bin_id'] == bin_id]['contig_name'].to_csv(out, index=False, header=False)
        logger.info(f"Generated bin files at {output_prefix}_*.txt")
    except Exception as e:
        logger.error(f"Failed to generate bin files: {e}")

def gen_seed_idx(seedURL: str, contig_id_list: List[str]) -> List[int]:
    """Return index list of seed contigs."""
    try:
        seed_list = [line.strip() for line in open(seedURL)]
    except FileNotFoundError:
        logger.error(f"Seed file not found: {seedURL}. No fixed membership applied.")
        return []
    name_map = {name: i for i, name in enumerate(contig_id_list)}
    return [name_map[s] for s in seed_list if s in name_map]

# -------- Leiden clustering --------
def run_leiden_clustering(
    output_file: str, namelist: List[str],
    sources: np.ndarray, targets: np.ndarray, weights: np.ndarray,
    length_weight: List[float], bandwidth: float, resolution_parameter: int,
    initial_list: Optional[List[Union[int, None]]],
    is_membership_fixed: Optional[bool]
):
    """Run Leiden clustering and write results."""
    logger.info(f"Running Leiden: bandwidth={bandwidth}, resolution={resolution_parameter}")

    g = Graph(len(namelist), list(zip(sources, targets)))
    if len(weights) != g.ecount():
        logger.error(f"Weight length ({len(weights)}) does not match edges ({g.ecount()}). Skipped.")
        return

    try:
        partition = leidenalg.RBERVertexPartition(
            g,
            weights=weights,
            initial_membership=initial_list,
            resolution_parameter=resolution_parameter,
            node_sizes=length_weight
        )
        optimiser = leidenalg.Optimiser()
        optimiser.optimise_partition(partition, is_membership_fixed=is_membership_fixed, n_iterations=-1)

        # Write clustering output
        with open(output_file, 'w') as f:
            for cluster_id, nodes in enumerate(partition):
                for node in nodes:
                    f.write(f"group{cluster_id}\t{namelist[node]}\n")
        logger.info(f"Written: {output_file}")

    except Exception as e:
        logger.error(f"Leiden failed: {e}")
        raise

# -------- Parallel batch runner --------
def run_all_clusterings(logger, input_dir, seed_file, contig_file,
                        output_dir, num_threads, partgraph_ratio,
                        bandwidth, max_edges_value):
    """Load data and run Leiden clustering in parallel."""
    logger.info(f"Starting batch: max_edges={max_edges_value}")

    output_cluster_dir = os.path.join(output_dir, 'cluster_res')
    os.makedirs(output_cluster_dir, exist_ok=True)

    # Load input data
    try:
        data = np.load(os.path.join(input_dir, 'extracted_edges.npz'))
        sources, targets, weights = data['sources'], data['targets'], data['weights']
        if not np.issubdtype(targets.dtype, np.integer):
            targets = targets.astype(np.int64)
        logger.info(f"Loaded {len(sources)} edges")
    except Exception as e:
        logger.error(f"Failed loading edges: {e}")
        return

    # Load metadata
    try:
        namelist = pd.read_csv(os.path.join(input_dir, 'namelist.txt'), header=None)[0].tolist()
        length_weight = pd.read_csv(os.path.join(input_dir, 'length_weight.txt'), header=None)[0].tolist()
    except Exception as e:
        logger.error(f"Failed loading metadata: {e}")
        return

    # Seed/fixed membership
    seed_idx = gen_seed_idx(seed_file, namelist)
    initial_list = list(range(len(namelist)))
    is_membership_fixed = [i in seed_idx for i in initial_list]

    # Parameter list
    resolution_parameter_list = [1, 5, 10, 30, 60, 90, 110,130,150,200]

    # Run in parallel
    async_results = []
    with multiprocessing.Pool(num_threads) as pool:
        for res_param in resolution_parameter_list:
            outname = (
                f"Leiden_bandwidth_{bandwidth}_res_maxedges{max_edges_value}"
                f"_respara_{res_param}_partgraph_ratio_{partgraph_ratio}.tsv"
            )
            outfile = os.path.join(output_cluster_dir, outname)

            if os.path.exists(outfile):
                logger.info(f"Exists, skipped: {outfile}")
                continue

            logger.info(f"Submitting: {outfile}")
            async_results.append(pool.apply_async(
                run_leiden_clustering,
                (outfile, namelist, sources, targets, weights,
                 length_weight, bandwidth, res_param,
                 initial_list, is_membership_fixed)
            ))

        for r in async_results:
            try:
                r.get()
            except Exception as e:
                logger.error(f"Task failed: {e}")

    logger.info(f"Completed batch: max_edges={max_edges_value}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run Leiden clustering in parallel.')
    parser.add_argument('--input_dir')
    parser.add_argument('--seed_file')
    parser.add_argument('--contig_file')
    parser.add_argument('--output_path')
    parser.add_argument('--num_threads', type=int, default=1)
    parser.add_argument('--bandwidth', type=float)
    parser.add_argument('--max_edges_values', type=str, default='100')
    parser.add_argument('--partgraph_ratio', type=int)

    args = parser.parse_args()

    # Parse max_edges list
    try:
        MAX_EDGES_LIST = [int(x) for x in args.max_edges_values.split(',') if x.strip()]
        if not MAX_EDGES_LIST:
            raise ValueError
    except:
        logger.error("Invalid --max_edges_values. Falling back to [60, 80, 100].")
        MAX_EDGES_LIST = [60, 80, 100]

    # Run all batches
    for max_e in MAX_EDGES_LIST:
        run_all_clusterings(
            logger, args.input_dir, args.seed_file, args.contig_file,
            args.output_path, args.num_threads,
            args.partgraph_ratio, args.bandwidth, max_e
        )

    logger.info("All clustering finished.")
