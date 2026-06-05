import hnswlib
import numpy as np
import pandas as pd
import time
import os
import logging
import argparse
import re # Import regular expression module
from sklearn.preprocessing import normalize
from typing import List, Optional, Union

# --- Logger Setup ---
logger = logging.getLogger('ArchLink_EdgeExtraction')
logger.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(message)s')
console_hdr = logging.StreamHandler()
console_hdr.setFormatter(formatter)
logger.addHandler(console_hdr)

# --- Utility Functions ---

def get_length(contig_file: str) -> dict:
    """
    Parses a FASTA file to extract contig names and their corresponding lengths
    by accumulating lengths from sequence lines. The contig ID is extracted as
    the first word after the '>' in the header line.

    :param contig_file: Path to the FASTA contig file.
    :return: A dictionary where keys are contig IDs (e.g., '1_k141_2012526')
             and values are their lengths calculated from sequence data.
    """
    logger.info(f"Extracting contig lengths from {contig_file} by accumulating sequence data.")
    lengths = {}
    try:
        with open(contig_file, 'r') as f:
            current_contig_id = None
            current_contig_length = 0
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    # If this is a new header, store the length of the previous contig
                    if current_contig_id is not None and current_contig_length > 0:
                        lengths[current_contig_id] = current_contig_length
                    
                    # Reset for the new contig. Extract ID as the first word after '>'
                    current_contig_id = line.lstrip('>').split(' ')[0]
                    current_contig_length = 0 # Reset length for the new contig
                else:
                    # Accumulate length from sequence lines
                    if current_contig_id is not None:
                        current_contig_length += len(line)
            
            # Store the length of the last contig after the loop finishes
            if current_contig_id is not None and current_contig_length > 0:
                lengths[current_contig_id] = current_contig_length

    except FileNotFoundError:
        logger.error(f"Contig file not found: {contig_file}. Please ensure the path is correct.")
        return {}
    except Exception as e:
        logger.error(f"Error reading contig file {contig_file}: {e}")
        return {}
    
    logger.info(f"Finished extracting lengths for {len(lengths)} contigs.")
    return lengths

def calculateN50(lengths: List[float]) -> float:
    """
    Calculates the N50 value for a list of lengths.
    N50 is a common metric used in genome assembly to describe contiguity.
    
    :param lengths: A list of contig or scaffold lengths.
    :return: The N50 value as a float. Returns 0.0 if the input list is empty.
    """
    if not lengths:
        return 0.0
    
    # Sort lengths in descending order
    sorted_lengths = sorted(lengths, reverse=True)
    
    # Calculate total length
    total_length = sum(sorted_lengths)
    
    cumulative_length = 0
    for length in sorted_lengths:
        cumulative_length += length
        # N50 is the smallest length L such that 50% of the total length
        # is contained in sequences of length L or greater.
        if cumulative_length >= total_length / 2:
            return float(length)
            
    return 0.0 # Should ideally not be reached if lengths is not empty

# --- Core HNSW and Edge Extraction Functions ---

def fit_hnsw_index(logger, features: np.ndarray, ef: int = 100, M: int = 16,
                   space: str = 'l2', save_index_file: bool = False) -> hnswlib.Index:
    """
    Fit an HNSW index with the given features using the HNSWlib library;
    Convenience function to create HNSW graph for efficient similarity search.

    :param logger: The logger object for logging messages.
    :param features: A numpy array containing the embeddings.
    :param ef: The ef parameter to tune the HNSW algorithm (default: 100).
    :param M: The M parameter to tune the HNSW algorithm (default: 16).
    :param space: The space in which the index operates ('l2', 'cosine', or 'ip', default: 'l2').
    :param save_index_file: Path to save the HNSW index file (optional).

    :return: The HNSW index created using the given features.
    """
    time_start = time.time()
    num_elements = len(features)
    labels_index = np.arange(num_elements)
    EMBEDDING_SIZE = features.shape[1]

    # Declaring index
    p = hnswlib.Index(space=space, dim=EMBEDDING_SIZE)

    # Initing index - the maximum number of elements should be known
    p.init_index(max_elements=num_elements, ef_construction=ef, M=M)

    # Element insertion
    p.add_items(features, labels_index)

    # Controlling the recall by setting ef; ef should always be > k (max_edges)
    p.set_ef(ef)

    # If you want to save the graph to a file
    if save_index_file:
        p.save_index(save_index_file)
    time_end = time.time()
    logger.info(f'Time cost for HNSW index fitting: {time_end - time_start:.2f}s')
    return p

def extract_edges(logger, args):
    """
    Extracts edges and their weights from embeddings using HNSW and saves them to files.
    This is the first step of the two-part process.
    """
    logger.info("Starting edge extraction.")

    emb_file = args.embeding_path
    contig_file = args.contig_file
    output_dir = args.output_path+'/binning/s_cluster/leiden_initial_edge'
    contig_len_threshold = args.contig_len
    not_l2normalize = args.not_l2normalize
    max_edges = args.max_edges
    partgraph_ratio = args.partgraph_ratio
    lmode = args.lmode
    bandwidth_for_edge_extraction = args.bandwidth_for_edge_extraction

    os.makedirs(output_dir, exist_ok=True)

    # Load embeddings and contig names
    try:
        embHeader = pd.read_csv(emb_file, sep='\t', nrows=1)
        embMat_full = pd.read_csv(emb_file, sep='\t', usecols=range(1, embHeader.shape[1])).values
        namelist_full = pd.read_csv(emb_file, sep='\t', usecols=range(1)).values[:, 0]
    except FileNotFoundError:
        logger.error(f"Embeddings file not found: {emb_file}. Exiting.")
        return
    except Exception as e:
        logger.error(f"Error loading embeddings file {emb_file}: {e}. Exiting.")
        return

    # Get contig lengths
    lengths_map = get_length(contig_file)
    # Ensure all contigs in namelist_full have a length, default to 0 if not found
    length_weight_full = np.array([lengths_map.get(seq_id, 0) for seq_id in namelist_full])

    # Filter based on contig length threshold
    mask = length_weight_full >= contig_len_threshold
    embMat_filtered = embMat_full[mask]
    namelist_filtered = namelist_full[mask]
    length_weight_filtered = length_weight_full[mask]

    logger.info(f"Number of contigs after length filtering: {len(namelist_filtered)}")
    if len(namelist_filtered) == 0:
        logger.error("No contigs left after length filtering. Cannot proceed with edge extraction. Exiting.")
        return

    N50 = calculateN50(list(length_weight_filtered))
    logger.info(f'N50 after filtering: {N50:.2f}')

    # Normalize embeddings if not skipped
    if not_l2normalize:
        norm_embeddings = embMat_filtered
    else:
        norm_embeddings = normalize(embMat_filtered)
    
    if len(norm_embeddings) == 0:
        logger.error("No embeddings left after filtering. Exiting.")
        return

    # Fit HNSW index
    # ef is set to max_edges * 10 for better recall, as suggested in original code
    p = fit_hnsw_index(logger, norm_embeddings, ef=max_edges * 10)

    # Perform KNN query
    time_start = time.time()
    # ann_neighbor_indices: (n_samples, k+1) where k is max_edges. First column is self-index.
    # ann_distances: (n_samples, k+1) where k is max_edges. First column is self-distance (0).
    ann_neighbor_indices, ann_distances = p.knn_query(norm_embeddings, max_edges + 1)
    time_end = time.time()
    logger.info(f'KNN query time cost: {time_end - time_start:.2f}s')

    # Process ANN results to get sources, targets, and initial weights
    # Repeat source index for each neighbor
    sources = np.repeat(np.arange(len(norm_embeddings)), max_edges)
    # Exclude the first column (self-loop) from neighbors and distances
    targets_indices = ann_neighbor_indices[:, 1:]
    targets = targets_indices.flatten()
    wei = ann_distances[:, 1:]
    wei = wei.flatten()

    # Apply partgraph_ratio cutoff to filter edges by distance percentile
    dist_cutoff = np.percentile(wei, partgraph_ratio)
    save_index = wei <= dist_cutoff

    sources = sources[save_index]
    targets = targets[save_index]
    wei = wei[save_index]

    # Apply distance mode transformation (l1 or l2) to convert distances to weights
    # This transformation is applied here so the saved weights are ready for Leiden.
    if lmode == 'l1':
        wei = np.sqrt(wei) # Convert squared L2 distance to L1 distance if lmode is 'l1'
        wei = np.exp(-wei / bandwidth_for_edge_extraction)
    elif lmode == 'l2':
        wei = np.exp(-wei / bandwidth_for_edge_extraction)
    else:
        logger.warning(f"Unknown lmode '{lmode}'. Weights are not transformed based on lmode.")

    # Remove duplicate edges (e.g., (u,v) and (v,u)) and self-loops (u,u)
    # Ensure only one representation (u,v with u < v) is kept for each unique edge.
    unique_edges = set()
    final_sources = []
    final_targets = []
    final_weights = []

    for i in range(len(sources)):
        u, v = sources[i], targets[i]
        w = wei[i]
        
        # Skip self-loops
        if u == v:
            continue
        
        # Ensure a canonical representation (u,v where u < v)
        if u < v:
            edge_tuple = (u, v)
        else:
            edge_tuple = (v, u)
        
        # Add the edge if it hasn't been added yet
        if edge_tuple not in unique_edges:
            unique_edges.add(edge_tuple)
            final_sources.append(u)
            final_targets.append(v) # Store as (u,v) or (v,u) as they came, but the tuple is canonical
            final_weights.append(w)
        # If edge_tuple is already in unique_edges, it means we've seen (u,v) or (v,u) before.
        # We keep the first one encountered. If weights could differ for (u,v) and (v,u) and
        # you need to average/sum, more complex logic would be needed.

    final_sources = np.array(final_sources)
    final_targets = np.array(final_targets)
    final_weights = np.array(final_weights)

    logger.info(f"Number of extracted edges (after filtering and deduplication): {len(final_sources)}")

    # Save extracted edges, namelist, and length_weight to the output directory
    np.savez(os.path.join(output_dir, 'extracted_edges.npz'),
             sources=final_sources, targets=final_targets, weights=final_weights)
    
    pd.DataFrame(namelist_filtered).to_csv(os.path.join(output_dir, 'namelist.txt'), index=False, header=False)
    pd.DataFrame(length_weight_filtered).to_csv(os.path.join(output_dir, 'length_weight.txt'), index=False, header=False)

    logger.info(f"Extracted edges, namelist, and length_weight saved to {output_dir}")
    logger.info("Edge extraction complete.")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Step 1: Extract edges and weights for clustering.')
    parser.add_argument('--emb_file', type=str, required=True, help='Path to the embeddings file.')
    parser.add_argument('--contig_file', type=str, required=True, help='Path to the contig file (for lengths).')
    parser.add_argument('--output_dir', type=str, required=True, 
                        help='Output directory to save extracted edges (extracted_edges.npz), namelist (namelist.txt), and length_weight (length_weight.txt).')
    parser.add_argument('--contig_len', type=int, default=1000, 
                        help='Minimum contig length to consider for clustering.')
    parser.add_argument('--not_l2normalize', action='store_true', 
                        help='Do not perform L2 normalization on embeddings before HNSW. By default, embeddings are L2 normalized.')
    parser.add_argument('--max_edges', type=int, default=100, 
                        help='Maximum number of neighbors to consider for HNSW query (k in k-NN).')
    parser.add_argument('--partgraph_ratio', type=int, default=50, 
                        help='Percentile cutoff for edge distances (e.g., 50 means keep edges with distances up to the 50th percentile).')
    parser.add_argument('--lmode', type=str, default='l2', choices=['l1', 'l2'], 
                        help='Distance mode for converting distances to weights (l1 or l2).')
    parser.add_argument('--bandwidth_for_edge_extraction', type=float, default=0.1, 
                        help='Bandwidth parameter used in the exponential transformation of distances to weights during edge extraction.')

    args = parser.parse_args()
    extract_edges(logger, args)
