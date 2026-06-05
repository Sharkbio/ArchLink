import sys
import itertools
from collections import defaultdict
import pickle
import os
import glob
import concurrent.futures
import threading

# Global lock to protect concurrent access to string_pool
string_pool_lock = threading.Lock()
# Global string pool for memory optimization
string_pool = {}

def get_cached_string_thread_safe(s):
    """
    Thread-safe string caching function.
    Ensures consistency when multiple threads access string_pool concurrently.
    """
    with string_pool_lock:
        return string_pool.setdefault(s, s)

def _process_single_file_for_parallel(file_path, num_files_total, file_idx):
    """
    Helper function to process a single cluster result file in a separate thread.
    Reads the file, builds a local inverted index for co-occurring contigs, and returns it.
    
    Args:
        file_path (str): Path to the cluster result file.
        num_files_total (int): Total number of files being processed.
        file_idx (int): Index of the current file (0-based) for progress printing.
        
    Returns:
        defaultdict: Local inverted index {contig_A: {contig_B: count}} for this file.
    """
    print(f"INFO: Processing file {file_idx+1}/{num_files_total}: {os.path.basename(file_path)}")
    bins = defaultdict(list)
    inverted_index_partial = defaultdict(lambda: defaultdict(int)) # Local inverted index

    # 1. Group contigs by bin_id
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                # Use thread-safe string caching
                bin_id = get_cached_string_thread_safe(parts[0])
                contig = get_cached_string_thread_safe(parts[1])
                bins[bin_id].append(contig)
    
    # 2. Generate all pairwise edges (co-occurrences) within each bin
    for contigs in bins.values():
        if len(contigs) > 1:
            # Sort and unique the contigs to ensure canonical edge definition
            sorted_contigs = sorted(set(contigs))
            for c1, c2 in itertools.combinations(sorted_contigs, 2):
                # Store the edge in a canonical order (smaller_contig, larger_contig)
                # The stored count is 1, indicating this edge occurred in this file.
                ordered_pair = tuple(sorted((c1, c2)))
                inverted_index_partial[ordered_pair[0]][ordered_pair[1]] = 1
    
    return inverted_index_partial

def process_files_and_generate_dict(input_files, output_file, num_threads=32):
    """
    Processes cluster result files, counts edge co-occurrences, calculates frequencies,
    and saves the resulting frequency index.
    
    Args:
        input_files (list): List of paths to individual cluster result files.
        output_file (str): Output pickle file path for the frequency index.
        num_threads (int): Number of threads to use for parallel file processing.
        
    Returns:
        dict: Nested dictionary frequency index {contig: {neighbor: frequency}}.
    """
    global string_pool

    # Clear the global string pool for each call to ensure a clean state
    string_pool.clear()

    # Final inverted index to store total co-occurrence counts across all files
    inverted_index = defaultdict(lambda: defaultdict(int))
    num_files = len(input_files)

    print(f"--- Step 1: Counting edge co-occurrences (Total files: {num_files}, Threads: {num_threads}) ---")
    
    tasks_args = [(file_path, num_files, i) for i, file_path in enumerate(input_files)]
    # Use ThreadPoolExecutor for parallel file processing
    with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as executor:
        futures = [executor.submit(_process_single_file_for_parallel, *args) for args in tasks_args]
        
        # Merge results from all threads
        for future in concurrent.futures.as_completed(futures):
            try:
                inverted_index_partial = future.result()
                # Merge the local inverted index into the main inverted index
                for c1, neighbors in inverted_index_partial.items():
                    for c2, count in neighbors.items():
                        # Count the number of files the edge appeared in (symmetric storage)
                        inverted_index[c1][c2] += count
                        inverted_index[c2][c1] += count # Add symmetric entry
            except Exception as exc:
                print(f'ERROR: File processing generated an exception: {exc}')
    
    print(f"INFO: Inverted index built, containing {len(inverted_index)} unique contigs with edges.")

    print("\n--- Step 2: Converting counts to normalized frequency ---")
    final_freq_index = defaultdict(dict)   
    # Calculate frequencies by dividing total count by the number of files
    for c1, neighbors in inverted_index.items():
        for c2, count in neighbors.items():
            frequency = count / num_files
            final_freq_index[c1][c2] = frequency

    # Save the optimized dictionary index file
    with open(output_file, 'wb') as pf:
        pickle.dump(final_freq_index, pf, protocol=4)

    print(f"INFO: Optimized dictionary index saved to {output_file}")
    print(f"INFO: Total unique contigs in frequency index: {len(final_freq_index)}")
    return final_freq_index

class EdgeFrequencyIndexDict:
    """Efficient edge frequency query using a dictionary inverted index."""
    def __init__(self, index_path):
        with open(index_path, 'rb') as f:
            self.index = pickle.load(f)
        print(f"INFO: Loaded dictionary index with {len(self.index)} contigs.")
    
    def query(self, contig1, contig2):
        """Query the frequency of an edge between two contigs."""
        # Check for both canonical orders due to symmetric storage
        if contig1 in self.index and contig2 in self.index[contig1]:
            return self.index[contig1][contig2]
        return 0.0
    
    def batch_query(self, edge_list):
        """Batch query frequencies for multiple edges."""
        results = []
        for a, b in edge_list:
            results.append(self.query(a, b))
        return results

import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Processes cluster result files within a directory to generate a frequency index.")
    # The script now expects a single directory path
    parser.add_argument('input_directory', 
                        help='Path to the base directory containing cluster result subfolders (e.g., /path/to/results/{B_VAL}). The script searches for files matching "*/result" within this directory.')
    parser.add_argument('output_file', 
                        help='Output pickle file path for the frequency index (e.g., /path/to/results/{B_VAL}/frequency.pkl).')
    parser.add_argument('--num_threads', type=int, default=32, 
                        help='Number of threads for parallel processing (default: 32).')
    
    args = parser.parse_args()
    
    # New logic: Find all 'result' files within subdirectories of input_directory
    result_files_pattern = os.path.join(args.input_directory, "*", "result")
    input_files = glob.glob(result_files_pattern)
    
    if not input_files:
        print(f"ERROR: No cluster result files found in {args.input_directory} using pattern '*/result'. Aborting.")
        sys.exit(1)

    print(f"INFO: Found {len(input_files)} cluster result files for processing.")

    # Ensure output directory exists (if not the current directory)
    output_dir = os.path.dirname(args.output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Call the processing function with thread count
    process_files_and_generate_dict(input_files, args.output_file, args.num_threads)
    print("INFO: Processing complete. Frequency data successfully saved.")