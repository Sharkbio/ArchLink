import pandas as pd
import re
import os
import sys
import shutil
from typing import List


def get_unique_top_methods(df: pd.DataFrame) -> List[str]:
    """
    Select the top 3 method IDs based purely on the 'sum' column. 
    'sum_cont5' is used as a tie-breaker if 'sum' scores are equal.

    Args:
        df (pd.DataFrame): Must contain 'sum', 'sum_cont5', and 'Binning_method'.

    Returns:
        List[str]: A list of up to 3 unique binning method IDs. 
                   Rank 1 (index 0) is considered the 'best'.
    """

    # 1. Sort by sum (primary, descending) and sum_cont5 (secondary, descending for tie-breaking)
    df_sorted = df.sort_values(by=['sum', 'sum_cont5'], ascending=[False, False])

    # 2. Select the top 3 method IDs
    # 提取前 3 个独特的方法 ID
    final_methods = df_sorted.head(3)['Binning_method'].tolist()

    return final_methods


def process_selected_binning_methods(b_val: str):
    """
    Workflow:
      1. Read estimate_res.txt
      2. Select best method IDs (Top 3 by sum, with the first one being the 'best')
      3. Create output folders
      4. Copy result files and swap the first two columns

    Args:
        b_val (str): Base directory name (e.g., 'urogenital').
    """

    # --- Step 1: Read estimate_res.txt ---

    estimate_file_path = os.path.expanduser(f"{b_val}/cluster_res/estimate_res.txt")
    source_base_dir = os.path.expanduser(f"{b_val}/cluster_res/")
    output_base_dir = os.path.expanduser(f"{b_val}/binning/frequency")


    os.makedirs(output_base_dir, exist_ok=True)

    if not os.path.exists(estimate_file_path):
        print(f"ERROR: estimate_res.txt not found - {estimate_file_path}", file=sys.stderr)
        sys.exit(1)

    try:
        df = pd.read_csv(estimate_file_path, sep='\t', low_memory=False)
    except Exception as e:
        print(f"ERROR: Error reading {estimate_file_path}: {e}", file=sys.stderr)
        sys.exit(1)

    required_cols = ['sum', 'sum_cont5', 'Binning_method']
    if not all(col in df.columns for col in required_cols):
        print(f"ERROR: Required columns {required_cols} missing.", file=sys.stderr)
        sys.exit(1)

    # --- Step 2: Select IDs using the independent rules (Top 3 by sum) ---
    # 现在提取前 3 个方法
    final_methods: List[str] = get_unique_top_methods(df)

    if len(final_methods) == 0:
        print("ERROR: No binning method IDs selected. Data may be insufficient.", file=sys.stderr)
        sys.exit(1)
    
    # 打印实际选择的数量（最多 3 个）
    print(f"INFO: Selected {len(final_methods)} methods (max 3): {final_methods}", file=sys.stderr)

    # --- Step 3: Create output folders and process files ---

    pattern = re.compile(
        r'Leiden_bandwidth_(\d+\.\d+)_res_maxedges(\d+)respara_(\d+)_partgraph_ratio_(\d+)\.tsv'
    )

    for i, original_id in enumerate(final_methods):

        # Convert long ID to a shorter format
        match = pattern.search(original_id)
        if match:
            converted_id = f"{match.group(1)}_{match.group(2)}_{match.group(3)}_{match.group(4)}"
        else:
            print(f"WARNING: Could not parse '{original_id}', using simplified form.", file=sys.stderr)
            converted_id = (
                original_id.replace('.tsv', '').replace('Leiden_bandwidth_', '')
            )

        # Mark the best (rank 1, which is the first element) with a suffix
        is_best = (i == 0)
        folder_name = converted_id + ('_best' if is_best else '')

        if is_best:
            print(" Best method identified (Rank 1 by sum).", file=sys.stderr)

        print(f" Original ID : {original_id}", file=sys.stderr)
        print(f" Folder name : {folder_name}", file=sys.stderr)

        target_folder = os.path.join(output_base_dir, folder_name)
        os.makedirs(target_folder, exist_ok=True)

        source_file = os.path.join(source_base_dir, original_id)
        destination_file = os.path.join(target_folder, "result")

        # Copy and process the file
        if os.path.exists(source_file):
            try:
                # Use standard copy
                shutil.copyfile(source_file, destination_file)
                print(f" Copied to: {destination_file}", file=sys.stderr)

                # Swap the first two columns (Contig ID <-> Bin ID)
                lines = []
                with open(destination_file, "r") as f:
                    lines = f.readlines()

                with open(destination_file, "w") as f:
                    for line in lines:
                        parts = line.strip().split('\t')
                        if len(parts) >= 2:
                            # Swap operation: parts[0] and parts[1]
                            parts[0], parts[1] = parts[1], parts[0]
                        f.write('\t'.join(parts) + '\n')

                print(f" Swapped first two columns in '{destination_file}'.", file=sys.stderr)

            except Exception as e:
                print(f"ERROR: File processing failed: {e}", file=sys.stderr)

        else:
            print(f"ERROR: Source file not found - {source_file}", file=sys.stderr)

    print(f"\nSUCCESS: Processing complete. Output: {output_base_dir}", file=sys.stderr)


if __name__ == "__main__":
    # Example usage: process_selected_binning_methods('./path/to/results')
    # This requires running it with an argument in a real environment.
    # For standalone test:
    # process_selected_binning_methods(sys.argv[1] if len(sys.argv) > 1 else 'test_dir')
    pass