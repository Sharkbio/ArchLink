import os
import logging
import argparse
import sys
from collections import defaultdict
import pandas as pd
import glob

from .filter_small_bins import filter_small_bins
import shutil # 导入 shutil 模块用于文件复制
import os, sys
PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

# 放入 sys.path
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)
# Set up logging configuration
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
from scripts.unitem_markers import Markers
from scripts.unitem_profile import Profile

def get_binstats(bin_contig_names, markers):
    """
    Calculates the completion and contamination of a single bin.
    :param bin_contig_names: A set of contig names in the bin.
    :param markers: An instance of the Markers object for bin quality calculation.
    :return: A tuple of completion and contamination.
    """
    _, comp, cont = markers.bin_quality(bin_contig_names)
    return comp, cont

def read_bins_nosequences(bin_dirs):
    """
    Reads binning results from files.
    - Assumes the input files are two-column TSV files: bin_id in the first column, contig_id in the second.
    :param bin_dirs: A dictionary where keys are method IDs and values are paths to the binning result files.
    :return: A dictionary containing bin data and a dictionary of contigs' assigned bins.
    """
    bins = defaultdict(lambda: defaultdict(set))
    contigs_in_bins = defaultdict(lambda: {})

    for method_id, bin_path in bin_dirs.items():
        try:
            # Modified to expect: Column 0 = bin_id, Column 1 = contig_id
            df = pd.read_csv(bin_path, sep='\t', header=None, names=['bin', 'contig'])
            
            for _, row in df.iterrows():
                # Extract data based on the new column names
                bin_id = str(row['bin']).strip()       # Column 0: bin_id
                contig = str(row['contig']).strip()    # Column 1: contig_id
                
                if contig and bin_id:  # Ensure values are not empty
                    bins[method_id][bin_id].add(contig)
                    contigs_in_bins[contig][method_id] = bin_id
        except Exception as e:
            logging.error(f"Error reading binning file {bin_path}: {str(e)}", exc_info=True)
            continue

    return bins, contigs_in_bins

def get_bin_quality(orig_bins, methods_sorted, markers):
    """
    Calculates quality metrics for all binning methods and determines the best method.
    :param orig_bins: Dictionary of original bin data.
    :param methods_sorted: A sorted list of method IDs.
    :param markers: An instance of the Markers object.
    :return: A dictionary with quality metrics for all bins and the ID of the best method.
    """
    bin_quality_dict = {}
    sum_list = []
    sumcont5_list = []

    for method_id in methods_sorted:
        if method_id not in orig_bins:
            logging.warning(f"No bins found for method: {method_id}")
            continue
            
        num_5010, num_7010, num_9010 = 0, 0, 0
        num_505, num_705, num_905 = 0, 0, 0

        for bin_id, contig_set in orig_bins[method_id].items():
            comp, cont = get_binstats(contig_set, markers)
            if comp > 50 and cont < 10: num_5010 += 1
            if comp > 70 and cont < 10: num_7010 += 1
            if comp > 90 and cont < 10: num_9010 += 1
            if comp > 50 and cont < 5: num_505 += 1
            if comp > 70 and cont < 5: num_705 += 1
            if comp > 90 and cont < 5: num_905 += 1

        metrics = {
            'num_5010': num_5010, 'num_7010': num_7010, 'num_9010': num_9010,
            'num_505': num_505, 'num_705': num_705, 'num_905': num_905,
            'sum': num_5010 + num_7010 + num_9010 + num_505 + num_705 + num_905,
            'sum_cont5': num_505 + num_705 + num_905
        }
        bin_quality_dict[method_id] = metrics
        sum_list.append(metrics['sum'])
        sumcont5_list.append(metrics['sum_cont5'])

    if not sum_list:
        raise ValueError("No bins found for any method to evaluate quality.")
        
    sum_max = max(sum_list)
    best_candidates = [method_id for method_id, q in bin_quality_dict.items() if q['sum'] == sum_max]
    
    if not best_candidates:
        raise ValueError("No candidate for best method found based on quality sum.")
        
    if len(best_candidates) > 1:
        best_cont5 = max([bin_quality_dict[method]['sum_cont5'] for method in best_candidates])
        # If multiple methods have the same best_cont5, choose the first one
        best_method = [method for method in best_candidates if bin_quality_dict[method]['sum_cont5'] == best_cont5][0]
    else:
        best_method = best_candidates[0]

    return bin_quality_dict, best_method

def save_high_quality_contigs(orig_bins, best_method, markers, outpath):
    """
    Saves high-quality contigs from the best binning result to a file.
    Creates two files: one for bins with >50% completion and <10% contamination,
    and another for bins with >50% completion and <5% contamination.
    :param orig_bins: Dictionary of original bin data.
    :param best_method: ID of the best method.
    :param markers: An instance of the Markers object.
    :param outpath: Output directory.
    """
    bin_count_5010, bin_count_5005 = 0, 0
    os.makedirs(outpath, exist_ok=True) # Ensure output directory exists
    
    f1_path = os.path.join(outpath, f"{best_method}_5010_contigs.tsv")
    f2_path = os.path.join(outpath, f"{best_method}_5005_contigs.tsv")
    
    if best_method not in orig_bins:
        logging.error(f"Best method '{best_method}' not found in bin data. Cannot save high-quality contigs.")
        return
        
    with open(f1_path, "w") as f1, open(f2_path, "w") as f2:
        for bin_id, contig_set in orig_bins[best_method].items():
            comp, cont = get_binstats(contig_set, markers)
            
            # Save bins with >50% completion and <10% contamination
            if comp > 50 and cont < 10:
                for contig in contig_set:
                    # Output format is: contig_id \t bin_id (This remains the same for consistency with downstream tools)
                    f1.write(f"{bin_id}\t{contig}\n") 
                bin_count_5010 += 1
                
            # Save bins with >50% completion and <5% contamination
            if comp > 50 and cont < 5:
                for contig in contig_set:
                    # Output format is: contig_id \t bin_id
                    f2.write(f"{bin_id}\t{contig}\n") 
                bin_count_5005 += 1
                
    logging.info(f"Saved {bin_count_5010} high-quality bins to {f1_path} (completion > 50%, contamination < 10%)")
    logging.info(f"Saved {bin_count_5005} high-quality bins to {f2_path} (completion > 50%, contamination < 5%)")

def write_quality_summary(bin_quality_dict, output_path):
    """
    Writes a quality summary of all methods to a file.
    :param bin_quality_dict: Dictionary containing quality metrics for all binning methods.
    :param output_path: Output directory.
    """
    output_file = os.path.join(output_path, "bin_quality_summary.tsv")
    with open(output_file, "w") as fout:
        fout.write("Method\t>50%<10%\t>70%<10%\t>90%<10%\t>50%<5%\t>70%<5%\t>90%<5%\tTotal\tTotal<5%\n")
        for method_id, metrics in bin_quality_dict.items():
            fout.write(f"{method_id}\t{metrics['num_5010']}\t{metrics['num_7010']}\t"
                        f"{metrics['num_9010']}\t{metrics['num_505']}\t{metrics['num_705']}\t"
                        f"{metrics['num_905']}\t{metrics['sum']}\t{metrics['sum_cont5']}\n")
    logging.info(f"Quality summary saved to {output_file}")

def read_fasta(filepath):
    """
    A simple FASTA file parser.
    :param filepath: Path to the FASTA file.
    :return: A dictionary containing contig IDs and sequences.
    """
    sequences = {}
    current_id = None
    current_seq_parts = []
    try:
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    if current_id:
                        sequences[current_id] = ''.join(current_seq_parts)
                    current_id = line[1:].split()[0] # Take only the first word as ID
                    current_seq_parts = []
                else:
                    current_seq_parts.append(line)
        if current_id: # Process the last contig
            sequences[current_id] = ''.join(current_seq_parts)
    except Exception as e:
        logging.error(f"Error reading FASTA file {filepath}: {str(e)}", exc_info=True)
        sys.exit(1)
    return sequences

def estimate_bins_quality(res_path, contig_file, bac_mg_table=None, ar_mg_table=None):
    """
    Main function for evaluating Leiden binning results.
    If marker gene tables are not provided, it will run unitem_profile to generate them.
    :param res_path: Directory containing Leiden result files.
    :param contig_file: Path to the contig FASTA file.
    :param bac_mg_table: Path to the bacterial marker gene stats table (optional).
    :param ar_mg_table: Path to the archaeal marker gene stats table (optional).
    :return: The ID of the selected best method.
    """
    markers = Markers()
    
    # Find all files starting with "Leiden" and ending with ".tsv"
    leiden_files = [f for f in os.listdir(res_path) if f.startswith("Leiden") and f.endswith(".tsv")]
    if not leiden_files:
        raise FileNotFoundError(f"No Leiden binning result files found in {res_path}. "
                                 f"Expected files starting with 'Leiden' and ending with '.tsv'.")
    
    # Build the bin_dirs dictionary, where keys are method IDs (filenames without extension) and values are full paths
    bin_dirs = {os.path.splitext(fname)[0]: os.path.join(res_path, fname) for fname in leiden_files}
    
    # Read binning results (contig-to-bin mapping)
    bins, _ = read_bins_nosequences(bin_dirs)
    # bins is a dict: {method_id: {bin_id: set(contig_ids)}}
    
    methods_sorted = sorted(bins.keys()) # Get a sorted list of all Leiden methods
    
    # Check if valid marker gene tables are provided, if not, generate them
    if not (bac_mg_table and os.path.exists(bac_mg_table) and ar_mg_table and os.path.exists(ar_mg_table)):
        logging.info("Bacterial or archaeal marker gene tables not provided or not found. "
                     "Running unitem_profile to generate them...")
        
        if not contig_file or not os.path.exists(contig_file):
            raise FileNotFoundError("Contig FASTA file is required to generate marker gene tables but was not provided or found.")
        
        logging.info(f"Reading contigs from {contig_file}...")
        contig_sequences = read_fasta(contig_file)
        logging.info(f"Read {len(contig_sequences)} contigs.")

        # Create a temporary output directory for unitem_profile
        profile_output_dir = os.path.join(res_path, "unitem_profile_temp_output")
        os.makedirs(profile_output_dir, exist_ok=True)
        
        # Select the first Leiden method to generate marker gene tables
        # unitem_profile requires binned FASTA files as input
        sample_method_for_profile = methods_sorted[0]

        temp_bin_fasta_dir = os.path.join(profile_output_dir, f"binned_fasta_for_{sample_method_for_profile}")
        os.makedirs(temp_bin_fasta_dir, exist_ok=True)
        
        logging.info(f"Writing binned FASTA files for method '{sample_method_for_profile}' to {temp_bin_fasta_dir} for profiling...")
        profile_input_bin_dirs = {}
        # Create a FASTA file for each bin of the selected_method
        for bin_id, contig_ids in bins[sample_method_for_profile].items():
            bin_fasta_path = os.path.join(temp_bin_fasta_dir, f"{bin_id}.fa")
            # print(bin_fasta_path)

            with open(bin_fasta_path, 'w') as f_out:
                
                for contig_id in contig_ids:
                    if contig_id in contig_sequences:

                        f_out.write(f">{contig_id}\n{contig_sequences[contig_id]}\n")
                    else:
                        logging.warning(f"Contig '{contig_id}' from bin '{bin_id}' of method '{sample_method_for_profile}' not found in contig file. Skipping.")
        
        profile_input_bin_dirs[sample_method_for_profile] = (temp_bin_fasta_dir, 'fa')
        # Run unitem_profile
        profiles = Profile(cpus=20) # Number of threads can be adjusted
        profiles.run(profile_input_bin_dirs, profile_output_dir)
        
        # Get the paths of the generated marker gene tables
        bac_mg_table = os.path.join(profile_output_dir, "binning_methods", sample_method_for_profile, "checkm_bac", "marker_gene_table.tsv")
        ar_mg_table = os.path.join(profile_output_dir, "binning_methods", sample_method_for_profile, "checkm_ar", "marker_gene_table.tsv")
        
        # Verify that the marker gene tables were generated
        if not (os.path.exists(bac_mg_table) and os.path.exists(ar_mg_table)):
            raise FileNotFoundError(f"Marker gene tables could not be generated or found at expected paths: "
                                     f"'{bac_mg_table}' and '{ar_mg_table}'. "
                                     "Please check unitem_profile output.")
        logging.info("Marker gene tables generated successfully.")
    else:
        logging.info("Using provided marker gene tables for quality assessment.")

    # Load marker gene tables into the Markers object
    markers.marker_gene_tables(bac_mg_table, ar_mg_table)
    
    # Evaluate the quality of all Leiden methods
    bin_quality_dict, best_method = get_bin_quality(bins, methods_sorted, markers)
    
    # Save high-quality contigs
    save_high_quality_contigs(bins, best_method, markers, res_path)
    
    # Write quality summary
    write_quality_summary(bin_quality_dict, res_path)
    
    logging.info(f"Best Leiden method selected: {best_method}")
    return best_method

def copy_best_result(best_method_path):
    """
    Copies the best result file and appends '_best' to its name.
    :param best_method_path: The full path to the best result file.
    """
    try:
        # Split path into root and extension
        root, ext = os.path.splitext(best_method_path)
        # Create the new path with '_best' appended
        best_path = f"{root}_best{ext}"
        
        shutil.copyfile(best_method_path, best_path)
        logging.info(f"Successfully copied the best result file to: {best_path}")
        return best_path
    except FileNotFoundError:
        logging.error(f"Error: The file '{best_method_path}' was not found. Cannot copy.")
        return None
    except Exception as e:
        logging.error(f"An error occurred while copying the file: {e}")
        return None

def main(args):
    
    args.__setattr__('res_path', os.path.join(args.output_path, 'binning/s_cluster/cluster_res/'))
    
    # 查找匹配的文件夹
    binning_dir = os.path.join(args.output_path, 'cluster_res/unitem_profile', 'binning_methods')
    
    # 使用glob查找匹配的文件
    pattern = os.path.join(binning_dir, 'weight_seed_kmeans_k_*_result.tsv')
    matching_folders = glob.glob(pattern)
    
    if not matching_folders:
        # 如果没有找到，尝试更灵活的模式
        pattern = os.path.join(binning_dir, 'weight_seed_kmeans_*_result.tsv')
        matching_folders = glob.glob(pattern)
    
    if matching_folders:
        # 获取第一个匹配的文件夹
        matched_folder = matching_folders[0]
        
        # 设置bac_mg_table和ar_mg_table路径
        bac_mg_path = os.path.join(matched_folder, 'checkm_bac', 'marker_gene_table.tsv')
        ar_mg_path = os.path.join(matched_folder, 'checkm_ar', 'marker_gene_table.tsv')
        
        # 检查路径是否存在
        if not os.path.exists(bac_mg_path):
            print(f"警告: 文件不存在: {bac_mg_path}")
        if not os.path.exists(ar_mg_path):
            print(f"警告: 文件不存在: {ar_mg_path}")
        
        args.__setattr__('bac_mg_table', bac_mg_path)
        args.__setattr__('ar_mg_table', ar_mg_path)
        
        print(f"已自动匹配文件夹: {os.path.basename(matched_folder)}")    
    
    
    log_file_path = os.path.join(args.res_path, "quality_assessment.log")
    # Reconfigure logging to output to both a file and the console
    logging.getLogger().handlers = [] # Clear old handlers
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s",
                         handlers=[logging.FileHandler(log_file_path), logging.StreamHandler(sys.stdout)])
    
    logging.info(f"Starting quality assessment for Leiden results in: {args.res_path}")
    logging.info(f"Using contig file: {args.contig_file}")
    if args.bac_mg_table and args.ar_mg_table:
        logging.info(f"Using provided bacterial marker table: {args.bac_mg_table}")
        logging.info(f"Using provided archaeal marker table: {args.ar_mg_table}")
    else:
        logging.info("Marker tables not provided, they will be generated using unitem_profile.")
    try:
        # Call the evaluation function
        best_method = estimate_bins_quality(
            res_path=args.res_path,
            contig_file=args.contig_file,
            bac_mg_table=args.bac_mg_table,
            ar_mg_table=args.ar_mg_table
        )
        
        # Get the original file path for the best method
        original_leiden_file_path = os.path.join(args.res_path, f"{best_method}.tsv")

        # Copy the best result file
        copy_best_result(original_leiden_file_path)

        logging.info("Filtering small bins based on contig lengths.")
        
        leiden_files_for_filter = [f for f in os.listdir(args.res_path) if f.startswith("Leiden") and f.endswith(".tsv")]
        if not leiden_files_for_filter:
            logging.error(f"No Leiden binning result files found in {args.res_path} for filtering. Cannot perform small bin filtering.")
            sys.exit(1) 
        
        original_bin_dirs_for_filter = {os.path.splitext(fname)[0]: os.path.join(args.res_path, fname) for fname in leiden_files_for_filter}
        original_leiden_file_for_filtering = original_bin_dirs_for_filter.get(best_method)

        if original_leiden_file_for_filtering and os.path.exists(original_leiden_file_for_filtering):
            logging.info(f"Applying small bin filtering to the original Leiden result: '{original_leiden_file_for_filtering}'")
            
            if not hasattr(args, 'output_path'):
                args.output_path = args.res_path
            else:
                if args.output_path != args.res_path:
                    logging.warning(f"args.output_path was already set to '{args.output_path}'. "
                                    f"Overriding with args.res_path: '{args.res_path}' for filter_small_bins compatibility.")
                    args.output_path = args.res_path


            filter_small_bins(
                logging.getLogger(),  # logger object
                args.contig_file,     # contig FASTA file path
                original_leiden_file_for_filtering, # Path to the original binning result file of the best Leiden method
                args              # argparse.Namespace object
            )
            logging.info(f"Small bin filtering complete for '{best_method}'. "
                         f"The output file path depends on the internal logic of filter_small_bins, usually saved in the '{args.output_path}' directory with a name like '{best_method}_filtered.tsv'.")
        else:
            logging.warning(f"Original Leiden file for best method '{best_method}' not found at '{original_leiden_file_for_filtering}'. "
                            "Skipping small bin filtering based on contig lengths.")

        logging.info("Assessment complete!")
        logging.info(f"Log file saved to: {log_file_path}")

    except Exception as e:
        logging.exception("An error occurred during processing. Please check the error message and log file for details.")
        sys.exit(1)

if __name__ == "__main__":
    main()