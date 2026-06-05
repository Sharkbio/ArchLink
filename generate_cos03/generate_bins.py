import pandas as pd
import os
import shutil
import sys

def filter_and_copy_bins(summary_file, source_fa_dir, output_dir):
    """
    Reads a summary file, filters for high-purity bins, and copies
    the corresponding FASTA files to an output directory.

    Args:
        summary_file (str): Path to the summary TSV file.
        source_fa_dir (str): Path to the directory containing source FASTA files.
        output_dir (str): Path to the output directory.
    """
    try:
        # Check if the summary file exists
        if not os.path.exists(summary_file):
            print(f"Error: Summary file '{summary_file}' not found.")
            return

        # Create output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)
        
        # Read the summary file using pandas
        df = pd.read_csv(summary_file, sep='\t')
        
        # Check for required columns
        required_cols = ['Name', 'Contamination', 'Genome_Size', 'Total_Contigs', 'Completeness']
        if not all(col in df.columns for col in required_cols):
            print(f"Error: The summary file must contain columns: {required_cols}")
            return

        # Filter the DataFrame for bins with Contamination < 10%
        high_purity_bins = df[df['Contamination'] < 10.0].copy()
        
        # Calculate 'Purity' from 'Contamination'
        high_purity_bins['Purity'] = 100 - high_purity_bins['Contamination']
        
        # Select and rename columns for the output summary file
        output_df = high_purity_bins[['Name', 'Total_Contigs', 'Genome_Size', 'Purity', 'Completeness']].copy()
        output_df.rename(columns={
            'Name': 'Bin_Name',
            'Total_Contigs': 'Contig_Count',
            'Genome_Size': 'Total_Size_bp',
            'Purity': 'Precision',
            'Completeness': 'Recall'
        }, inplace=True)
        
        # Define the path for the output summary file
        output_summary_path = os.path.join(output_dir, "high_purity_bins_summary.tsv")
        output_df.to_csv(output_summary_path, sep='\t', index=False)
        print(f"Successfully filtered {len(high_purity_bins)} high-purity bins.")
        print(f"High-purity bins summary saved to: {output_summary_path}")

        # Create individual subdirectories and copy the corresponding FASTA files
        print("Creating subdirectories and copying high-purity bin FASTA files...")
        
        copied_count = 0
        for bin_name in high_purity_bins['Name']:
            # Create a new subdirectory for each bin
            bin_dir = os.path.join(output_dir, str(bin_name))
            os.makedirs(bin_dir, exist_ok=True)
            
            source_file = os.path.join(source_fa_dir, f"{bin_name}.fa")
            destination_file = os.path.join(bin_dir, f"{bin_name}.fasta")

            if os.path.exists(source_file):
                shutil.copyfile(source_file, destination_file)
                copied_count += 1
                print(f"  - Copied '{os.path.basename(source_file)}' to '{os.path.basename(bin_dir)}'")
            else:
                print(f"  - Warning: Source file '{source_file}' not found. Skipping.")

        print(f"Finished copying {copied_count} FASTA files.")
        
    except FileNotFoundError as e:
        print(f"Error: A file or directory was not found. Please check your paths. Details: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

def main():
    if len(sys.argv) != 4:
        print("Usage: python filter_bins.py <summary_file> <source_fa_dir> <output_dir>")
        print("Example: python filter_bins.py bin_summary.tsv ./all_bins ./high_purity_results")
        sys.exit(1)

    summary_file = sys.argv[1]
    source_fa_dir = sys.argv[2]
    output_dir = sys.argv[3]

    filter_and_copy_bins(summary_file, source_fa_dir, output_dir)

if __name__ == "__main__":
    main()
