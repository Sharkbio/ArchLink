import sys
from pathlib import Path
import os
def process_single_file_by_bin(file_path, output_dir_path):
    """
    Reads a single tab-separated file, processes the data by binID, and
    generates a connect_graph.txt file for each bin.

    The input file format is expected to be:
    contig1  contig1_direction  contig2  contig2_direction  weight  binID
    
    Args:
        file_path (str): Path to the single input file.
        output_dir_path (str): The base directory to create sub-directories for each bin.
    """
    
    # Dictionary to store SEG and JUNC lines for each bin
    # The key is the binID, the value is another dictionary with 'seg' and 'junc' lists
    bin_data = {}

    try:
        with open(file_path, 'r') as input_file:
            for line_number, line in enumerate(input_file):
                line = line.strip()
                if not line or line.startswith('#'):  # Skip empty or commented lines
                    continue

                parts = line.split('\t')
                
                # Check for the correct number of columns
                if len(parts) < 6:
                    print(f"warn:  {line_number + 1} error, skip: {line}")
                    continue
                
                # Parse the data from the new file format
                contig1, dir1, contig2, dir2, weight_str, bin_id = parts[0], parts[1], parts[2], parts[3], parts[4], parts[5]
                
                try:
                    # Convert weight string to a float
                    weight = float(weight_str)
                except ValueError:
                    print(f"warn:  {line_number + 1} not skip: {weight_str}")
                    continue

                # Initialize data structures for the bin if it doesn't exist
                if bin_id not in bin_data:
                    bin_data[bin_id] = {'seg': set(), 'junc': []}

                # Add contigs to the SEG list for the current bin. Using a set
                # to automatically handle duplicates.
                bin_data[bin_id]['seg'].add(contig1)
                bin_data[bin_id]['seg'].add(contig2)

                # Format the JUNC line and add it to the JUNC list
                junc_line = f"JUNC\t{contig1}\t{dir1}\t{contig2}\t{dir2}\t{weight}"
                bin_data[bin_id]['junc'].append(junc_line)

    except FileNotFoundError:
        print(f"ERROR: not found: {file_path}")
        return
    except Exception as e:
        print(f"ERROE: {e}")
        return

    print("Writing file...")
    
    # Write the SEG and JUNC lines to their respective bin files
    for bin_id, data in bin_data.items():
        # Create the full output path
        output_path = os.path.join(output_dir_path, bin_id, "connect_graph_dir.txt")
        
        # Ensure the sub-directory exists
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        
        try:
            with open(output_path, 'w') as out_file:
                # First, write all SEG lines using the graph format expected downstream.
                for contig in sorted(list(data['seg'])):
                    out_file.write(f"SEG\t{contig}\t*\n")
                
                # Then, write all JUNC lines
                for junc_line in data['junc']:
                    out_file.write(junc_line + '\n')
            print(f"已为 binID '{bin_id}' 保存图文件至 {output_path}")
        except IOError as e:
            print(f"错误: 无法写入输出文件 {output_path}. 错误: {e}")

    print("任务完成！")
if __name__ == "__main__":
    process_single_file_by_bin()
