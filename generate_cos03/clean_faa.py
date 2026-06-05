import os
import sys

def remove_asterisks_from_file(input_filepath: str, output_filepath: str):
    """
    Reads an input file, removes all asterisk (*) characters from each line,
    and writes the modified content to an output file.

    Args:
        input_filepath (str): The path to the input file.
        output_filepath (str): The path to the output file where
                               the modified content will be saved.
    """
    try:
        # Check if the input file exists
        if not os.path.exists(input_filepath):
            print(f"Error: Input file not found at '{input_filepath}'")
            return

        print(f"Reading from: {input_filepath}")
        print(f"Writing to: {output_filepath}")

        with open(input_filepath, 'r', encoding='utf-8') as infile:
            with open(output_filepath, 'w', encoding='utf-8') as outfile:
                for line_num, line in enumerate(infile, 1):
                    # Remove all occurrences of '*' from the line
                    modified_line = line.replace('*', '')
                    outfile.write(modified_line)
        print("File processing complete. Asterisks removed successfully!")

    except IOError as e:
        print(f"IO Error: Could not read/write file. Details: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == "__main__":
    # --- Configuration ---
    # Use sys.argv to get input and output file paths from command line arguments.
    # Expected usage: python your_script_name.py <input_file_path> <output_file_path>
    if len(sys.argv) != 3:
        print("Usage: python clean_file.py <input_file_path> <output_file_path>")
        print("Example: python clean_file.py /path/to/your_input.fasta /path/to/your_output_cleaned.fasta")
        sys.exit(1) # Exit with an error code

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    # Call the function to perform the operation
    remove_asterisks_from_file(input_file, output_file)


