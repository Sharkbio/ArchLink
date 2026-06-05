from . import make_fa2
from .generate_graph_dir import *
from . import generate_predict_data_dir1
from .random_forest_predict_dir1 import *
import multiprocessing
import os
from pathlib import Path
import subprocess
import shutil
import sys


def clean_connect_graph(BINNING_DIR):
    """Remove connect_graph3.txt from each bin directory under linking/bins_0.9."""

    linking_dir = Path(BINNING_DIR) / "linking" / "bins_0.9"

    if not linking_dir.exists():
        print(f"Directory not found: {linking_dir}")
        return

    # Iterate over first-level subdirectories (bins)
    for bin_dir in linking_dir.iterdir():
        if bin_dir.is_dir():
            BIN_ID = bin_dir.name
            print(f"Entering: {bin_dir} (Bin ID: {BIN_ID})")

            graph_file = bin_dir / "connect_graph3.txt"
            if graph_file.is_file():
                graph_file.unlink()
                print(f"Deleted: {graph_file}")


def matching(BINNING_DIR, SCRIPT_PATH):
    """Run the external 'matching' executable for each bin directory."""

    bins_dir = Path(BINNING_DIR) / "linking" / "bins_0.9"

    if not bins_dir.exists():
        print(f"Error: Directory not found: {bins_dir}", file=sys.stderr)
        sys.exit(1)

    # Traverse each bin directory
    for bin_dir in bins_dir.iterdir():
        if not bin_dir.is_dir():
            continue

        BIN_ID = bin_dir.name
        print(f"Entering: {bin_dir} (Bin ID: {BIN_ID})")

        graph_file = bin_dir / "connect_graph_dir.txt"
        r_file = bin_dir / "connect_dir.r"
        c_file = bin_dir / "connect_dir.c"

        if not graph_file.is_file():
            print(f"Warning: connect_graph_dir.txt not found in {bin_dir}, skipping.")
            print("----------------------------------------")
            continue

        print(f"Found graph file: {graph_file}")
        print("Running matching...")

        cmd = [
            str(Path(SCRIPT_PATH) / "matching"),
            "-g", str(graph_file),
            "-r", str(r_file),
            "-c", str(c_file),
            "--verbose", "1",
            "-i", "10",
            "-b", "1",
            "--ignore_copy"
        ]

        try:
            result = subprocess.run(
                cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
            )
            print(f"Matching executed successfully for {BIN_ID}.")
        except subprocess.CalledProcessError as e:
            print(f"Error: matching failed for {BIN_ID}.", file=sys.stderr)
            print(f"Return code: {e.returncode}", file=sys.stderr)
            print(f"Error output:\n{e.stderr}", file=sys.stderr)
            sys.exit(1)
        except FileNotFoundError:
            print(f"Error: matching executable not found: {cmd[0]}", file=sys.stderr)
            sys.exit(1)

        print("----------------------------------------")


def make_output_dir(BINNING_DIR_ROOT):
    """Create output directories and copy/produce FASTA files for binning and connection."""

    bins_09_dir = Path(BINNING_DIR_ROOT) / "linking" / "bins_0.9"

    if not bins_09_dir.exists():
        print(f"Error: bins_0.9 not found: {bins_09_dir}", file=sys.stderr)
        sys.exit(1)

    parent_dir = bins_09_dir.parent
    binning_target_dir = parent_dir / "binning"
    connect_target_dir = parent_dir / "connect"

    binning_target_dir.mkdir(exist_ok=True)
    connect_target_dir.mkdir(exist_ok=True)

    for bin_dir in bins_09_dir.iterdir():
        if not bin_dir.is_dir():
            continue

        BIN_ID = bin_dir.name
        print(f"Entering: {bin_dir} (Bin ID: {BIN_ID})")

        # 1. Copy original FASTA to binning
        original_fasta = bin_dir / f"{BIN_ID}.fasta"
        target_binning_fasta = binning_target_dir / f"{BIN_ID}_g.fa"

        if original_fasta.is_file():
            print(f"Found original FASTA: {original_fasta}")
            shutil.copyfile(original_fasta, target_binning_fasta)
            print(f"Copied to binning: {target_binning_fasta}")
        else:
            print(f"Warning: Original FASTA not found: {original_fasta}")

        # 2. Process connect FASTA using make_fa2 (or fallback to original)
        graph_file = bin_dir / "connect_graph_dir.txt"
        r_file = bin_dir / "connect_dir.r"
        out_fasta = bin_dir / "connect.fasta"
        target_connect_fasta = connect_target_dir / f"{BIN_ID}_c.fa"
        fallback_fasta = connect_target_dir / f"{BIN_ID}_g.fa"

        def copy_original_to_connect():
            """Fallback: copy original FASTA to connect directory."""
            if original_fasta.is_file():
                shutil.copyfile(original_fasta, fallback_fasta)
                print(f"Copied original FASTA to connect: {fallback_fasta}")
            else:
                print("Warning: Original FASTA missing, cannot copy.")

        if not graph_file.is_file():
            print("Warning: connect_graph_dir.txt not found, skipping make_fa.py.")
            copy_original_to_connect()
        else:
            print(f"Found: {graph_file}")
            print("Running make_fa.py...")

            if not (original_fasta.is_file() and r_file.is_file()):
                print("Warning: Missing FASTA or .r file, cannot generate connect FASTA.")
                copy_original_to_connect()
            else:
                make_fa2.main(original_fasta, r_file, out_fasta)
                try:
                    print(f"make_fa.py executed successfully for {BIN_ID}.")

                    if out_fasta.is_file():
                        shutil.copyfile(out_fasta, target_connect_fasta)
                        print(f"Copied connect FASTA: {target_connect_fasta}")
                    else:
                        print("Error: make_fa.py completed but no connect.fasta produced.")
                        copy_original_to_connect()

                except subprocess.CalledProcessError as e:
                    print(f"Error: make_fa.py failed for {BIN_ID}.", file=sys.stderr)
                    print(f"Error output:\n{e.stderr}", file=sys.stderr)
                    print("Falling back to original FASTA.")
                    copy_original_to_connect()

        print("----------------------------------------")


def connect_main(args):

    # 1. Generate prediction input data
    bins_0_9 = args.output_path + '/linking/bins_0.9'
    bam_graph_dir = args.output_path + '/bam.graph'
    raw_cos = args.output_path + '/linking/cosine/cosine_model_features_softmax.pkl'
    outs = args.output_path + '/linking/combine_connect3_dir1.pkl'
    generate_predict_data_dir1.main(bins_0_9, bam_graph_dir, raw_cos, outs)

    # 2. Predict edges using Random Forest model
    # feature_path = args.linking_path + '/save_models/feature_columns_gas_connect_COMB_A_weight1_A_weight2_C1_cosine_C2_cosine3.pkl'
    # model_path = args.linking_path + '/save_models/best_random_forest_model_gas_connect_COMB_A_weight1_A_weight2_C1_cosine_C2_cosine3.pkl'
    save_models_dir = os.path.join(args.linking_path, 'save_models')
    feature_path1 = os.path.join(save_models_dir, 'feature_columns_gas_connect_COMB_A_weight1_A_weight23.pkl')
    model_path1 = os.path.join(save_models_dir, 'best_random_forest_model_gas_connect_COMB_A_weight1_A_weight23.pkl')
    feature_path2 = os.path.join(save_models_dir, 'feature_columns_gas_connect_COMB_C1_cosine_C2_cosine3.pkl')
    model_path2 = os.path.join(save_models_dir, 'best_random_forest_model_gas_connect_COMB_C1_cosine_C2_cosine3.pkl')
    in_out_dir = args.output_path + '/linking/'
    output_basename = 'predictions_A1A2C1C2.tsv'
    output_file = os.path.join(in_out_dir, output_basename)

    use_dual_model = all(
        os.path.exists(path) for path in (model_path1, feature_path1, model_path2, feature_path2)
    )

    if use_dual_model:
        predict_and_output_edges(
            args,
            in_out_dir,
            output_file,
            model_path1,
            feature_path1,
            model_path2,
            feature_path2,
        )
    else:
        raise FileNotFoundError(
            "Required linking random-forest model artifacts were not found in save_models/."
        )
    # 3. Clean old graph files
    # clean_connect_graph(args.output_path)

    # 4. Split predictions into bin-level files
    process_single_file_by_bin(output_file,bins_0_9)
    # 5. Run graph matching
    matching(args.output_path, SCRIPT_PATH=args.linking_path + '/save_models/')

    # 6. Generate final output directories
    make_output_dir(args.output_path)


if __name__ == "__main__":
    connect_main()
