from .getid import *
from .extract_initial_edge import *
from . import best_feature_fre 
from .random_forest_predict_bd0 import *
from .multi_leiden_test import *
from scripts.clustering_config import get_final_clustering_grid
from . import score_cluster
from . import transfer_fa
import argparse
import logging
import os
import shlex
import subprocess


def resolve_checkm2_command(args):
    """Resolve the CheckM2 executable from config or PATH."""
    checkm2_bin = getattr(args, "checkm2_bin", None)
    if checkm2_bin:
        return checkm2_bin

    checkm2_path = getattr(args, "checkm2_path", None)
    if checkm2_path:
        candidate = os.path.join(checkm2_path, "bin", "checkm2")
        if os.path.exists(candidate):
            return candidate
        return checkm2_path

    return "checkm2"

def binning_init(args,logger):
    output_path = args.output_path
    
    
    estimate_file_path = args.output_path+'/cluster_res/estimate_res.txt'
    
    best_id_for_threshold,other_five_ids_string=get_best_and_other_ids_for_shell(estimate_file_path)
    print(best_id_for_threshold)
    print(other_five_ids_string)
    
    bandwith,maxedges,parameters,ratio = best_id_for_threshold.split('_')
    args.bandwidth_for_edge_extraction = float(bandwith)
    args.partgraph_ratio = float(ratio)

    extract_edges(logger, args)
    
    best_feature_fre.main(args)
    
    skip_codon_features = args.skip_codon_features
    project_dir = args.output_path+'/binning/s_cluster'
    output_enhanced_dir = project_dir+'/leiden_enhanced_edge'
    model_override_dir = args.linking_path +'/save_models'
    enhance_edges_with_rf_prediction(args,project_dir, output_enhanced_dir, skip_codon_features,model_override_dir)
    
    
    final_grid = get_final_clustering_grid(
        getattr(args, "clustering_mode", "full")
    )
    MAX_EDGES_LIST = final_grid["max_edges"]
    resolution_parameter_list = final_grid["resolution_parameters"]
    logger.info(
        "Final Leiden grid (%s mode): %d jobs",
        getattr(args, "clustering_mode", "full"),
        len(MAX_EDGES_LIST) * len(resolution_parameter_list),
    )
    input_dir = output_enhanced_dir
    seed_file = args.seed_file
    contig_file = args.contig_file
    output_dir = project_dir
    num_threads = args.num_threads
    partgraph_ratio = ratio
    bandwith  = bandwith
    

    for max_e in MAX_EDGES_LIST:
        run_all_clusterings(
            logger,
            input_dir,
            seed_file,
            contig_file,
            output_dir,
            num_threads,
            partgraph_ratio,
            bandwith,
            max_e,
            resolution_parameter_list=resolution_parameter_list,
        )
    score_cluster.main(args)
    args.output_path = output_path
    transfer_fa.main(args,logger)
    
    checkm2_cmd = resolve_checkm2_command(args)
    command = [
        checkm2_cmd,
        "predict",
        "--input",
        f"{args.output_path}/binning/bins",
        "--output-directory",
        f"{args.output_path}/binning/checkm2_bins",
        "-x",
        "fa",
        "-t",
        str(args.num_threads),
    ]
    logger.info("Running CheckM2: %s", " ".join(shlex.quote(part) for part in command))
    subprocess.run(command, check=True)

