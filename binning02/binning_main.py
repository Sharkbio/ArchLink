from .getid import *
from .extract_initial_edge import *
from . import best_feature_fre 
from .random_forest_predict_bd0 import *
from .multi_leiden_test import *
from . import score_cluster
from . import transfer_fa
import argparse
import logging

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
    
    
    MAX_EDGES_LIST = [60,70,80, 85, 100]
    input_dir = output_enhanced_dir
    seed_file = args.seed_file
    contig_file = args.contig_file
    output_dir = project_dir
    num_threads = args.num_threads
    partgraph_ratio = ratio
    bandwith  = bandwith
    

    for max_e in MAX_EDGES_LIST:
        run_all_clusterings(logger,input_dir,seed_file,contig_file,output_dir,num_threads,partgraph_ratio,bandwith, max_e)
    score_cluster.main(args)
    args.output_path = output_path
    transfer_fa.main(args,logger)
    
    # checkm
    os.system(f'{args.checkm2_path}/bin/checkm2 predict --input {args.output_path}/binning/bins --output-directory {args.output_path}/binning/checkm2_bins -x fa -t {args.num_threads}')

