import os
import sys
from .train_CLmodel import train_CLmodel
from .cluster import cluster
from .utils import gen_seed
from .get_final_result import run_get_final_result
from .generate_augfasta_and_saveindex import run_gen_augfasta
from .gen_cov import run_gen_cov
from .gen_var import run_gen_cov_var


def contrastive_learning_main(args, logger):

    # --- Step 1: Data Augmentation ---
    logger.info('\n--- Step 1/4: Generating Augmentation Data ---')
    data_aug_path = os.path.join(args.output_path, 'aug') 

    if not os.path.exists(data_aug_path):
        logger.info(f"Augmented data path '{data_aug_path}' not found. Generating augmented data.")
        os.makedirs(args.output_path, exist_ok=True)
        # Generate augmented FASTA and auxiliary metadata
        run_gen_augfasta(logger, args)
        run_gen_cov(logger, args)
        run_gen_cov_var(logger, args)
        
    else:
        logger.info(f"Augmented data path '{data_aug_path}' already exists. Skipping data generation.")
    # Ensure output directory exists

    # Generate augmented FASTA and auxiliary metadata
    # --- Step 2: Contrastive Learning Model Training ---
    logger.info('\n--- Step 2/4: Training Contrastive Learning Model ---')
    args.emb_file = os.path.join(args.output_path, 'embeddings.tsv')
    if not os.path.exists(args.emb_file):
        logger.info(f"embedding path '{args.emb_file}' not found. Generating embedding.")
        # train_CLmodel expects data generated in args.data
        train_CLmodel(logger, args)

        # Expected output embedding file (path depends on the training script)
        logger.info(f"Assuming embeddings will be saved to: {args.emb_file}")
    else:
        logger.info(f"embedding path '{args.emb_file}' already exists. Skipping embedding generation.")

    # --- Step 3: Binning / Clustering ---
    logger.info('\n--- Step 3/4: Clustering (Binning) ---')

    num_threads = args.num_threads

    # gen_seed prepares marker gene info and returns contig count or file path
    seed_num = gen_seed(
        args,
        logger,
        args.contig_file,
        num_threads,
        args.contig_len,
        marker_name="bacar_marker",
        quarter="2quarter"
    )

    # Perform clustering
    cluster(logger, args)

    # --- Step 4: Generate Final Results ---
    logger.info('\n--- Step 4/4: Generating Final Binning Results ---')

    # Regenerate seed info to ensure consistency
    seed_num = gen_seed(
        args,
        logger,
        args.contig_file,
        num_threads,
        args.contig_len,
        marker_name="bacar_marker",
        quarter="2quarter"
    )

    # Ensure clustering results path is set
    logger.info(f"Using binning results from: {args.binning_res_path}")

    # Generate final binning output
    run_get_final_result(
        logger,
        args,
        seed_num,
        num_threads,
        ignore_kmeans_res=True
    )
    
    logger.info("\nContrastive_Learning finished all steps successfully.")
