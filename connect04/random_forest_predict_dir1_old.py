import pandas as pd
import numpy as np
import pickle
import os
import sys
import argparse
import time
from sklearn.metrics.pairwise import cosine_similarity
from scipy.spatial.distance import euclidean
from multiprocessing import Pool, cpu_count
from functools import partial

# --- 全局变量，将在所有工作进程之间共享 ---
full_data_dict = None
contig_vectors_dict = None
loaded_best_model = None
loaded_feature_columns = None


def init_worker(data_dict, contig_vectors, model, features):
    """
    每个工作进程的初始化函数。
    将大型数据加载到每个工作进程的全局作用域中。
    """
    global full_data_dict, contig_vectors_dict, loaded_best_model, loaded_feature_columns
    full_data_dict = data_dict
    contig_vectors_dict = contig_vectors
    loaded_best_model = model
    loaded_feature_columns = features


def process_chunk(chunk_data):
    """
    Compute features for a set of contig pairs, apply filtering and imputation, and perform prediction.

    Args:
        chunk_data (list of dicts): A list of contig pairs with directionality and their features.

    Returns:
        pandas.DataFrame: DataFrame containing predicted contig pairs with predicted class (0 or 1) and probability.
    """
    # Access the global data dictionaries loaded by the initializer.
    global full_data_dict, contig_vectors_dict, loaded_best_model, loaded_feature_columns

    # Pre-allocate a list to store feature rows
    rows_list = []

    # Compute all features for each pair in the chunk
    for data_point in chunk_data:
        c1_name = data_point['contig1']
        c2_name = data_point['contig2']
        c1_dir = data_point['dir1']
        c2_dir = data_point['dir2']

        # Retrieve specific features used in model training (A_weight1/2, C1/C2_cosine)
        row = {
            'contig1': c1_name,
            'contig2': c2_name,
            'contig1_direction': c1_dir,
            'contig2_direction': c2_dir,
            'binID': data_point.get('binID', np.nan),

            # --- Match feature names used in the training script ---
            'A_weight1': data_point.get('A_weight1', np.nan),
            'A_weight2': data_point.get('A_weight2', np.nan),
            'C1_cosine': data_point.get('C1_cosine', np.nan),
            'C2_cosine': data_point.get('C2_cosine', np.nan),
            # --- Keep D feature, loaded and filled by main function ---
            'Edge_Weight_D': data_point.get('Edge_Weight_D', np.nan)
        }

        # --- E feature computation (mirrors training script logic) ---
        # Compute E features only if the model requires them
        is_e_required = any(f in ['E1_vector_diff', 'E2_cosine_sim'] for f in loaded_feature_columns)
        
        row['E1_vector_diff'] = np.nan
        row['E2_cosine_sim'] = np.nan

        if is_e_required and c1_name in contig_vectors_dict and c2_name in contig_vectors_dict:
            vec1 = contig_vectors_dict[c1_name]
            vec2 = contig_vectors_dict[c2_name]
            # E1: Euclidean distance
            row['E1_vector_diff'] = euclidean(vec1, vec2)
            # E2: Cosine similarity
            row['E2_cosine_sim'] = cosine_similarity(vec1.reshape(1, -1), vec2.reshape(1, -1))[0][0]

        rows_list.append(row)

    # Create a DataFrame for this chunk
    df_chunk = pd.DataFrame(rows_list)

    # --- Data processing (mirrors training script logic) ---
    # 1. Filter out pairs missing 'A' features. Training script used 'A_weight1' and 'A_weight2'
    # Filtering criterion: if model uses any A features, filter by them
    a_features_to_check = [f for f in ['A_weight1', 'A_weight2'] if f in loaded_feature_columns]

    # Skip filtering if the model doesn't use A features
    if a_features_to_check:
        present_a_features = [f for f in a_features_to_check if f in df_chunk.columns]
        if present_a_features:
            df_filtered = df_chunk.dropna(subset=present_a_features).copy()
        else:
            df_filtered = df_chunk.copy()
    else:
        df_filtered = df_chunk.copy()
    
    # Handle case where all rows are dropped after filtering
    if df_filtered.empty:
        # Create an empty DataFrame with correct columns to avoid concat issues
        return pd.DataFrame(columns=['contig1', 'contig2', 'contig1_direction', 'contig2_direction', 'binID', 'prediction', 'probability'])

    # 2. Fill remaining missing values in all features with 0
    all_feature_columns = [
        'A_weight1', 'A_weight2', 'C1_cosine', 'C2_cosine', 
        'Edge_Weight_D', 'E1_vector_diff', 'E2_cosine_sim'
    ]
    imputable_cols = [col for col in all_feature_columns if col in df_filtered.columns]
    
    # Fill any NaN values with 0
    df_filtered[imputable_cols] = df_filtered[imputable_cols].fillna(0)

    # Select only the feature columns used by the model
    X_predict = df_filtered[loaded_feature_columns]
    
    # --- Fix warning: convert to NumPy array for prediction ---
    X_predict_array = X_predict.values

    # 3. Make predictions
    predictions = loaded_best_model.predict(X_predict_array)
    probabilities = loaded_best_model.predict_proba(X_predict_array)[:, 1] # probability for positive class (1)

    # Add predictions and probabilities back to DataFrame
    df_filtered['prediction'] = predictions
    df_filtered['probability'] = probabilities

    # Return only positive predictions
    return df_filtered[df_filtered['prediction'] == 1].copy()


def predict_and_output_edges(args, input_dir, output_file, model_file, features_file):
    """
    Main function: load data and model, perform parallel prediction, and output results.
    
    Args:
        input_dir (str): Directory containing raw feature data.
        output_file (str): Path to save the final prediction results.
        model_file (str): Full path to the pre-trained model file (.pkl).
        features_file (str): Full path to the model feature file (list of feature column names) (.pkl).
    """
    print("="*80)
    print(f"--- Starting prediction for model: {os.path.basename(model_file)} ---")
    print(f"--- Output will be saved to: {output_file} ---")
    start_time_total = time.time()

    # --- Step 1: Load all required files ---
    output_feature_file = os.path.join(input_dir, 'combine_connect3_dir1.pkl')
    contig_vectors_file = os.path.join(input_dir, 'codon.pkl')
    
    initial_edge_dir = args.output_path + '/cluster_res/leiden_initial_edge'
    namelist_path = os.path.join(initial_edge_dir, 'namelist.txt')
    extracted_edges_path = os.path.join(initial_edge_dir, 'extracted_edges.npz')

    # Load model and feature columns
    try:
        with open(model_file, 'rb') as f:
            loaded_model = pickle.load(f)
        print(f"Successfully loaded model from '{model_file}'.")

        with open(features_file, 'rb') as f:
            loaded_features = pickle.load(f)
        print(f"Successfully loaded feature list from '{features_file}'. Features: {loaded_features}")

    except FileNotFoundError as e:
        print(f"Error: Required model/feature file not found: {e}")
        return
    except Exception as e:
        print(f"Error loading model or features: {e}")
        return

    # Load raw data
    try:
        with open(output_feature_file, 'rb') as f:
            raw_data = pickle.load(f)
        print(f"Successfully loaded raw features from '{output_feature_file}'. Total pairs: {len(raw_data)}")

        directional_data_points = []
        for canonical_directional_key, features_dict in raw_data.items():
            c1_name, c1_dir, c2_name, c2_dir = canonical_directional_key
            
            data_point = {
                'contig1': c1_name,
                'contig2': c2_name,
                'dir1': c1_dir,
                'dir2': c2_dir,
                **features_dict,  # Expand feature dictionary (includes A_weight1/2, C1/C2_cosine, binID, etc.)
            }
            directional_data_points.append(data_point)

        # Load contig vectors if any E features are required
        contig_vectors = {}
        if any(f in ['E1_vector_diff', 'E2_cosine_sim'] for f in loaded_features):
            with open(contig_vectors_file, 'rb') as f:
                contig_vectors = pickle.load(f)
            print(f"Successfully loaded contig vectors from '{contig_vectors_file}'.")
        else:
            print("E features are not in the model, skipping contig vector loading.")

        # Load 'D' feature if required
        if 'Edge_Weight_D' in loaded_features:
            if not os.path.exists(extracted_edges_path) or not os.path.exists(namelist_path):
                print(f"Warning: D feature required, but data files ({extracted_edges_path} or {namelist_path}) not found. Skipping D feature loading.")
            else:
                loaded_data = np.load(extracted_edges_path)
                sources = loaded_data['sources']
                targets = loaded_data['targets']
                weights = loaded_data['weights']

                namelist_df = pd.read_csv(namelist_path, header=None)
                namelist = namelist_df[0].tolist()
                
                # Create a convenient dictionary for D feature (undirected edge weights)
                edge_weights_D = {}
                for i in range(len(sources)):
                    u_idx, v_idx = sources[i], targets[i]
                    u_name, v_name = namelist[u_idx], namelist[v_idx]
                    pair = tuple(sorted((u_name, v_name)))
                    edge_weights_D[pair] = weights[i]
                
                # Add D feature to each flattened data point (using canonical undirected key)
                for dp in directional_data_points:
                    canonical_pair = tuple(sorted((dp['contig1'], dp['contig2'])))
                    dp['Edge_Weight_D'] = edge_weights_D.get(canonical_pair, np.nan)
                
                print(f"Successfully loaded and merged D features.")
        else:
            print("D feature is not in the model, skipping extracted edge loading.")

    except Exception as e:
        print(f"Error loading data files: {e}")
        return

    # --- Step 2: Prepare data for prediction ---
    all_data_points = directional_data_points

    if not all_data_points:
        print("No contig pairs found to predict. Skipping.")
        return

    # --- Step 3: Parallel prediction ---
    num_processes = cpu_count()
    print(f"\nUsing {num_processes} processes for parallel prediction.")

    # Split all data points into chunks for each process
    chunk_size = len(all_data_points) // num_processes + 1
    chunks = [all_data_points[i:i + chunk_size] for i in range(0, len(all_data_points), chunk_size)]

    # Use Pool for parallel execution
    with Pool(num_processes, initializer=init_worker,
              initargs=({}, contig_vectors, loaded_model, loaded_features)) as pool:

        print(f"Spawning {len(chunks)} prediction workers...")
        results = pool.map(process_chunk, chunks)

    # Concatenate results from all processes
    df_all_predictions = pd.concat(results, ignore_index=True)

    # Sort predictions by probability for easier inspection
    df_all_predictions = df_all_predictions.sort_values(by='probability', ascending=False)

    print(f"\nTotal positive predictions found for this model: {len(df_all_predictions)}")

    # --- Step 4: Save results ---
    if not df_all_predictions.empty:
        df_output = pd.DataFrame()
        df_output['contig1'] = df_all_predictions['contig1']
        df_output['contig1_direction'] = df_all_predictions['contig1_direction']
        df_output['contig2'] = df_all_predictions['contig2']
        df_output['contig2_direction'] = df_all_predictions['contig2_direction']
        df_output['prediction_probability'] = df_all_predictions['probability']

        if 'binID' in df_all_predictions.columns:
            df_output['binID'] = df_all_predictions['binID']
        
        # Add used feature list to output
        df_output['Model_Features'] = "+".join(loaded_features)
        
        # Ensure output directory exists
        output_dir = os.path.dirname(output_file)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)

        df_output.to_csv(
            output_file,
            sep='\t',
            index=False,
            header=True
        )
        print(f"Successfully saved {len(df_output)} positive predictions to '{output_file}'.")
    else:
        print("No edges were predicted as '1'. Output file was not created.")

    end_time_total = time.time()
    print(f"--- Prediction finished successfully in {end_time_total - start_time_total:.2f} seconds. ---")
    print("="*80)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Predict contig pairs using a specified pre-trained model and feature list.')
    
    # 新增参数: 模型文件和特征文件路径 (必需)
    parser.add_argument('--model_file', type=str, required=True,
                        help='Path to the pre-trained model file (.pkl).')
    parser.add_argument('--features_file', type=str, required=True,
                        help='Path to the feature columns list file (.pkl) used by the model.')

    # 保持原有数据输入和输出目录参数 (必需)
    parser.add_argument('--input_dir', type=str, required=True,
                        help='Path to the base directory containing all necessary input data files (e.g., combine_connect3_dir2.pkl, codon.pkl, and leiden_initial_edge subdirectory).')
    parser.add_argument('--output_dir', type=str, default='.',
                        help='The directory where the prediction output file will be saved. Default is the current directory.')
    
    args = parser.parse_args()
    
    # 构造输出文件名: 使用模型文件名作为基础
    model_basename = os.path.basename(args.model_file).replace('.pkl', '')
    output_basename = f'predictions_A1A2C1C2.tsv'
    output_file = os.path.join(args.output_dir, output_basename)

    print(f"--- Starting Single Prediction Task ---")
    print(f"  Model File: {args.model_file}")
    print(f"  Features File: {args.features_file}")
    print(f"  Input Data Directory: {args.input_dir}")
    print(f"  Final Output File: {output_file}")


    # 调用预测主函数
    predict_and_output_edges(
        args.input_dir, 
        output_file, 
        args.model_file, 
        args.features_file
    )

    print("\n--- The prediction task has been completed. ---")