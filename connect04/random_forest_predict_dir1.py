# This script consolidates prediction logic to calculate the A_union_C model result.
# It loads two distinct models (A and C), runs predictions using their respective
# feature sets and filtering rules, calculates the logical union of their positive
# predictions, and outputs only the final positive (prediction == 1) results.
# D and E features have been removed from the processing pipeline as requested.

import pandas as pd
import numpy as np
import pickle
import os
import sys
import argparse
import time
from multiprocessing import Pool, cpu_count

# --- D/E features related global variables removed ---

def calculate_features_chunk(chunk_data):
    """
    为一组contig对计算所有潜在特征 (仅A, C)。
    这个函数不执行预测，只进行特征计算。

    Args:
        chunk_data (list of dicts): 带有方向性的contig对及其基础特征的列表。

    Returns:
        pandas.DataFrame: 包含所有计算出的特征和原始键的DataFrame。
    """

    rows_list = []

    # 为块中的每一对计算所有特征
    for data_point in chunk_data:
        c1_name = data_point['contig1']
        c2_name = data_point['contig2']
        c1_dir = data_point['dir1']
        c2_dir = data_point['dir2']

        row = {
            # 基础信息
            'contig1': c1_name,
            'contig2': c2_name,
            'contig1_direction': c1_dir,
            'contig2_direction': c2_dir,
            'binID': data_point.get('binID', np.nan),
            'label': data_point.get('label', np.nan), # 原始标签 (供参考，最终输出中可能移除)

            # A, C 特征 (直接从原始数据加载)
            'A_weight1': data_point.get('A_weight1', np.nan),
            'A_weight2': data_point.get('A_weight2', np.nan),
            'C1_cosine': data_point.get('C1_cosine', np.nan),
            'C2_cosine': data_point.get('C2_cosine', np.nan),
            
            # D and E features (Edge_Weight_D, E1_vector_diff, E2_cosine_sim) have been removed.
        }
        
        rows_list.append(row)

    return pd.DataFrame(rows_list)


def load_model_and_features(model_file, features_file):
    """加载 Random Forest 模型和对应的特征列表。"""
    try:
        with open(model_file, 'rb') as f:
            loaded_model = pickle.load(f)
        with open(features_file, 'rb') as f:
            loaded_features = pickle.load(f)
        return loaded_model, loaded_features
    except FileNotFoundError as e:
        print(f"Error: Model or feature file not found: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"Error loading model or features: {e}")
        sys.exit(1)


def load_raw_data(input_dir):
    """加载原始特征数据。D和E特征相关文件加载已被移除。"""
    print("--- 1. Loading Input Data ---")
    
    # 文件路径定义
    raw_feature_file = os.path.join(input_dir, 'combine_connect3_dir1.pkl')

    # 1. 加载原始特征数据 (所有对)
    try:
        with open(raw_feature_file, 'rb') as f:
            raw_data = pickle.load(f)
        print(f"Loaded raw features from '{raw_feature_file}'. Total pairs: {len(raw_data)}")
    except Exception as e:
        print(f"Error loading raw features file '{raw_feature_file}': {e}")
        sys.exit(1)

    # 2. 扁平化数据点
    directional_data_points = []
    for canonical_directional_key, features_dict in raw_data.items():
        c1_name, c1_dir, c2_name, c2_dir = canonical_directional_key

        data_point = {
            'contig1': c1_name,
            'contig2': c2_name,
            'dir1': c1_dir,
            'dir2': c2_dir,
            **features_dict, # 展开特征字典 (包含 A_weight1/2, C1/C2_cosine, binID, label等)
        }
        directional_data_points.append(data_point)

    # 3. D and E features loading and merging logic removed.
    print("D and E features loading skipped as requested.")

    # 返回扁平化的数据点列表
    return directional_data_points


def apply_prediction_logic(df_all_features, loaded_model, loaded_features, model_name):
    """
    应用模型特定的预处理（过滤和填充）和预测。
    """
    df_filtered = df_all_features.copy()
    
    # 1. 过滤: 如果模型需要 A-features，则移除缺少它们的对。
    a_features_to_check = [f for f in ['A_weight1', 'A_weight2'] if f in loaded_features]

    if a_features_to_check:
        df_filtered = df_filtered.dropna(subset=a_features_to_check).copy()
    
    if df_filtered.empty:
        return None

    # 2. 填充: 对所有特征的缺失值进行 0 填充。
    all_feature_cols = [
        'A_weight1', 'A_weight2', 'C1_cosine', 'C2_cosine',
    ]
    imputable_cols = [col for col in all_feature_cols if col in df_filtered.columns]
    df_filtered[imputable_cols] = df_filtered[imputable_cols].fillna(0)

    # 3. 预测
    X_predict = df_filtered[loaded_features]
    
    predictions = loaded_model.predict(X_predict)
    probabilities = loaded_model.predict_proba(X_predict)[:, 1] # 获取正类（1）的概率

    df_filtered['prediction'] = predictions
    df_filtered['probability'] = probabilities
    df_filtered['model_name'] = model_name
    
    # 返回关键列
    return df_filtered[['contig1', 'contig2', 'contig1_direction', 'contig2_direction', 
                        'binID', 'label', 'prediction', 'probability']].copy()


def predict_and_output_edges(args, input_dir, output_file, model_A_file, features_A_file, model_C_file, features_C_file):
    """
    核心预测逻辑函数：加载模型，计算特征，执行预测，计算 A_union_C 结果，并保存输出。
    直接将结果保存到指定的 output_file 路径。
    """
    final_output_file = output_file
    input_directory = input_dir

    print("--- Starting A_union_C Prediction (Excluding D & E Features) ---")
    start_time_total = time.time()

    # --- Step 1: 加载模型和特征 ---
    loaded_model_A, loaded_features_A = load_model_and_features(model_A_file, features_A_file)
    loaded_model_C, loaded_features_C = load_model_and_features(model_C_file, features_C_file)
    print("Successfully loaded both Model A and Model C.")

    # --- Step 2: 加载和准备原始数据 ---
    all_data_points = load_raw_data(input_directory)
    
    if not all_data_points:
        print("No contig pairs found to predict. Exiting.")
        return

    # --- Step 3: 并行计算所有特征 (仅 A 和 C) ---
    num_processes = cpu_count()
    print(f"\n--- 2. Calculating A and C Features in Parallel (Using {num_processes} processes) ---")

    # 将所有数据点分成多个块
    chunk_size = len(all_data_points) // num_processes + 1
    chunks = [all_data_points[i:i + chunk_size] for i in range(0, len(all_data_points), chunk_size)]

    # 使用 Pool 来并行执行特征计算
    with Pool(num_processes) as pool:
        results = pool.map(calculate_features_chunk, chunks)

    # 合并所有进程的结果 (包含 A 和 C 特征的完整 DataFrame)
    df_all_features = pd.concat(results, ignore_index=True)
    print(f"Total pairs with calculated features: {len(df_all_features)}")

    # --- Step 4: 应用模型 A 和模型 C 的预测逻辑 ---
    print("\n--- 3. Applying Model A Prediction Logic ---")
    df_pred_A = apply_prediction_logic(df_all_features, loaded_model_A, loaded_features_A, 'A1A2')
    
    print("\n--- 4. Applying Model C Prediction Logic ---")
    df_pred_C = apply_prediction_logic(df_all_features, loaded_model_C, loaded_features_C, 'C1C2')
    
    if df_pred_A is None or df_pred_C is None:
        print("Error: One or both models failed to produce predictions after filtering. Cannot compute union. Exiting.")
        return

    # --- Step 5: 计算 A_union_C 并集结果 ---
    print("\n--- 5. Calculating A_union_C Result ---")
    
    key_cols = ['contig1', 'contig2', 'contig1_direction', 'contig2_direction']
    
    # 提取关键信息并重命名 prediction 列
    df_A_subset = df_pred_A[key_cols + ['prediction']].rename(columns={'prediction': 'pred_A'})
    df_C_subset = df_pred_C[key_cols + ['prediction']].rename(columns={'prediction': 'pred_C'})
    
    # 合并两个结果 (以内连接保证只考虑两个模型都处理过的对)
    df_merged = pd.merge(df_A_subset, df_C_subset, on=key_cols, how='inner')
    
    # 提取原始的 binID 和 label 列 (从 A 的结果中，因为它们应该是相同的)
    df_metadata = df_pred_A[key_cols + ['binID', 'label']].drop_duplicates()
    df_merged = pd.merge(df_merged, df_metadata, on=key_cols, how='left')

    # 计算并集预测 (A_union_C): 当且仅当 A 或 C 模型预测为 1 时，并集为 1
    df_merged['prediction_union'] = ((df_merged['pred_A'] == 1) | (df_merged['pred_C'] == 1)).astype(int)
    
    # 筛选出预测结果为 1 的行
    df_positive_union = df_merged[df_merged['prediction_union'] == 1].copy()
    
    print(f"Total positive predictions found for A_union_C model: {len(df_positive_union)}")

    # --- Step 6: 保存最终输出 ---
    
    if not df_positive_union.empty:
        df_output = pd.DataFrame()
        df_output['contig1'] = df_positive_union['contig1']
        df_output['contig1_direction'] = df_positive_union['contig1_direction']
        df_output['contig2'] = df_positive_union['contig2']
        df_output['contig2_direction'] = df_positive_union['contig2_direction']
        
        # 并集模型没有直接的概率，统一设置为 1.0
        df_output['prediction_probability'] = 1.0 
        df_output['binID'] = df_positive_union['binID']
        
        # 记录是哪个模型贡献了正向预测
        is_A = (df_positive_union['pred_A'] == 1)
        is_C = (df_positive_union['pred_C'] == 1)

        conditions = [
            (is_A & is_C),    # A+C: A and C both predicted 1
            (is_A & ~is_C),   # A: Only A predicted 1
            (~is_A & is_C)    # C: Only C predicted 1
        ]
        
        choices = [
            'A+C',
            'A',
            'C'
        ]
        
        df_output['Source_Model'] = np.select(conditions, choices, default='ERROR')
        
        # 确保输出文件所在的目录存在
        output_directory = os.path.dirname(final_output_file)
        if output_directory and not os.path.exists(output_directory):
            os.makedirs(output_directory, exist_ok=True) # 使用 exist_ok=True 避免目录已存在时的错误

        # 直接保存到指定的 output_file 路径
        df_output.to_csv(
            final_output_file,
            sep='\t',
            index=False,
            header=True
        )
        print(f"\nSuccessfully saved {len(df_output)} positive predictions to '{final_output_file}'.")
    else:
        print("No edges were predicted as '1' by the A_union_C model. Output file was not created.")

    end_time_total = time.time()
    print(f"\n--- A_union_C Prediction finished successfully in {end_time_total - start_time_total:.2f} seconds. ---")


def main():
    """
    主函数: 负责解析命令行参数并调用核心预测逻辑函数。
    """
    parser = argparse.ArgumentParser(description='Predict contig pairs using two specified pre-trained model sets (A and C) and take the union of positive predictions.')
    
    # 模型 A 参数
    parser.add_argument('--model_A_file', type=str, required=True,
                        help='Path to the first model file (A1).')
    parser.add_argument('--features_A_file', type=str, required=True,
                        help='Path to the feature columns list file (A2) used by the first model.')

    # 模型 C 参数
    parser.add_argument('--model_C_file', type=str, required=True,
                        help='Path to the second model file (C1).')
    parser.add_argument('--features_C_file', type=str, required=True,
                        help='Path to the feature columns list file (C2) used by the second model.')

    # 数据输入和输出参数
    parser.add_argument('--input_dir', type=str, required=True,
                        help='Path to the base directory containing core input data files (e.g., combine_connect3_dir1.pkl).')
    # 更改为 --output_file 参数，指定完整的输出文件路径
    parser.add_argument('--output_file', type=str, required=True,
                        help='The full path to the prediction output file (e.g., /path/to/predictions_A_union_C.tsv).')
    
    # 辅助路径参数 (保留以兼容原始接口，尽管 D features 已移除)
    parser.add_argument('--output_path', type=str, required=True,
                        help='Base output path used to locate the cluster_res/leiden_initial_edge subdirectory for D features (e.g., the root of the main pipeline output).')

    args = parser.parse_args()
    
    # 调用核心预测主函数
    predict_and_output_edges(
        args=args, # 传递 args 访问 output_path
        input_dir=args.input_dir,
        output_file=args.output_file, # 传递 output_file
        model_A_file=args.model_A_file,
        features_A_file=args.features_A_file,
        model_C_file=args.model_C_file,
        features_C_file=args.features_C_file
    )


if __name__ == '__main__':
    main()