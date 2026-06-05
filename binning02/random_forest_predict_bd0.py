import pandas as pd
import numpy as np
import pickle
import os
import sys
from collections import defaultdict
import argparse
from sklearn.metrics.pairwise import cosine_similarity
from scipy.spatial.distance import euclidean
from multiprocessing import Pool, cpu_count, Manager # 引入 Manager

# 全局变量将保存代理对象（针对大数据）或实际对象（针对小数据/模型）
# 在工作进程的全局范围内。
raw_data = None
contig_vectors = None
edge_weights_D_lookup = None
loaded_best_model = None
loaded_feature_columns = None
calculate_e_features = True

def init_worker(raw_data_proxy, contig_vectors_proxy, edge_weights_D_lookup_proxy, model, features, load_e_features):
    """
    每个工作进程的初始化函数。
    此函数将轻量级代理对象、模型和特征加载到每个工作进程的全局范围。
    """
    global raw_data, contig_vectors, edge_weights_D_lookup, loaded_best_model, loaded_feature_columns, calculate_e_features
    
    # 将 Manager 代理对象分配给工作进程中的全局变量
    raw_data = raw_data_proxy
    contig_vectors = contig_vectors_proxy
    edge_weights_D_lookup = edge_weights_D_lookup_proxy
    
    # 模型和特征列很小，可以安全地通过值传递/复制
    loaded_best_model = model
    loaded_feature_columns = features
    calculate_e_features = load_e_features

def process_chunk(chunk_data):
    """
    计算一段 contig 对的特征并进行预测。
    它通过 Manager 代理（全局变量）访问大型数据结构。
    """
    global raw_data, contig_vectors, edge_weights_D_lookup, loaded_best_model, loaded_feature_columns, calculate_e_features
    
    # 根据过滤后的 contig 对计算所有特征
    b_frequency = [raw_data.get(pair, {}).get('B_frequency', 0.0) for pair in chunk_data]
    c_cosine_max = [raw_data.get(pair, {}).get('C_cosine_max', np.nan) for pair in chunk_data]
    edge_weight_d = [edge_weights_D_lookup.get(pair, np.nan) for pair in chunk_data]

    # 确定是否计算 E 特征
    e1_vector_diff = []
    e2_cosine_sim = []
    
    if calculate_e_features:
        for pair in chunk_data:
            c1_name, c2_name = pair
            # 通过代理访问 contig 向量
            if c1_name in contig_vectors and c2_name in contig_vectors:
                vector_c1 = contig_vectors[c1_name]
                vector_c2 = contig_vectors[c2_name]
                e1_vector_diff.append(euclidean(vector_c1, vector_c2))
                # 调整形状以适应 cosine_similarity，要求 2D 数组：(1, n_features)
                e2_cosine_sim.append(cosine_similarity(vector_c1.reshape(1, -1), vector_c2.reshape(1, -1))[0][0])
            else:
                e1_vector_diff.append(np.nan)
                e2_cosine_sim.append(np.nan)
    else:
        # 如果不计算 E 特征，则用 NaNs 填充
        e1_vector_diff = [np.nan] * len(chunk_data)
        e2_cosine_sim = [np.nan] * len(chunk_data)

    # 为此块创建 DataFrame。
    df_chunk = pd.DataFrame({
        'contig_pair': chunk_data,
        'B_frequency': b_frequency,
        'C_cosine_max': c_cosine_max,
        'Edge_Weight_D': edge_weight_d,
        'E1_vector_diff': e1_vector_diff,
        'E2_cosine_sim': e2_cosine_sim
    })
    
    # 过滤并准备数据进行预测。
    df_for_prediction = df_chunk.dropna(subset=loaded_feature_columns).copy()

    if df_for_prediction.empty:
        # 返回具有预期列的空结果
        return pd.DataFrame(columns=['contig_pair', 'prediction_proba'])

    X_for_prediction = df_for_prediction[loaded_feature_columns]

    # 确保模型已加载并准备就绪
    if loaded_best_model is None:
        raise RuntimeError("模型未在工作进程中初始化。")

    # 预测类别 1（正确边）的概率
    y_pred_proba_correct = loaded_best_model.predict_proba(X_for_prediction)[:, 1]
    df_for_prediction['prediction_proba'] = y_pred_proba_correct

    return df_for_prediction[['contig_pair', 'prediction_proba']]

# 修正后的函数签名，以匹配外部调用者传递的 5 个参数：(args, project_dir, output_enhanced_dir, skip_codon_features, model_override_dir)
def enhance_edges_with_rf_prediction(args_from_caller, project_dir, output_enhanced_dir, skip_codon_features, model_override_dir):
    """
    加载预训练的随机森林模型，用它来预测“正确”的 contig 对，
    增强这些预测边的权重，并保存修改后的边数据。
    
    此函数已修改为使用 multiprocessing.Manager 共享大型数据结构，解决潜在的内存溢出问题。
    """
    
    # 核心修复: 从 args_from_caller 对象中提取 num_threads
    try:
        num_threads = args_from_caller.num_threads
    except AttributeError:
        # 如果调用者传递的 args 对象没有 num_threads 属性，则使用默认值
        print("警告: 在传入的 'args' 对象中未找到 'num_threads'。回退到 CPU 核心数。")
        num_threads = cpu_count() 
    
    model_suffix = "_D_B"
    
    # --- 0. 从参数中派生路径变量 ---
    # project_dir 现在是正确的字符串
    effective_model_dir = model_override_dir if model_override_dir else project_dir
    input_dir = os.path.join(project_dir, 'leiden_initial_edge') 

    print(f"--- 启动项目根目录 {project_dir} 的 RF 预测边增强 ---")
    print(f"原始边的输入目录: {input_dir}")
    print(f"模型和特征文件的目录（有效）: {effective_model_dir}")
    print(f"增强后的边的输出目录: {output_enhanced_dir}")
    print(f"使用的模型后缀: {model_suffix}")
    print(f"使用 {num_threads} 个工作进程。") # 打印确认线程数

    # --- 1. 定义特定文件路径 ---
    loaded_model_filename = os.path.join(
        effective_model_dir, f'best_random_forest_model_focus0{model_suffix}2.pkl'
    )
    loaded_feature_columns_filename = os.path.join(
        effective_model_dir, f'feature_columns_focus0{model_suffix}2.pkl'
    )
    raw_feature_data_file = os.path.join(project_dir, 'combine_binning.pkl')
    contig_vectors_file = os.path.join(project_dir, 'codon.pkl')

    os.makedirs(output_enhanced_dir, exist_ok=True)

    # --- 2. 检查并加载模型和特征列名称 ---
    raw_data_temp = {}
    contig_vectors_temp = {}
    
    try:
        # 加载模型和特征列（小对象 - 安全地正常加载）
        with open(loaded_model_filename, 'rb') as f:
            loaded_best_model = pickle.load(f)
        print(f"成功加载模型: {loaded_model_filename}")

        with open(loaded_feature_columns_filename, 'rb') as f:
            loaded_feature_columns = pickle.load(f)
        print(f"成功加载特征列名称: {loaded_feature_columns_filename}")

        # 加载原始特征数据 (大型对象 1)
        with open(raw_feature_data_file, 'rb') as f:
            raw_data_temp = pickle.load(f)
        print(f"成功加载原始特征数据文件: {raw_feature_data_file}")

        # 加载 contig 向量 (大型对象 2)
        if not skip_codon_features:
            with open(contig_vectors_file, 'rb') as f:
                contig_vectors_temp = pickle.load(f)
            print(f"成功加载所有 contig 向量: {contig_vectors_file}")
        else:
            print("跳过 codon.pkl 加载。")

    except Exception as e:
        print(f"加载模型、特征列名称或数据文件时出错: {e}。退出。")
        return

    # --- 3. 加载原始边数据和 contig 信息 ---
    
    namelist_path = os.path.join(input_dir, 'namelist.txt')
    extracted_edges_path = os.path.join(input_dir, 'extracted_edges.npz')
    length_weight_path = os.path.join(input_dir, 'length_weight.txt')

    if not os.path.exists(namelist_path) or \
       not os.path.exists(extracted_edges_path) or \
       not os.path.exists(length_weight_path):
        print(f"错误: 在派生的输入目录 '{input_dir}' 中找不到所需的边文件。退出。")
        return

    try:
        namelist_df = pd.read_csv(namelist_path, header=None)
        original_namelist = namelist_df[0].tolist()
        name_to_idx = {name: i for i, name in enumerate(original_namelist)} 
        print(f"成功加载包含 {len(original_namelist)} 个 contig 的原始命名列表。")

        loaded_data = np.load(extracted_edges_path)
        original_sources = loaded_data['sources']
        original_targets = loaded_data['targets']
        original_weights = loaded_data['weights']
        print(f"成功加载 {len(original_sources)} 条原始边（特征 'D'）。")

        length_weight_df = pd.read_csv(length_weight_path, header=None)
        original_length_weight = length_weight_df[0].tolist()
        print(f"成功加载原始长度权重。")

    except Exception as e:
        print(f"从 '{input_dir}' 加载原始边数据或 contig 信息时出错: {e}。退出。")
        return

    # 创建 D 特征查找表 (大型对象 3)
    edge_weights_D_lookup_temp = {}
    for i in range(len(original_sources)):
        u_idx, v_idx = original_sources[i], original_targets[i]
        u_name, v_name = original_namelist[u_idx], original_namelist[v_idx]
        edge_weights_D_lookup_temp[tuple(sorted((u_name, v_name)))] = original_weights[i]

    print(f"为 {len(edge_weights_D_lookup_temp)} 对创建了 D 特征查找。")

    # --- 4. 使用 Manager 并行化特征计算和预测 ---
    
    all_d_pairs = list(edge_weights_D_lookup_temp.keys())
    print(f"总 D 特征对数: {len(all_d_pairs)}")
    
    filtered_contig_pairs = all_d_pairs
    print(f"使用所有 D 特征对作为基础: {len(filtered_contig_pairs)} 对")

    # 使用 multiprocessing.Manager 的上下文管理器
    with Manager() as manager:
        print("\n--- 将大型数据结构转换为 Manager 代理以共享内存 ---")
        
        # 将大型字典转换为由 Manager 管理的共享字典
        shared_raw_data = manager.dict(raw_data_temp)
        shared_edge_weights_D_lookup = manager.dict(edge_weights_D_lookup_temp)
        shared_contig_vectors = manager.dict(contig_vectors_temp) 

        # 释放主进程内存中的原始大型对象
        del raw_data_temp
        del edge_weights_D_lookup_temp
        del contig_vectors_temp
        print("已从主进程内存中释放原始大型数据对象。")

        # 对预过滤的数据进行分块
        chunk_size = len(filtered_contig_pairs) // num_threads + 1
        chunks = [filtered_contig_pairs[i:i + chunk_size] for i in range(0, len(filtered_contig_pairs), chunk_size)]

        print(f"\n--- 使用 {num_threads} 个进程并行化特征计算和预测 ---")
        print(f"总共要并行处理的对数: {len(filtered_contig_pairs)}")
        print(f"块的数量: {len(chunks)}")
        
        # 使用带有共享代理对象的进程池
        try:
            with Pool(num_threads, 
                      initializer=init_worker, 
                      initargs=(shared_raw_data, 
                                shared_contig_vectors, 
                                shared_edge_weights_D_lookup, 
                                loaded_best_model, 
                                loaded_feature_columns, 
                                not skip_codon_features)) as pool:
                
                # map 函数只需要传递较小的块列表
                chunk_results = pool.map(process_chunk, chunks)
            
            # 将所有进程的结果连接成一个 DataFrame。
            df_all_features = pd.concat(chunk_results, ignore_index=True)

            print("\n--- 所有特征数据预览（前 5 行）---")
            print(df_all_features.head())

        except Exception as e:
            print(f"多进程期间出错: {e}。退出。")
            return

    # --- 5. 组合和过滤数据进行增强（在 Manager 上下文之外） ---
    df_for_enhancement = pd.DataFrame(df_all_features, columns=['contig_pair', 'prediction_proba'])
    df_for_enhancement = df_for_enhancement.dropna(subset=['prediction_proba']).set_index('contig_pair')
    
    if df_for_enhancement.empty:
        print("\n错误: 所有数据点因缺少特征值而被删除。无法执行预测。")
        return

    print(f"\n过滤后用于预测的样本数: {len(df_for_enhancement)}。")

    # --- 6. 根据预测增强原始边权重 ---
    print("\n--- 启动边增强 ---")
    prediction_lookup = df_for_enhancement['prediction_proba'].to_dict()
    
    enhanced_sources = []
    enhanced_targets = []
    enhanced_weights = []

    prediction_threshold = 0.5
    
    # 遍历原始边并应用权重调整
    for i in range(len(original_sources)):
        u_idx, v_idx = original_sources[i], original_targets[i]
        u_name, v_name = original_namelist[u_idx], original_namelist[v_idx]
        original_weight = original_weights[i]

        pair = tuple(sorted((u_name, v_name)))
        
        prediction_prob = prediction_lookup.get(pair)
        
        # 核心修改: 根据预测概率调整权重
        if prediction_prob is not None and prediction_prob < prediction_threshold:
            # 当预测为 0（不正确）时，将权重降低到十分之一
            new_weight = original_weight / 10.0
        else:
            # 当预测为 1（正确）或没有找到预测时，保持原始权重
            new_weight = original_weight
        
        enhanced_sources.append(u_idx)
        enhanced_targets.append(v_idx)
        enhanced_weights.append(new_weight)
    
    enhanced_weights = np.array(enhanced_weights, dtype=np.float32)

    print(f"\n增强了 {len(enhanced_weights)} 条边。正在保存增强后的数据。")

    # --- 7. 将增强后的数据保存到输出目录 ---
    np.savez(os.path.join(output_enhanced_dir, 'extracted_edges.npz'),
             sources=np.array(enhanced_sources, dtype=np.int32),
             targets=np.array(enhanced_targets, dtype=np.int32),
             weights=enhanced_weights)

    pd.DataFrame(original_namelist).to_csv(os.path.join(output_enhanced_dir, 'namelist.txt'), index=False, header=False)
    pd.DataFrame(original_length_weight).to_csv(os.path.join(output_enhanced_dir, 'length_weight.txt'), index=False, header=False)

    print(f"\n增强后的边数据已保存到 '{output_enhanced_dir}'。")
    print("此目录现在可以用作 'run_clustering.py' 的 --input_dir。")
    print("--- 边增强完成 ---")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Enhance edge weights based on Random Forest prediction using a unified project directory structure.')
    
    parser.add_argument('--project_dir', type=str, required=True,
                        help='The root directory of the project (e.g., "my_project_data/binning").')
    
    parser.add_argument('--model_dir', type=str, required=False, default=None,
                        help='Optional: Specify a directory if model files are located outside the project_dir.')

    parser.add_argument('--output_enhanced_dir', type=str, required=True, # Made required for clarity
                        help='Path to the new directory where enhanced edge data will be saved (e.g., "my_project_data/binning/enhanced_with_rf").')
    
    # ADDED: Argument for controlling the number of threads/processes
    parser.add_argument('--num_threads', type=int, default=cpu_count(),
                        help=f'Number of processes to use for parallel feature calculation. Default is CPU count ({cpu_count()}). Reducing this helps with OOM.')
    
    parser.add_argument('--skip_codon_features', action='store_true',
                        help='If this flag is present, codon.pkl will not be loaded and E-features will be ignored.')
    
    args = parser.parse_args()
    
    # 使用与外部代码相同的调用模式（增强后）
    enhance_edges_with_rf_prediction(args, 
                                     args.project_dir, 
                                     args.output_enhanced_dir, 
                                     args.skip_codon_features, 
                                     args.model_dir)
