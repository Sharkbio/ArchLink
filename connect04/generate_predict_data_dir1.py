# -*- coding: utf-8 -*-
import sys
import pickle
import random
import os
from collections import defaultdict
import numpy as np

# --- Helper Functions for Canonical Directional Pairs ---
def _flip_direction(d):
    """ Flips the direction symbol ('+' or '-'). """
    return '+' if d == '-' else '-'

def _get_canonical_directional_pair_key(c1, d1, c2, d2):
    """
    Standardizes a directional contig pair based on lexicographical order and
    the equivalence rules: (A+, B+) <=> (B-, A-) etc.
    """
    # 确保 contig ID 排序是 c1 <= c2
    if c1 > c2:
        return (c2, _flip_direction(d2), c1, _flip_direction(d1))
    return (c1, d1, c2, d2)

def _get_canonical_non_directional_pair_key(c1, c2):
    """ Standardizes a non-directional contig pair by sorting lexicographically. """
    return tuple(sorted((c1, c2)))

# --- 1. Load Bin Information from FASTA Files ---
def load_contigs_from_fastas(bins_dir):
    """
    Loads contig IDs and their assigned bin IDs from FASTA files within the bins directory.
    """
    contig_to_bin = {}
    print(f"正在从目录 '{bins_dir}' 中加载 bin 和 contig...")
    try:
        for item in os.listdir(bins_dir):
            bin_id = item
            folder_path = os.path.join(bins_dir, item)
            # print(folder_path)
            if os.path.isdir(folder_path):
                fasta_path = os.path.join(folder_path, f"{bin_id}.fasta")
                if os.path.exists(fasta_path):
                    
                    print(f"     > 正在读取 bin: {bin_id}")
                    with open(fasta_path, 'r') as f:
                        for line in f:
                            if line.startswith('>'):
                                # 从 >contigID [其他信息] 中提取 contigID
                                contig_id = line.strip()[1:].split()[0]
                                contig_to_bin[contig_id] = bin_id
                else:
                    print(f"     > 警告: bin '{bin_id}' 的 FASTA 文件未找到，跳过。")
    except FileNotFoundError:
        print(f"错误: Bins 目录 '{bins_dir}' 未找到。")
        sys.exit(1)
    except Exception as e:
        print(f"读取 bins 目录时出错: {e}")
        sys.exit(1)
    return contig_to_bin

# --- 2. Process JUNC File (Modified to read w1 and w2 as floats, without +1 correction) ---
def process_junc_file(filepath, contig_to_bin):
    """
    Processes JUNC file to aggregate A_weight1 (w1) and A_weight2 (w2)
    for canonical directional pairs. Only includes pairs within the same bin.
    """
    # 核心修改: 存储 A_weight1 (w1) 和 A_weight2 (w2)
    directional_junc_data = defaultdict(lambda: {'A_weight1': 0.0, 'A_weight2': 0.0})
    # non_directional_junc_data 现在存储总权重 (w1+w2)
    non_directional_junc_data = defaultdict(lambda: {'A_total_weight': 0.0, 'A_num_junc': 0})
    all_contig_pairs = set()
    all_directional_pairs = set()
    # print(contig_to_bin)
    try:
        with open(filepath, 'r') as f:
            for line in f:
                if line.startswith("JUNC"):
                    parts = line.strip().split()
                    # print(parts)
                    # 预期至少 7 个部分: JUNC c1 d1 c2 d2 w1 w2
                    if len(parts) >= 7:
                        contig1 = parts[1]
                        dir1 = parts[2]
                        contig2 = parts[3]
                        dir2 = parts[4]

                        try:
                            weight1 = float(parts[5]) # 直接读取 w1，作为浮点数
                            weight2 = float(parts[6]) # 直接读取 w2，作为浮点数
                        except ValueError:
                            print(f"警告: 忽略 JUNC 行由于无效的权重值 (非数字): {line.strip()}")
                            continue

                        # 仅处理属于同一个 bin 的 contig 对
                        if contig1 in contig_to_bin and contig2 in contig_to_bin and contig_to_bin[contig1] == contig_to_bin[contig2]:
                            canonical_directional_key = _get_canonical_directional_pair_key(contig1, dir1, contig2, dir2)
                            canonical_non_directional_key = _get_canonical_non_directional_pair_key(contig1, contig2)

                            # 存储方向性特征 A_weight1 (w1) 和 A_weight2 (w2)
                            directional_junc_data[canonical_directional_key]['A_weight1'] += weight1
                            directional_junc_data[canonical_directional_key]['A_weight2'] += weight2

                            # 存储非方向性总权重 (w1+w2) 和次数
                            non_directional_junc_data[canonical_non_directional_key]['A_total_weight'] += (weight1 + weight2)
                            non_directional_junc_data[canonical_non_directional_key]['A_num_junc'] += 1

                            all_contig_pairs.add(canonical_non_directional_key)
                            all_directional_pairs.add(canonical_directional_key)
                        # else: # 如果 contig 不在同一个 bin，则跳过 (静默)
                            
                    else:
                        print(f"警告: 忽略 JUNC 行，预期至少 7 个字段但只有 {len(parts)} 个: {line.strip()}")

    except FileNotFoundError:
        print(f"错误: JUNC 文件 '{filepath}' 未找到。")
        sys.exit(1)
    except Exception as e:
        print(f"读取 JUNC 文件时出错: {e}")
        sys.exit(1)

    return directional_junc_data, non_directional_junc_data, all_contig_pairs, all_directional_pairs

# --- 3. Load C Feature (cosine.pkl) (Modified to read C1 and C2) ---
def load_cosine_feature(filepath):
    """
    Loads the cosine similarity feature file and processes it into
    both canonical directional (C1, C2) and non-directional max (C1_max, C2_max) pairs.
    """
    # 存储 {'C1': score_1, 'C2': score_2}
    directional_cosine_index = {}
    print(f"正在读取余弦特征文件: {filepath}")
    # 存储 {'C1_max': max_c1, 'C2_max': max_c2}，用负无穷初始化以正确追踪最大值
    non_directional_cosine_index = defaultdict(lambda: {'C1_max': -float('inf'), 'C2_max': -float('inf')})

    # --- 修复核心：更新 direction_map 以匹配数据中的 ++, +-, -+, -- 符号 ---
    # 之前使用的是 'head-to-head' 等全称，导致无法匹配数据
    direction_map = {
        '++': ('+', '+'), # head-to-head
        '+-': ('+', '-'), # head-to-tail
        '-+': ('-', '+'), # tail-to-head
        '--': ('-', '-')  # tail-to-tail
    }
    # -------------------------------------------------------------

    try:
        with open(filepath, 'rb') as f:
            cosine_index = pickle.load(f)
        
        print(f"余弦索引包含 {len(cosine_index)} 个非方向对。")

        for (c1, c2), dir_data in cosine_index.items():
            canonical_non_directional_key = _get_canonical_non_directional_pair_key(c1, c2)
            
            # 预期 dir_data.items() 返回 (direction_str, [score_1, score_2])
            for direction_str, cosine_scores in dir_data.items():
                if not isinstance(cosine_scores, list) or len(cosine_scores) != 2:
                    # print(f"警告: 密钥 ({(c1, c2, direction_str)}) 的余弦分数不是包含两个值的列表: {cosine_scores}. 将跳过。")
                    continue
                    
                try:
                    c1_val, c2_val = float(cosine_scores[0]), float(cosine_scores[1])
                except ValueError:
                    # print(f"警告: 密钥 ({(c1, c2, direction_str)}) 的余弦分数无法转换为浮点数. 将跳过。")
                    continue

                # --- 修复：使用更新后的 map，并处理找不到的情况 ---
                d1, d2 = direction_map.get(direction_str, (None, None))
                
                if d1 is None:
                    # 如果方向字符串无法识别（例如不是 ++, +-, -+, --），则跳过，避免生成错误的 key
                    # print(f"跳过无法识别的方向: {direction_str}")
                    continue

                canonical_directional_key = _get_canonical_directional_pair_key(c1, d1, c2, d2)

                # 存储方向性 C1 和 C2
                directional_cosine_index[canonical_directional_key] = {
                    'C1': c1_val,
                    'C2': c2_val
                }
                
                # 追踪非方向性 C1_max 和 C2_max
                current_maxes = non_directional_cosine_index[canonical_non_directional_key]
                current_maxes['C1_max'] = max(current_maxes['C1_max'], c1_val)
                current_maxes['C2_max'] = max(current_maxes['C2_max'], c2_val)

        # 最终清理: 将未被更新的 -inf 转换为 None
        for key, maxes in non_directional_cosine_index.items():
            if maxes['C1_max'] == -float('inf'):
                maxes['C1_max'] = None
            if maxes['C2_max'] == -float('inf'):
                maxes['C2_max'] = None

        print(f"已加载并处理余弦相似度特征文件: {filepath}")
        print(f"提取了 {len(directional_cosine_index)} 个有效的方向性余弦特征。")
        return directional_cosine_index, non_directional_cosine_index
    except FileNotFoundError:
        print(f"错误: 余弦相似度特征文件 '{filepath}' 未找到。")
        return {}, {}
    except Exception as e:
        print(f"读取余弦相似度特征文件时出错: {e}")
        sys.exit(1)

# --- Normalization Helper ---
def normalize_feature_list(values):
    """
    Min-Max 归一化特征列表，保留 None 值。
    """
    valid_values = [v for v in values if v is not None]
    if not valid_values:
        # 如果没有有效值，所有归一化值均为 0.0 (如果非 None) 或 None
        return [0.0 if v is not None else None for v in values]
    min_val = min(valid_values)
    max_val = max(valid_values)
    if max_val == min_val:
        # 如果所有有效值都相同，返回 0.0 (或 None)
        return [0.0 if v is not None else None for v in values]

    # 应用归一化，保留原始的 None 值
    return [(v - min_val) / (max_val - min_val) if v is not None else None for v in values]


# --- 4. Main ---
def main(bins_0_9, bam_graph_dir, raw_cos, outs):
    bins_folder_path = bins_0_9
    junc_file = bam_graph_dir
    c_feature_cosine_file = raw_cos
    output_feature_file = outs

    print("--- Loading bin information from FASTA files ---")
    contig_to_bin = load_contigs_from_fastas(bins_folder_path)
    
    print("\n--- Processing JUNC file (now reading w1 and w2 as floats, removing +1 adjustment) ---")
    directional_junc_data, non_directional_junc_data, all_contig_pairs, all_directional_pairs = process_junc_file(junc_file, contig_to_bin)
    print(f"JUNC file contains {len(all_contig_pairs)} contig pairs from the same bin.")
    print(f"JUNC file contains {len(all_directional_pairs)} directional pairs from the same bin.")

    print("\n--- Loading cosine similarity feature file (reading C1 and C2) ---")
    directional_cosine_index, non_directional_cosine_index = load_cosine_feature(c_feature_cosine_file)

    print("\n--- Extracting features (new feature set) ---")
    all_features_and_bins = {}

    for canonical_directional_key in all_directional_pairs:
        contigA, d1, contigB, d2 = canonical_directional_key
        canonical_non_directional_key = _get_canonical_non_directional_pair_key(contigA, contigB)

        # Verify that the pair belongs to the same bin
        if contigA in contig_to_bin and contigB in contig_to_bin and contig_to_bin[contigA] == contig_to_bin[contigB]:
            bin_id = contig_to_bin[contigA]
            
            # --- Extract non-directional A features ---
            non_dir_a_data = non_directional_junc_data.get(canonical_non_directional_key, {})
            a_num_junc_val = non_dir_a_data.get('A_num_junc', 0)      # integer count
            a_total_weight_val = non_dir_a_data.get('A_total_weight', 0.0)  # float sum of w1+w2
            
            # --- Extract directional A features (A_weight1=w1, A_weight2=w2) ---
            dir_a_data = directional_junc_data.get(canonical_directional_key, {})
            a_weight1_val = dir_a_data.get('A_weight1', 0.0)
            a_weight2_val = dir_a_data.get('A_weight2', 0.0)
            
            # --- Extract directional C features (C1_cosine, C2_cosine) ---
            # 现在 key 的格式已经统一，应该能取到值了
            dir_c_data = directional_cosine_index.get(canonical_directional_key, {})
            c1_cosine_val = dir_c_data.get('C1', None)
            c2_cosine_val = dir_c_data.get('C2', None)

            # --- Extract non-directional C features (C1_cosine_max, C2_cosine_max) ---
            non_dir_c_data = non_directional_cosine_index.get(canonical_non_directional_key, {})
            c1_cosine_max_val = non_dir_c_data.get('C1_max', None)
            c2_cosine_max_val = non_dir_c_data.get('C2_max', None)

            all_features_and_bins[canonical_directional_key] = {
                'binID': bin_id,
                'A_num_junc': a_num_junc_val,
                'A_total_weight': a_total_weight_val,
                'A_weight1': a_weight1_val,
                'A_weight2': a_weight2_val,
                'C1_cosine': c1_cosine_val,
                'C2_cosine': c2_cosine_val,
                'C1_cosine_max': c1_cosine_max_val,
                'C2_cosine_max': c2_cosine_max_val,
            }

    print(f"\nProcessed {len(all_features_and_bins)} directional pairs and extracted features.")

    # --- Normalization ---
    print("\n--- Normalizing features ---")
    
    # 1. Collect numeric feature lists (8 numeric features)
    a_num_junc_values = [v['A_num_junc'] for v in all_features_and_bins.values()]
    a_total_weight_values = [v['A_total_weight'] for v in all_features_and_bins.values()]
    a_weight1_values = [v['A_weight1'] for v in all_features_and_bins.values()]
    a_weight2_values = [v['A_weight2'] for v in all_features_and_bins.values()]
    c1_cosine_max_values = [v['C1_cosine_max'] for v in all_features_and_bins.values()]
    c2_cosine_max_values = [v['C2_cosine_max'] for v in all_features_and_bins.values()]
    c1_cosine_values = [v['C1_cosine'] for v in all_features_and_bins.values()]
    c2_cosine_values = [v['C2_cosine'] for v in all_features_and_bins.values()]

    # 2. Normalize each feature list
    a_num_junc_norm = normalize_feature_list(a_num_junc_values)
    a_total_weight_norm = normalize_feature_list(a_total_weight_values)
    a_weight1_norm = normalize_feature_list(a_weight1_values)
    a_weight2_norm = normalize_feature_list(a_weight2_values)
    c1_cosine_max_norm = normalize_feature_list(c1_cosine_max_values)
    c2_cosine_max_norm = normalize_feature_list(c2_cosine_max_values)
    c1_cosine_norm = normalize_feature_list(c1_cosine_values)
    c2_cosine_norm = normalize_feature_list(c2_cosine_values)

    # 3. Reconstruct normalized feature dictionary
    normalized_features_and_bins = {}
    
    # Create iterators
    a_num_junc_iter = iter(a_num_junc_norm)
    a_total_weight_iter = iter(a_total_weight_norm)
    a_weight1_iter = iter(a_weight1_norm)
    a_weight2_iter = iter(a_weight2_norm)
    c1_cosine_max_iter = iter(c1_cosine_max_norm)
    c2_cosine_max_iter = iter(c2_cosine_max_norm)
    c1_cosine_iter = iter(c1_cosine_norm)
    c2_cosine_iter = iter(c2_cosine_norm)

    for pair, features in all_features_and_bins.items():
        normalized_features_and_bins[pair] = {
            'binID': features['binID'],
            'A_num_junc': next(a_num_junc_iter),
            'A_total_weight': next(a_total_weight_iter),
            'A_weight1': next(a_weight1_iter),
            'A_weight2': next(a_weight2_iter),
            'C1_cosine_max': next(c1_cosine_max_iter),
            'C2_cosine_max': next(c2_cosine_max_iter),
            'C1_cosine': next(c1_cosine_iter),
            'C2_cosine': next(c2_cosine_iter),
        }
    
    # --- 新增分析: 特征存在性计数 ---
    print("\n--- 特征存在性分析 (基于筛选后的数据) ---")
    
    total_pairs = len(normalized_features_and_bins)
    both_a_and_c = 0 # 既有 A 特征也有 C 特征
    only_a = 0       # 只有 A 特征 (即 JUNC 存在，但 Cosine 缺失)
    only_c = 0       # 只有 C 特征 (在当前脚本中不可能发生)

    for pair, features in normalized_features_and_bins.items():
        # 检查 C 特征 (Cosine) 的存在性：只要 C1 或 C2 不为 None 即认为 C 特征存在。
        has_C = features['C1_cosine'] is not None or features['C2_cosine'] is not None

        if has_C:
            both_a_and_c += 1
        else:
            only_a += 1

    # --- 仅修改此处的打印输出 ---
    print(f"筛选后的方向性连接对总数: {total_pairs}")
    print(f"包含 A (JUNC) 和 C (Cosine) 两种特征的连接对数: {both_a_and_c}")
    print(f"只包含 A (JUNC) 单个特征的连接对数: {only_a}")
    print(f"只包含 C (Cosine) 单个特征的连接对数: {only_c}")

    # --- 打印随机归一化数据样本 ---
    print("\n--- 3 个随机归一化数据样本 ---")
    if normalized_features_and_bins:
        all_pairs = list(normalized_features_and_bins.keys())
        for i, pair in enumerate(random.sample(all_pairs, min(3, len(all_pairs)))):
            print(f"\n--- 样本 {i+1} ---")
            print(f"Contig Pair (规范化方向对): {pair}")
            for feature_name, value in normalized_features_and_bins[pair].items():
                print(f"     {feature_name}: {value}")
    else:
        print("数据集为空，无样本可打印。")

    print("\n--- Saving feature file ---")
    try:
        with open(output_feature_file, 'wb') as f:
            pickle.dump(normalized_features_and_bins, f, protocol=pickle.HIGHEST_PROTOCOL)
        print(f"Saved to: {output_feature_file}")
    except Exception as e:
        print(f"Error saving file: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()
