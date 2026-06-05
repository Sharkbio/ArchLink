import pandas as pd
import re
import os
import sys
from typing import List

def get_unique_top_methods(df: pd.DataFrame) -> List[str]:
    """
    根据新的逻辑选择最佳的 6 个方法 ID: 
    1. sum 最高的 3 个。
    2. sum_cont5 最高的 3 个 (必须是与 sum 选出的不同的方法)。
    
    Args:
        df (pd.DataFrame): 包含 'sum', 'sum_cont5', 'Binning_method' 的 DataFrame。
        
    Returns:
        List[str]: 包含最多 6 个唯一方法 ID 的列表。
    """
    
    # 1. 忽略组概念，对整个数据集排序 (按 sum)
    df_sum_sorted = df.sort_values(by=['sum', 'sum_cont5'], ascending=[False, False]).copy()
    
    # 2. 选取 sum 最高的 3 个
    sum_top_3: List[str] = df_sum_sorted.head(3)['Binning_method'].tolist()
    
    if not sum_top_3:
        return []

    # 3. 选取 sum_cont5 最高的 3 个 (必须是与 sum 选出的不同的方法)
    # 按 'sum_cont5' 降序排序，次要排序键为 'sum'
    df_cont5_sorted = df.sort_values(by=['sum_cont5', 'sum'], ascending=[False, False]).copy()
    
    cont5_unique_top_3: List[str] = []
    
    # 用于快速查找 sum 已经选过的方法
    sum_selected_set = set(sum_top_3)
    
    # 遍历 sum_cont5 排序后的列表，选取独有的方法
    for method_id in df_cont5_sorted['Binning_method']:
        if method_id not in sum_selected_set:
            cont5_unique_top_3.append(method_id)
            if len(cont5_unique_top_3) >= 3:
                break
    
    # 4. 组合最终方法列表 (sum_top_3 在前，保持 sum 中最好的方法排第一)
    final_methods = sum_top_3 + cont5_unique_top_3
    
    return final_methods[:6]


def get_best_and_other_ids_for_shell(estimate_file_path: str):
    """
    读取 estimate_res.txt 文件，进行排序，并应用新的 ID 选择逻辑：
    1. 选取 sum 最高的 3 个方法 ID。
    2. 选取 sum_cont5 最高的 3 个独特的 ID。
    3. 将选定的 ID 转换为简化的格式供 shell 脚本使用。
    4. 打印最佳 ID (sum 排名第一) 和其他 ID。
    """
    if not os.path.exists(estimate_file_path):
        print(f"ERROR: estimate_res.txt not found - {estimate_file_path}", file=sys.stderr)
        sys.exit(1)

    try:
        # 使用 low_memory=False 避免 DtypeWarning
        df = pd.read_csv(estimate_file_path, sep='\t', low_memory=False) 
    except Exception as e:
        print(f"ERROR: Error reading {estimate_file_path}: {e}", file=sys.stderr)
        sys.exit(1)

    required_cols = ['sum', 'sum_cont5', 'Binning_method']
    if not all(col in df.columns for col in required_cols):
        print("ERROR: Required columns ('sum', 'sum_cont5', or 'Binning_method') not found in estimate_res.txt.", file=sys.stderr)
        sys.exit(1)

    # --- Step 1 & 2 & 3: ID 选择逻辑 ---
    final_methods = get_unique_top_methods(df)

    if not final_methods:
        print("ERROR: Could not select any binning method IDs.", file=sys.stderr)
        sys.exit(1)
        
    print(f"INFO: Successfully selected {len(final_methods)} methods: {final_methods}", file=sys.stderr)

    # --- Step 4: 转换为简化格式 ---
    # 示例：0.15_100_1_100
    conversion_regex = re.compile(r'Leiden_bandwidth_(\d+\.\d+)_res_maxedges(\d+)respara_(\d+)_partgraph_ratio_(\d+)\.tsv')

    converted_ids = []
    for original_id in final_methods:
        match = conversion_regex.search(original_id)
        if match:
            # 组索引：(1: bandwidth, 2: maxedges, 3: respara, 4: partgraph_ratio)
            converted_ids.append(f"{match.group(1)}_{match.group(2)}_{match.group(3)}_{match.group(4)}")
        else:
            # 对于不匹配特定模式的 ID 进行回退处理
            converted_ids.append(original_id.replace('.tsv', '').replace('Leiden_bandwidth_', ''))
            
            
    # 5. 打印结果供 shell 脚本使用
    # # 在第一行打印最佳 ID (sum 排名第一的那个)
    # print(converted_ids[0])
    # # 在第二行打印其余 ID，以空格分隔
    # if len(converted_ids) > 1:
    #     print(" ".join(converted_ids[1:]))
    # else:
    #     print("") # 确保打印第二行，即使为空
    if len(converted_ids) > 1:
        return converted_ids[0],converted_ids[1:]
    else:
        return converted_ids[0],None
        

if __name__ == "__main__":
    # 脚本期望文件路径作为第二个参数
    if len(sys.argv) != 2:
        print("Usage: python get_ids_helper.py <path_to_estimate_file>", file=sys.stderr)
        sys.exit(1)
    
    # b_val 参数被保留用于兼容性
    estimate_file_path = sys.argv[1]
    
    get_best_and_other_ids_for_shell(estimate_file_path)
