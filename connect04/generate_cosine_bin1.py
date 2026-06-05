import pickle
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
import sys
import os
import time
import random
import json
import multiprocessing
from multiprocessing import Manager
import uuid
from collections import defaultdict # 导入 defaultdict
import itertools # 导入 itertools

# =============== 模型部分（来自当前脚本） ===============
class GenePositionalEncoding(nn.Module):
    def __init__(self, d_model, max_genes=5):
        super().__init__()
        self.position_embed = nn.Embedding(max_genes, d_model)

    def forward(self, x):
        # x shape: (seq_len, batch_size, d_model)
        positions = torch.arange(x.size(0), device=x.device).expand(x.size(1), x.size(0))
        return x + self.position_embed(positions).permute(1, 0, 2)

class GeneTransformer(nn.Module):
    def __init__(self, gene_dim, d_model, nhead, num_layers, num_classes, dropout=0.1):
        super().__init__()
        self.gene_proj = nn.Linear(gene_dim, d_model)
        self.pos_encoder = GenePositionalEncoding(d_model)
        encoder_layer = nn.TransformerEncoderLayer(
            d_model=d_model,
            nhead=nhead,
            dim_feedforward=d_model * 4,
            dropout=dropout,
            batch_first=False # 确保 batch_first=False 以匹配 permute
        )
        self.transformer = nn.TransformerEncoder(encoder_layer, num_layers=num_layers)
        
        # --- 修改：使用注意力池化 ---
        self.pool = nn.Linear(d_model, 1)
        
        # --- (新) 修复：添加 cls_head 以匹配预训练模型 ---
        # 我们添加这个只是为了成功加载 state_dict，
        # forward 方法仍将返回 pool 后的嵌入
        self.cls_head = nn.Sequential(
            nn.LayerNorm(d_model),
            nn.Dropout(dropout),
            nn.Linear(d_model, num_classes),
        )

    def forward(self, x):
        # 输入 x 假设 shape: (batch_size, seq_len, gene_dim)
        # 调整为 (seq_len, batch_size, gene_dim) 以适应 TransformerEncoder
        x = x.permute(1, 0, 2) 
        x = self.gene_proj(x)
        x = self.pos_encoder(x)
        x = self.transformer(x) # (seq_len, batch_size, d_model)
        
        # --- 修改：应用注意力池化 ---
        # (seq_len, batch_size, d_model) -> (seq_len, batch_size, 1)
        attn_weights = torch.softmax(self.pool(x).squeeze(-1), dim=0) 
        # (batch_size, d_model)
        pooled = (x * attn_weights.unsqueeze(-1)).sum(dim=0) 
        
        return pooled # 返回池化后的嵌入 (batch_size, d_model)

# --- (新) 辅助函数：反转 f5 嵌入 (来自当前脚本) ---
def reverse_f5_embedding(f5_emb, gene_dim=2560, num_genes=5):
    """
    将 f5 嵌入 (12800,) 反转为 (5, 2560) 的基因顺序。
    """
    try:
        # 确保 f5_emb 是扁平的
        f5_emb_flat = f5_emb.reshape(-1)
        expected_len = gene_dim * num_genes
        if len(f5_emb_flat) != expected_len:
            print(f"警告: f5 嵌入 shape 异常, 期望 {expected_len}, 得到 {len(f5_emb_flat)}")
            # 尝试填充或截断 (这里我们假设它只是末尾的填充/截断)
            if len(f5_emb_flat) > expected_len:
                f5_emb_flat = f5_emb_flat[:expected_len]
            else:
                f5_emb_flat = np.pad(f5_emb_flat, (0, expected_len - len(f5_emb_flat)), 'constant')

        chunks = [f5_emb_flat[i*gene_dim:(i+1)*gene_dim] for i in range(num_genes)]
        # 反转基因顺序
        reversed_f5_emb = np.concatenate(chunks[::-1])
        return reversed_f5_emb.reshape(num_genes, gene_dim)
    except Exception as e:
        print(f"反转嵌入时出错: {e}. 嵌入 shape: {f5_emb.shape}")
        return np.zeros((num_genes, gene_dim), dtype=np.float32)

# --- (新) 辅助函数：计算余弦相似度 (来自当前脚本) ---
def calculate_cosine(e1, e2):
    """
    计算两个 numpy 向量的余弦相似度。
    """
    norm_e1 = np.linalg.norm(e1)
    norm_e2 = np.linalg.norm(e2)
    
    if norm_e1 == 0 or norm_e2 == 0:
        return 0.0
    
    cosine = np.dot(e1, e2) / (norm_e1 * norm_e2)
    # 裁剪到 [-1, 1] 范围以防浮点数误差
    return np.clip(cosine, -1.0, 1.0)


# --- (新) 数据加载函数 (来自原始脚本) ---
def load_contig_embeddings(contig_id, index_map, embedding_dir):
    """
    Load embedding data for a single contig based on the index map.
    """
    if contig_id not in index_map:
        # 注意：在多进程中减少打印，只在关键时刻打印
        # print(f"警告：Contig ID '{contig_id}' 未在索引中找到。")
        return None
    
    file_info = index_map[contig_id]
    file_path = os.path.join(embedding_dir, file_info["file"])

    if not os.path.exists(file_path):
        print(f"错误: 索引中引用的文件 '{file_path}' 未找到。")
        return None
    
    try:
        # 优化：可以考虑缓存加载的文件内容，但目前按原逻辑
        with open(file_path, 'rb') as f:
            data = pickle.load(f)
        return data.get(contig_id) # 返回该 contig 的数据
    except Exception as e:
        print(f"加载文件 '{file_path}' 时出错: {e}")
        return None

def process_junc_file(filepath):
    """
    Read the TSV file, parse all JUNC lines, and extract contig pairs.
    """
    relevant_pairs = set()
    try:
        with open(filepath, 'r') as f:
            for line in f:
                if line.startswith("JUNC"):
                    parts = line.strip().split()
                    if len(parts) >= 4:
                        contig1 = parts[1]
                        contig2 = parts[3]
                        # Normalize the contig pair to treat (A, B) and (B, A) as the same
                        normalized_pair = tuple(sorted((contig1, contig2)))
                        relevant_pairs.add(normalized_pair)
    except FileNotFoundError:
        print(f"错误: JUNC文件 '{filepath}' 未找到。")
        sys.exit(1)
    except Exception as e:
        print(f"读取JUNC文件时发生错误: {e}")
        sys.exit(1)
    
    return relevant_pairs

# --- (新) 合并的 Worker function ---
def worker_process(
    task_pairs, 
    index_map_path, 
    embedding_dir, 
    model_path, # 新增
    device_id, 
    results_dir,
    completed_tasks
):
    """
    子进程工作函数：
    1. 加载模型到指定 GPU
    2. 遍历 Contig 对
    3. 按需加载原始嵌入
    4. 通过模型运行嵌入
    5. 计算余弦相似度
    6. 保存临时结果
    """
    
    # 1. Setup device
    device = torch.device(f"cuda:{device_id}")
    print(f"进程 {os.getpid()} 已分配给设备 {device}")
    
    # 2. Setup constants
    GENE_DIM = 2560
    D_MODEL = 128
    
    # 3. Load Model
    model = GeneTransformer(
        gene_dim=GENE_DIM,
        d_model=D_MODEL,
        nhead=8,
        num_layers=3,
        num_classes=23 # 保留以匹配
    ).to(device)
    
    try:
        model.load_state_dict(torch.load(model_path, map_location=device))
    except Exception as e:
        print(f"[PID {os.getpid()}] 加载模型 '{model_path}' 时出错: {e}")
        return # 此进程无法继续
    model.eval()

    # 4. Load Index Map
    try:
        with open(index_map_path, 'r') as f:
            index_map = json.load(f)
    except Exception as e:
        print(f"[PID {os.getpid()}] 加载索引 '{index_map_path}' 时出错: {e}")
        return

    # 5. Init results and cache
    process_feature_dict = {}
    embedding_cache = {} # 缓存原始加载的数据

    # 6. Iterate tasks
    for contig1, contig2 in task_pairs:
        
        # 7. Load embeddings (using cache)
        try:
            if contig1 not in embedding_cache:
                embedding_cache[contig1] = load_contig_embeddings(contig1, index_map, embedding_dir)
            if contig2 not in embedding_cache:
                embedding_cache[contig2] = load_contig_embeddings(contig2, index_map, embedding_dir)
        except Exception as e:
            print(f"[PID {os.getpid()}] 加载嵌入时出错 {contig1} or {contig2}: {e}")
            completed_tasks.value += 1
            continue # Skip pair

        emb1_data = embedding_cache[contig1]
        emb2_data = embedding_cache[contig2]
        
        if emb1_data is None or emb2_data is None:
            completed_tasks.value += 1
            continue
        
        # 8. Process and get model output for all 8 embeddings
        try:
            # 辅助函数，用于处理单个原始嵌入
            def get_model_emb(raw_emb_data, key_name, is_f5):
                raw_emb = raw_emb_data.get(key_name)
                
                # 检查 None 或全零 (np.all)
                if raw_emb is None or not np.any(raw_emb):
                    # print(f"警告: Contig {raw_emb_data.get('contigid', 'N/A')} 缺少或全零键 {key_name}")
                    return None # 返回 None
                
                if is_f5:
                    processed_emb = reverse_f5_embedding(raw_emb, GENE_DIM, 5)
                else:
                    processed_emb = raw_emb.reshape(5, GENE_DIM)
                
                # 再次检查
                if not np.any(processed_emb):
                    return None

                # 转换为 tensor, 添加 batch dim, 送入模型
                tensor = torch.from_numpy(processed_emb).float().unsqueeze(0).to(device)
                with torch.no_grad():
                    model_output = model(tensor)
                return model_output.squeeze(0).cpu().numpy()

            # 获取所有 8 个模型输出
            emb1 = {
                'l5': get_model_emb(emb1_data, 'l5_embedding', is_f5=False),
                'l52': get_model_emb(emb1_data, 'l5_embedding2', is_f5=False),
                'rev_f5': get_model_emb(emb1_data, 'f5_embedding', is_f5=True),
                'rev_f52': get_model_emb(emb1_data, 'f5_embedding2', is_f5=True)
            }
            emb2 = {
                'l5': get_model_emb(emb2_data, 'l5_embedding', is_f5=False),
                'l52': get_model_emb(emb2_data, 'l5_embedding2', is_f5=False),
                'rev_f5': get_model_emb(emb2_data, 'f5_embedding', is_f5=True),
                'rev_f52': get_model_emb(emb2_data, 'f5_embedding2', is_f5=True)
            }
            
            # 检查是否有任何嵌入处理失败
            if any(v is None for v in emb1.values()) or any(v is None for v in emb2.values()):
                # print(f"跳过 {contig1}/{contig2}，因为存在无效的输入嵌入。")
                completed_tasks.value += 1
                continue

        except Exception as e:
            print(f"[PID {os.getpid()}] 在 {contig1}/{contig2} 上处理模型嵌入时出错: {e}")
            completed_tasks.value += 1
            continue

        # 9. 计算余弦相似度 (使用辅助函数)
        contig_pair_features = {}
        
        # ++: (l5, rev_f52) 和 (l52, rev_f5)
        cos1 = calculate_cosine(emb1['l5'], emb2['rev_f52'])
        cos2 = calculate_cosine(emb1['l52'], emb2['rev_f5'])
        contig_pair_features['++'] = [cos1, cos2]

        # --: (rev_f5, l52) 和 (rev_f52, l5)
        cos1 = calculate_cosine(emb1['rev_f5'], emb2['l52'])
        cos2 = calculate_cosine(emb1['rev_f52'], emb2['l5'])
        contig_pair_features['--'] = [cos1, cos2]

        # +-: (l5, l52) 和 (l52, l5)
        cos1 = calculate_cosine(emb1['l5'], emb2['l52'])
        cos2 = calculate_cosine(emb1['l52'], emb2['l5'])
        contig_pair_features['+-'] = [cos1, cos2]

        # -+: (rev_f5, rev_f52) 和 (rev_f52, rev_f5)
        cos1 = calculate_cosine(emb1['rev_f5'], emb2['rev_f52'])
        cos2 = calculate_cosine(emb1['rev_f52'], emb2['rev_f5'])
        contig_pair_features['-+'] = [cos1, cos2]

        process_feature_dict[tuple(sorted((contig1, contig2)))] = contig_pair_features
        
        # 更新共享计数器
        completed_tasks.value += 1
    
    # 10. 保存临时文件
    temp_filename = f"results_{os.getpid()}_{uuid.uuid4().hex}.pkl"
    temp_filepath = os.path.join(results_dir, temp_filename)
    try:
        with open(temp_filepath, 'wb') as f:
            pickle.dump(process_feature_dict, f)
    except Exception as e:
        print(f"进程 {os.getpid()} 保存临时文件时出错: {e}")
    
    print(f"进程 {os.getpid()} 完成。处理了 {len(task_pairs)} 对中的 {len(process_feature_dict)} 对。")


# =============== 主处理逻辑 (修改) ===============
def main(embedding_dir, junc_file, model_path, output_dir):
    start_time = time.time()
    os.makedirs(output_dir, exist_ok=True)
    
    index_map_path = os.path.join(embedding_dir, "contig_index_map.json")
    
    # 检查模型路径是否存在
    if not os.path.exists(model_path):
        print(f"错误: 指定的模型文件 '{model_path}' 未找到。")
        sys.exit(1)
    
    print("\n--- 载入JUNC文件并生成 contig 对 ---")
    all_relevant_pairs = process_junc_file(junc_file)
    print(f"基于JUNC文件，总共需要处理 {len(all_relevant_pairs)} 个独特的contig对。")

    # 预加载索引图并过滤
    try:
        with open(index_map_path, 'r') as f:
            index_map = json.load(f)
    except FileNotFoundError:
        print(f"错误: 索引文件 '{index_map_path}' 未找到。")
        sys.exit(1)
    
    print("\n--- 过滤 contig 对，仅保留同时存在于JUNC和索引中的对 ---")
    filtered_pairs = []
    for contig1, contig2 in all_relevant_pairs:
        # 数据问题，假设挑选所有的对
        filtered_pairs.append((contig1, contig2))
        
        # if contig1 in index_map and contig2 in index_map:
        #     filtered_pairs.append((contig1, contig2))
    
    print(f"过滤后，实际将处理 {len(filtered_pairs)} 个contig对。")
    if not filtered_pairs:
        print("没有可以处理的contig对，程序退出。")
        sys.exit(0)

    num_gpus = torch.cuda.device_count()
    if num_gpus == 0:
        print("错误: 未检测到可用的GPU。此脚本需要GPU。")
        sys.exit(1)

    print(f"检测到 {num_gpus} 个可用GPU。将使用多进程进行计算。")

    # 将任务对分成等于 GPU 数量的块
    chunks = np.array_split(filtered_pairs, num_gpus)
    
    # 创建 Manager 和共享计数器
    manager = multiprocessing.Manager()
    completed_tasks = manager.Value('i', 0)
    
    # 为进程结果创建临时目录
    temp_results_dir = os.path.join(output_dir, f"temp_results_{uuid.uuid4().hex}")
    os.makedirs(temp_results_dir, exist_ok=True)

    processes = []
    
    # 启动每个进程
    for i, chunk in enumerate(chunks):
        if not chunk.size:
            continue
        p = multiprocessing.Process(
            target=worker_process, 
            args=(
                chunk.tolist(), 
                index_map_path, 
                embedding_dir, 
                model_path, # 传递模型路径
                i, # device_id
                temp_results_dir, 
                completed_tasks
            )
        )
        processes.append(p)
        p.start()

    # --- 进度条显示 ---
    total_tasks = len(filtered_pairs)
    while any(p.is_alive() for p in processes):
        completed = completed_tasks.value
        progress = (completed / total_tasks) * 100
        # 使用 \r 返回行首以动态更新
        sys.stdout.write(f"\r进度: {completed}/{total_tasks} 对已处理 ({progress:.2f}%)")
        sys.stdout.flush()
        time.sleep(1) # 每秒更新
    
    # 确保最终进度条为 100%
    completed = completed_tasks.value # 获取最终值
    progress = (completed / total_tasks) * 100
    sys.stdout.write(f"\r进度: {completed}/{total_tasks} 对已处理 ({progress:.2f}%)\n")
    sys.stdout.flush()
    # --- 进度条结束 ---

    # 等待所有进程完成
    for p in processes:
        p.join()

    # --- 合并所有临时文件的结果 ---
    print("\n--- 所有子进程已完成，正在合并模型特征结果... ---")
    feature_dict = {} # 初始化最终字典
    temp_files = [f for f in os.listdir(temp_results_dir) if f.endswith(".pkl")]
    
    for temp_file in temp_files:
        temp_filepath = os.path.join(temp_results_dir, temp_file)
        try:
            with open(temp_filepath, 'rb') as f:
                feature_part = pickle.load(f)
                feature_dict.update(feature_part)
            os.remove(temp_filepath) # 合并后删除临时文件
        except Exception as e:
            print(f"合并临时文件 '{temp_file}' 时出错: {e}")
    try:
        os.rmdir(temp_results_dir) # 删除临时目录
    except OSError as e:
        print(f"删除临时目录时出错: {e}")

    print(f"总共计算了 {len(feature_dict)} 对 contig 的特征。")

    print("\n--- 打印生成数据集中的任意3个数据样本 ---")
    if feature_dict:
        all_pairs = list(feature_dict.keys())
        num_samples_to_print = min(3, len(all_pairs))
        if num_samples_to_print > 0:
            random_samples = random.sample(all_pairs, num_samples_to_print)
            for i, pair in enumerate(random_samples):
                print(f"\n--- 样本 {i+1} ---")
                print(f"Contig 对: {pair}")
                for direction, cosine_values in feature_dict[pair].items():
                    print(f"    方向 {direction}: [{cosine_values[0]:.4f}, {cosine_values[1]:.4f}]")
        else:
            print("生成的余弦相似度数据集为空，无法打印样本。")
    else:
        print("生成的余弦相似度数据集为空，无法打印样本。")

    # 保存结果到 pkl 文件
    output_file = os.path.join(output_dir, "cosine_model_features.pkl")

    print(f"\nSaving results to {output_file}...")
    try:
        with open(output_file, 'wb') as f:
            pickle.dump(feature_dict, f, protocol=4)
        
        print(f"Total time: {time.time()-start_time:.2f}s")
        total_directional_combinations = len(feature_dict) * 4
        total_scores = total_directional_combinations * 2
        print(f"Saved features for {len(feature_dict)} unique contig pairs.")
        print(f"Total {total_directional_combinations} directional combinations, {total_scores} total cosine scores.")
    except Exception as e:
        print(f"保存特征文件时发生错误: {e}")
        sys.exit(1)


if __name__ == "__main__":
    # 必须在 if __name__ == "__main__": 块内设置多进程启动方法
    try:
        multiprocessing.set_start_method("spawn", force=True)
    except RuntimeError as e:
        # 忽略已经设置的错误
        if "context has already been set" not in str(e):
            print(f"警告：无法设置 multiprocessing start method: {e}")
        pass

    if len(sys.argv) != 5:
        print("用法: python 01.generate_cosine_model_features_v2.py <embedding_dir> <junc_file.tsv> <model_path.pth> <output_dir>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
