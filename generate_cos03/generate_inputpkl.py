import pandas as pd
import numpy as np
import pickle
import sys
import os
import torch
import esm
from esm import FastaBatchedDataset
from tqdm import tqdm
from multiprocessing import Pool, Manager, get_context
from itertools import repeat
import json

# 设置内存分配优化参数，有助于防止 CUDA 内存碎片化
os.environ["PYTORCH_CUDA_ALLOC_CONF"] = "max_split_size_mb:512"

def read_junc_contig_pairs(junc_file_path):
    """
    读取JUNC文件以获取所有contig配对。
    返回一个包含(contig1, contig2)元组的集合。
    """
    contig_pairs = set()
    try:
        with open(junc_file_path, 'r') as f:
            for line in f:
                parts = line.strip().split()
                # 检查是否是JUNC行且有足够的列
                if len(parts) >= 4 and parts[0] == 'JUNC':
                    contig_pairs.add((parts[1], parts[3]))
    except FileNotFoundError:
        print(f"错误：JUNC文件未找到： {junc_file_path}")
        sys.exit(1)
    except Exception as e:
        print(f"读取JUNC文件时出错: {e}")
        sys.exit(1)
    return contig_pairs

def read_file_A(file_path):
    """
    读取文件A (gene index文件)并返回所有contig及其基因列表。
    格式: contig_id \t gGeneA;gGeneB;...
    """
    contig_dict = {}
    try:
        with open(file_path, 'r') as f:
            for line in f:
                if not line.strip():
                    continue
                parts = line.strip().split('\t')
                if len(parts) != 2:
                    continue
                contig_id, genes = parts
                # 提取基因ID，移除前缀（如'g'或'c'）
                gene_list = [gene[1:] for gene in genes.split(';') if gene]
                contig_dict[contig_id] = gene_list
    except FileNotFoundError:
        print(f"错误：Index文件未找到： {file_path}")
        sys.exit(1)
    return contig_dict

def read_fasta_subset(fasta_file, gene_ids):
    """
    从FASTA文件中读取指定基因ID的序列。
    返回一个字典，键为基因ID，值为其序列。
    """
    gene_sequences = {}
    current_gene = None
    current_sequence = []
    
    gene_ids_set = set(gene_ids)

    def replace_J(sequence):
        # 替换 J 为 I，这是 ESM 模型要求的标准操作
        return sequence.replace('J', 'I')

    try:
        with open(fasta_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    # 处理前一个序列
                    if current_gene and current_sequence:
                        gene_id_from_header = current_gene.split()[0]
                        if gene_id_from_header in gene_ids_set:
                            sequence = "".join(current_sequence)
                            # 过滤掉含有未知氨基酸 'X' 的序列
                            if 'X' not in sequence:
                                gene_sequences[gene_id_from_header] = replace_J(sequence)
                    # 开始新的序列
                    current_gene = line[1:]
                    current_sequence = []
                else:
                    current_sequence.append(line)
            
            # 添加最后一个序列
            if current_gene and current_sequence:
                gene_id_from_header = current_gene.split()[0]
                if gene_id_from_header in gene_ids_set:
                    sequence = "".join(current_sequence)
                    if 'X' not in sequence:
                        gene_sequences[gene_id_from_header] = replace_J(sequence)
    
    except FileNotFoundError:
        print(f"错误：FASTA文件未找到： {fasta_file}")
        sys.exit(1)
    return gene_sequences

def generate_embeddings_worker(gene_ids_chunk, fasta_file, gpu_id, output_dir):
    """
    工作进程函数：在特定 GPU 上为基因块生成 ESM-2 嵌入。
    返回 {gene_id: embedding_vector} 字典。
    """
    torch.cuda.set_device(gpu_id)
    print(f"进程 {gpu_id} 已在 GPU {gpu_id} 上启动.")

    # 在工作进程中临时生成一个FASTA文件，以供 FastaBatchedDataset 使用
    temp_fasta_file = os.path.join(output_dir, f"temp_genes_gpu_{gpu_id}.fasta")
    
    # 读取FASTA文件并写入子集
    # 注意：这里的 gene_ids_chunk 已经只包含 Contig 两端的基因
    filtered_gene_sequences = read_fasta_subset(fasta_file, gene_ids_chunk)
    
    if not filtered_gene_sequences:
        if os.path.exists(temp_fasta_file): os.remove(temp_fasta_file)
        return {}
        
    # 重写临时文件以匹配过滤后的序列
    with open(temp_fasta_file, 'w') as f:
        for gene_id, sequence in filtered_gene_sequences.items():
            f.write(f">{gene_id}\n{sequence}\n")
    
    toks_per_batch = 2000
    dataset = FastaBatchedDataset.from_file(temp_fasta_file)
    batches = dataset.get_batch_indices(toks_per_batch, extra_toks_per_seq=1)

    model_name = "esm2_t36_3B_UR50D"
    # 加载模型并转换为半精度（half precision）以节省显存
    model, vocab = esm.pretrained.load_model_and_alphabet(model_name)
    model = model.cuda().half()
    model.eval()

    data_loader = torch.utils.data.DataLoader(
        dataset,
        collate_fn=vocab.get_batch_converter(),
        batch_sampler=batches
    )

    sequence_representations = {}
    with torch.inference_mode():
        for batch_idx, (labels, strs, toks) in tqdm(enumerate(data_loader), total=len(data_loader), desc=f"GPU {gpu_id} 正在生成嵌入"):
            try:
                toks = toks.cuda(non_blocking=True).long()
                # 截断序列长度到 ESM2 支持的最大值
                toks = toks[:, :12288]
                
                # 提取第33层的残基表示
                results = model(toks, repr_layers=[33], return_contacts=False)
                token_representations = results["representations"][33]
                
                for i, label in enumerate(labels):
                    # 计算平均池化（排除 CLS 和 EOS 令牌）
                    truncate_len = min(12288, len(strs[i]))
                    seq_emb = token_representations[i, 1: truncate_len + 1].mean(0)
                    sequence_representations[label] = seq_emb.cpu().numpy()
            
            except RuntimeError as e:
                # 处理 CUDA 内存不足错误
                print(f"处理 GPU {gpu_id} 上的批次 {batch_idx} 时出错: {e}")
                torch.cuda.empty_cache()
                continue
            
            # 清理 CUDA 内存
            del toks, results, token_representations
            torch.cuda.empty_cache()
    
    # 删除临时文件
    if os.path.exists(temp_fasta_file):
        os.remove(temp_fasta_file)
    
    print(f"进程 {gpu_id} 运行结束.")
    return sequence_representations

def pad_embeddings(embeddings, target_length=12800, pad_at_head=False):
    """
    填充或截断基因嵌入向量到目标长度。
    用于将多个基因嵌入（f5/l5等）连接后统一维度。
    """
    if not embeddings:
        return np.zeros(target_length, dtype=np.float32)
    
    # 过滤掉 None 或空数组
    embeddings = [e for e in embeddings if e is not None and e.size > 0]
    if not embeddings:
        return np.zeros(target_length, dtype=np.float32)

    try:
        embedding = np.concatenate(embeddings)
    except ValueError as e:
        print(f"连接嵌入时出错: {e}")
        return np.zeros(target_length, dtype=np.float32)

    if len(embedding) < target_length:
        # 填充到目标长度
        padding = np.zeros(target_length - len(embedding), dtype=np.float32)
        if pad_at_head:
            # 在头部填充 (L5)
            embedding = np.concatenate([padding, embedding])
        else:
            # 在尾部填充 (F5)
            embedding = np.concatenate([embedding, padding])
    
    # 截断到目标长度
    return embedding[:target_length]

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("用法: python generate_inputpkl_bam.py <JUNC_file> <INDEX_file> <FASTA_file> <output_dir>")
        sys.exit(1)

    junc_file, index_file, fasta_file, output_dir = sys.argv[1:]
    num_gpus = 4 # 可根据您的环境进行调整
    chunk_size = 5000 # 每个文件保存的contig数量

    os.makedirs(output_dir, exist_ok=True)
    
    print("第1步: 读取JUNC文件获取所有contig配对...")
    all_junc_pairs = read_junc_contig_pairs(junc_file)
    # 提取所有JUNC文件中出现的唯一Contig ID
    all_junc_contigs = set()
    for contig1, contig2 in all_junc_pairs:
        all_junc_contigs.add(contig1)
        all_junc_contigs.add(contig2)
    print(f"共找到 {len(all_junc_pairs)} 对contig配对，涉及 {len(all_junc_contigs)} 个唯一contig。")

    print("第2步: 读取索引文件获取contig-基因对应关系...")
    contig_dict = read_file_A(index_file)
    print(f"共找到 {len(contig_dict)} 个contig的基因数据。")
    
    print("第3步-A: 基于'基因数 >= 3'和'出现在JUNC文件中'进行双重过滤...")
    
    # 第一次过滤：基因数 >= 3
    gene_count_filtered_contig_dict = {
        contig_id: genes for contig_id, genes in contig_dict.items() if len(genes) >= 3
    }
    
    # 第二次过滤：必须出现在JUNC配对中
    filtered_contig_dict = {}
    for contig_id, genes in gene_count_filtered_contig_dict.items():
        if contig_id in all_junc_contigs:
            filtered_contig_dict[contig_id] = genes

    print(f"经过基因数量和JUNC文件双重过滤后，剩余 {len(filtered_contig_dict)} 个contig满足条件。")

    # 确定最终用于模型训练的contig配对
    target_contig_pairs = set()
    for contig1, contig2 in all_junc_pairs:
        if contig1 in filtered_contig_dict and contig2 in filtered_contig_dict:
            sorted_pair = tuple(sorted((contig1, contig2)))
            target_contig_pairs.add(sorted_pair)
    
    print(f"最终需要处理的contig配对有 {len(target_contig_pairs)} 对。")

    # --- 核心修改部分：第三重过滤 ---
    print("\n第3步-B: 确认需要生成嵌入的唯一基因（仅收集每个contig的前6和后6个基因）...")
    all_target_genes = set()

    for contig_id, genes in filtered_contig_dict.items():
        # 提取前6个基因
        first_6 = genes[:6]
        # 提取后6个基因。如果 Contig 基因数 < 12，则与 first_6 会有重叠。
        last_6 = genes[-6:]
        
        all_target_genes.update(first_6)
        all_target_genes.update(last_6)
        
    print(f"经过三轮过滤，共找到 {len(all_target_genes)} 个需要生成嵌入的唯一基因（仅包含 contig 两端基因）。")
    # -------------------------------

    all_target_genes_list = list(all_target_genes)
    
    # 将基因列表分块，传递给每个 GPU 进程
    # 确保每个 GPU 至少分配到一个基因块，即使总数小于 num_gpus
    chunk_size_per_gpu = max(1, len(all_target_genes_list) // num_gpus)
    gene_chunks = [all_target_genes_list[i:i + chunk_size_per_gpu] for i in range(0, len(all_target_genes_list), chunk_size_per_gpu)]
    actual_num_processes = min(len(gene_chunks), num_gpus)
    
    print(f"\n第4步: 启动多进程任务 ({actual_num_processes} 个进程)，只为两端基因生成嵌入...")
    
    # 使用 'spawn' 上下文确保 PyTorch/CUDA 兼容性
    ctx = get_context('spawn')
    all_gene_embeddings = {}
    with ctx.Pool(processes=actual_num_processes) as pool:
        try:
            embedding_dicts = pool.starmap(generate_embeddings_worker, zip(
                gene_chunks,
                repeat(fasta_file),
                list(range(actual_num_processes)), # 确保 GPU ID 数量与进程数匹配
                repeat(output_dir)
            ))
            
            # 在主进程中合并所有结果
            for d in embedding_dicts:
                all_gene_embeddings.update(d)
                
            print("所有异步进程已成功完成，嵌入数据已合并。")
        except Exception as e:
            print(f"\n🚨 错误: 多进程期间发生意外错误: {e}")
            pool.terminate()
            pool.join()
            sys.exit(1)
            
    print(f"总共生成了 {len(all_gene_embeddings)} 个基因的嵌入。")

    print("\n第5步: 在主进程中整合 Contig 嵌入数据并分块保存...")
    
    output_data_chunk = {}
    contig_idx_map = {}
    contig_list = list(filtered_contig_dict.keys()) 
    
    for i, contig_id in tqdm(enumerate(contig_list), total=len(contig_list), desc="整合contig数据"):
        gene_list = filtered_contig_dict[contig_id]
        
        # 提取用于聚合的基因列表 (F5, L5, F5-1, L5-1)
        # 因为我们在第3步-B已经确保了这些基因的嵌入都已生成（只要它们不含'X'）
        first_5_genes = gene_list[:5]
        last_5_genes = gene_list[-5:]
        first_6_genes = gene_list[:6]
        f5_minus_1_genes = first_6_genes[1:] # 2-6号基因
        last_6_genes = gene_list[-6:]
        l5_minus_1_genes = last_6_genes[:-1] # 倒数 6-2号基因
        
        # 聚合和填充嵌入
        # 注意：这里的 .get(g) 是安全的，因为如果某个基因（例如 FASTA 中有 'X'）的嵌入没有生成，它会返回 None，
        # 并在 pad_embeddings 中被过滤掉，然后用零向量填充。
        
        # F5
        f5_embedding = pad_embeddings([all_gene_embeddings.get(g) for g in first_5_genes if g in all_gene_embeddings], pad_at_head=False)
        # L5
        l5_embedding = pad_embeddings([all_gene_embeddings.get(g) for g in last_5_genes if g in all_gene_embeddings], pad_at_head=True)
        # F5-1 (第2到6个)
        f5_embedding2 = pad_embeddings([all_gene_embeddings.get(g) for g in f5_minus_1_genes if g in all_gene_embeddings], pad_at_head=False)
        # L5-1 (倒数第2到6个)
        l5_embedding2 = pad_embeddings([all_gene_embeddings.get(g) for g in l5_minus_1_genes if g in all_gene_embeddings], pad_at_head=True)
        
        output_data_chunk[contig_id] = {
            "f5_embedding": f5_embedding,
            "l5_embedding": l5_embedding,
            "f5_embedding2": f5_embedding2,
            "l5_embedding2": l5_embedding2,
            "f_labels": [], # 预留字段
            "l_labels": []  # 预留字段
        }
        
        # 每达到chunk_size或在循环结束时保存
        if (i + 1) % chunk_size == 0 or i == len(contig_list) - 1:
            file_index = i // chunk_size
            file_name = f"contig_embeddings_{file_index:04d}.pkl"
            file_path = os.path.join(output_dir, file_name)
            
            try:
                with open(file_path, 'wb') as f:
                    pickle.dump(output_data_chunk, f, protocol=4) 
                print(f"已将 {len(output_data_chunk)} 个contig数据保存到 {file_path}")
                
                # 记录这些 Contig 所在的 pkl 文件信息
                for c_id in output_data_chunk:
                    contig_idx_map[c_id] = {"file": file_path, "index": file_index}
                
                output_data_chunk = {} # 重置，准备下一个块
            except Exception as e:
                print(f"保存数据块文件时出错: {e}")
                sys.exit(1)

    # 将索引表保存到文件，用于后续查询 Contig 嵌入
    index_file_path = os.path.join(output_dir, "contig_index_map.json")
    with open(index_file_path, 'w') as f:
        json.dump(contig_idx_map, f)
    
    print(f"\n✅ 索引表已成功生成: {index_file_path}")
    print("\n整个流程已完成。")
