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
from collections import defaultdict 
import itertools 
GENE_DIM = 2560
D_MODEL = 128
NUM_CLASSES = 23
# ==============================================================================
# Model Definition (Modified to output Softmax Probabilities)
# ==============================================================================

class GenePositionalEncoding(nn.Module):
    """
    Gene-aware Positional Encoding.
    基因感知的位置编码。
    """
    def __init__(self, d_model, max_genes=5):
        super().__init__()
        self.position_embed = nn.Embedding(max_genes, d_model)

    def forward(self, x):
        # x shape: (seq_len, batch_size, d_model)
        # Positions tensor: (batch_size, seq_len)
        positions = torch.arange(x.size(0), device=x.device).expand(x.size(1), x.size(0))
        # Add position embedding to input, ensuring correct dimension alignment via permute
        return x + self.position_embed(positions).permute(1, 0, 2)

class GeneTransformer(nn.Module):
    """
    Transformer model for gene embeddings.
    MODIFIED: Now returns the Softmax probability vector as the feature vector.
    已修改：现在返回 Softmax 概率向量作为特征向量。
    """
    def __init__(self, gene_dim, d_model, nhead, num_layers, num_classes, dropout=0.1):
        super().__init__()
        self.gene_proj = nn.Linear(gene_dim, d_model)
        self.pos_encoder = GenePositionalEncoding(d_model)
        
        encoder_layer = nn.TransformerEncoderLayer(
            d_model=d_model,
            nhead=nhead,
            dim_feedforward=d_model * 4,
            dropout=dropout,
            batch_first=False 
        )
        self.transformer = nn.TransformerEncoder(encoder_layer, num_layers=num_layers)
        
        # Simple linear layer for attention pooling 
        self.pool = nn.Linear(d_model, 1)
        
        # Classification Head (Outputting Logits)
        self.cls_head = nn.Sequential(
            nn.LayerNorm(d_model),
            nn.Dropout(dropout),
            nn.Linear(d_model, num_classes) 
        )

    def forward(self, x):
        # Input x shape assumed: (batch_size, seq_len, gene_dim)
        # Permute to (seq_len, batch_size, gene_dim)
        x = x.permute(1, 0, 2) 
        x = self.gene_proj(x)
        x = self.pos_encoder(x)
        x = self.transformer(x) # (seq_len, batch_size, d_model)
        
        # Apply Attention Pooling
        # (seq_len, batch_size, d_model) -> (seq_len, batch_size, 1)
        attn_weights = torch.softmax(self.pool(x).squeeze(-1), dim=0) 
        # Weighted sum: (batch_size, d_model)
        pooled = (x * attn_weights.unsqueeze(-1)).sum(dim=0) 
        
        # 1. Get Logits
        logits = self.cls_head(pooled) # (batch_size, num_classes)

        # --- MODIFIED: Return Probabilities (Softmax output) ---
        # 2. Apply Softmax to get the probability vector
        # Apply softmax along the feature dimension (dim=-1)
        probabilities = F.softmax(logits, dim=-1)
        
        return probabilities # Returns the probability vector for cosine similarity

# ==============================================================================
# Auxiliary Functions (No Change)
# ==============================================================================

def reverse_f5_embedding(f5_emb, gene_dim=2560, num_genes=5):
    """
    Reverses the gene order of a flattened f5 embedding to (5, 2560).
    反转展平的 f5 嵌入的基因顺序。
    """
    try:
        f5_emb_flat = f5_emb.reshape(-1)
        expected_len = gene_dim * num_genes
        if len(f5_emb_flat) != expected_len:
            # Padding or truncating logic (simplified for robustness)
            if len(f5_emb_flat) > expected_len:
                f5_emb_flat = f5_emb_flat[:expected_len]
            else:
                f5_emb_flat = np.pad(f5_emb_flat, (0, expected_len - len(f5_emb_flat)), 'constant')

        chunks = [f5_emb_flat[i*gene_dim:(i+1)*gene_dim] for i in range(num_genes)]
        # Reverse the list of chunks and concatenate
        reversed_f5_emb = np.concatenate(chunks[::-1])
        return reversed_f5_emb.reshape(num_genes, gene_dim)
    except Exception as e:
        # print(f"Error reversing embedding: {e}. Embedding shape: {f5_emb.shape}")
        return np.zeros((num_genes, gene_dim), dtype=np.float32)

def calculate_cosine(e1, e2):
    """
    Calculates the cosine similarity between two numpy vectors.
    计算两个 numpy 向量之间的余弦相似度。
    """
    norm_e1 = np.linalg.norm(e1)
    norm_e2 = np.linalg.norm(e2)
    
    if norm_e1 == 0 or norm_e2 == 0:
        return 0.0
    
    cosine = np.dot(e1, e2) / (norm_e1 * norm_e2)
    # Clip to [-1, 1] range to prevent floating point errors
    return np.clip(cosine, -1.0, 1.0)


def load_contig_embeddings(contig_id, index_map, embedding_dir):
    """
    Load embedding data for a single contig based on the index map.
    根据索引图加载单个 contig 的嵌入数据。
    """
    if contig_id not in index_map:
        return None
    
    file_info = index_map[contig_id]
    file_path = os.path.join(embedding_dir, file_info["file"])

    if not os.path.exists(file_path):
        print(f"Error: File '{file_path}' referenced in index not found.")
        return None
    
    try:
        with open(file_path, 'rb') as f:
            data = pickle.load(f)
        return data.get(contig_id) 
    except Exception as e:
        print(f"Error loading file '{file_path}': {e}")
        return None

def process_junc_file(filepath):
    """
    Read the TSV file, parse all JUNC lines, and extract contig pairs.
    读取 TSV 文件，解析所有 JUNC 行，提取 contig 对。
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
        print(f"Error: JUNC file '{filepath}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading JUNC file: {e}")
        sys.exit(1)
    
    return relevant_pairs

# ==============================================================================
# Worker Process (Now returns Softmax Probabilities)
# ==============================================================================

def worker_process(
    task_pairs, 
    index_map_path, 
    embedding_dir, 
    model_path, 
    device_id, 
    results_dir,
    completed_tasks
):
    """
    Worker function for child processes.
    子进程工作函数。
    """
    
    # 1. Setup device
    device = torch.device(f"cuda:{device_id}")
    print(f"Process {os.getpid()} assigned to device {device}")
    
    # 2. Setup constants (Retain original script constants)
    GENE_DIM = 2560
    D_MODEL = 128
    NUM_CLASSES = 23 # Assuming the model was trained with 23 classes
    
    # 3. Load Model
    model = GeneTransformer(
        gene_dim=GENE_DIM,
        d_model=D_MODEL,
        nhead=8,
        num_layers=3,
        num_classes=NUM_CLASSES 
    ).to(device)
    
    try:
        # Load the weights
        model.load_state_dict(torch.load(model_path, map_location=device))
    except Exception as e:
        print(f"[PID {os.getpid()}] Error loading model '{model_path}': {e}")
        return 
    model.eval()

    # 4. Load Index Map
    try:
        with open(index_map_path, 'r') as f:
            index_map = json.load(f)
    except Exception as e:
        print(f"[PID {os.getpid()}] Error loading index '{index_map_path}': {e}")
        return

    # 5. Init results and cache
    process_feature_dict = {}
    embedding_cache = {}

    # 6. Iterate tasks
    for contig1, contig2 in task_pairs:
        
        # 7. Load embeddings (using cache)
        try:
            if contig1 not in embedding_cache:
                embedding_cache[contig1] = load_contig_embeddings(contig1, index_map, embedding_dir)
            if contig2 not in embedding_cache:
                embedding_cache[contig2] = load_contig_embeddings(contig2, index_map, embedding_dir)
        except Exception as e:
            print(f"[PID {os.getpid()}] Error loading embeddings for {contig1} or {contig2}: {e}")
            completed_tasks.value += 1
            continue 

        emb1_data = embedding_cache[contig1]
        emb2_data = embedding_cache[contig2]

        if emb1_data is None or emb2_data is None:
            completed_tasks.value += 1
            continue
        
        # 8. Process and get model output (now probabilities/features)
        try:
            # Helper function to process a single raw embedding
            def get_model_emb(raw_emb_data, key_name, is_f5):
                raw_emb = raw_emb_data.get(key_name)
                
                if raw_emb is None or not np.any(raw_emb):
                    return None 
                
                if is_f5:
                    processed_emb = reverse_f5_embedding(raw_emb, GENE_DIM, 5)
                else:
                    processed_emb = raw_emb.reshape(5, GENE_DIM)
                
                if not np.any(processed_emb):
                    return None

                # Convert to tensor, add batch dim, feed into model
                tensor = torch.from_numpy(processed_emb).float().unsqueeze(0).to(device)
                with torch.no_grad():
                    # model_output is now the Softmax probability vector (1, NUM_CLASSES)
                    model_output = model(tensor) 
                return model_output.squeeze(0).cpu().numpy() # (NUM_CLASSES,)

            # Get all 8 model outputs
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
            
            if any(v is None for v in emb1.values()) or any(v is None for v in emb2.values()):
                completed_tasks.value += 1
                continue

        except Exception as e:
            print(f"[PID {os.getpid()}] Error processing model embeddings on {contig1}/{contig2}: {e}")
            completed_tasks.value += 1
            continue

        # 9. Calculate Cosine Similarity (using the new probability features)
        contig_pair_features = {}
        
        # ++: (l5, rev_f52) and (l52, rev_f5)
        cos1 = calculate_cosine(emb1['l5'], emb2['rev_f52'])
        cos2 = calculate_cosine(emb1['l52'], emb2['rev_f5'])
        contig_pair_features['++'] = [cos1, cos2]

        # --: (rev_f5, l52) and (rev_f52, l5)
        cos1 = calculate_cosine(emb1['rev_f5'], emb2['l52'])
        cos2 = calculate_cosine(emb1['rev_f52'], emb2['l5'])
        contig_pair_features['--'] = [cos1, cos2]

        # +-: (l5, l52) and (l52, l5)
        cos1 = calculate_cosine(emb1['l5'], emb2['l52'])
        cos2 = calculate_cosine(emb1['l52'], emb2['l5'])
        contig_pair_features['+-'] = [cos1, cos2]

        # -+: (rev_f5, rev_f52) and (rev_f52, rev_f5)
        cos1 = calculate_cosine(emb1['rev_f5'], emb2['rev_f52'])
        cos2 = calculate_cosine(emb1['rev_f52'], emb2['rev_f5'])
        contig_pair_features['-+'] = [cos1, cos2]

        process_feature_dict[tuple(sorted((contig1, contig2)))] = contig_pair_features
        
        # Update shared counter
        completed_tasks.value += 1
    
    # 10. Save temporary file
    temp_filename = f"results_{os.getpid()}_{uuid.uuid4().hex}.pkl"
    temp_filepath = os.path.join(results_dir, temp_filename)
    try:
        with open(temp_filepath, 'wb') as f:
            pickle.dump(process_feature_dict, f)
    except Exception as e:
        print(f"Process {os.getpid()} failed to save temp file: {e}")
    
    print(f"Process {os.getpid()} finished. Processed {len(task_pairs)} pairs, saved features for {len(process_feature_dict)}.")


# ==============================================================================
# Main Processing Logic (Output file name updated)
# ==============================================================================
def main(embedding_dir, junc_file, model_path, output_dir):
    start_time = time.time()
    os.makedirs(output_dir, exist_ok=True)
    
    index_map_path = os.path.join(embedding_dir, "contig_index_map.json")
    
    if not os.path.exists(model_path):
        print(f"Error: Specified model file '{model_path}' not found.")
        sys.exit(1)
    
    print("\n--- Loading JUNC file and generating contig pairs ---")
    all_relevant_pairs = process_junc_file(junc_file)
    print(f"Total unique contig pairs from JUNC file: {len(all_relevant_pairs)}")

    try:
        with open(index_map_path, 'r') as f:
            index_map = json.load(f)
    except FileNotFoundError:
        print(f"Error: Index file '{index_map_path}' not found.")
        sys.exit(1)
    
    print("\n--- Filtering contig pairs present in both JUNC and index ---")
    filtered_pairs = []
    for contig1, contig2 in all_relevant_pairs:
        if contig1 in index_map and contig2 in index_map:
            filtered_pairs.append((contig1, contig2))
    
    print(f"Filtered pairs to process: {len(filtered_pairs)}")
    if not filtered_pairs:
        print("No contig pairs left to process. Exiting.")
        sys.exit(0)

    num_gpus = torch.cuda.device_count()
    if num_gpus == 0:
        print("Error: No available GPUs detected. This script requires a GPU.")
        sys.exit(1)

    print(f"Detected {num_gpus} available GPUs. Using multiprocessing for computation.")

    # Convert to array for splitting, then back to list for worker args
    chunks = np.array_split(filtered_pairs, num_gpus)
    
    manager = multiprocessing.Manager()
    completed_tasks = manager.Value('i', 0)
    
    temp_results_dir = os.path.join(output_dir, f"temp_results_{uuid.uuid4().hex}")
    os.makedirs(temp_results_dir, exist_ok=True)

    processes = []
    
    for i, chunk in enumerate(chunks):
        if not chunk.size:
            continue
        p = multiprocessing.Process(
            target=worker_process, 
            args=(
                chunk.tolist(), 
                index_map_path, 
                embedding_dir, 
                model_path, 
                i, # device_id
                temp_results_dir, 
                completed_tasks
            )
        )
        processes.append(p)
        p.start()

    # --- Progress Bar Display ---
    total_tasks = len(filtered_pairs)
    while any(p.is_alive() for p in processes):
        completed = completed_tasks.value
        progress = (completed / total_tasks) * 100
        sys.stdout.write(f"\rProgress: {completed}/{total_tasks} pairs processed ({progress:.2f}%)")
        sys.stdout.flush()
        time.sleep(1) 
    
    completed = completed_tasks.value 
    progress = (completed / total_tasks) * 100
    sys.stdout.write(f"\rProgress: {completed}/{total_tasks} pairs processed ({progress:.2f}%)\n")
    sys.stdout.flush()
    # --- End of Progress Bar ---

    for p in processes:
        p.join()

    # --- Merge Results ---
    print("\n--- All child processes finished, merging model features... ---")
    feature_dict = {} 
    temp_files = [f for f in os.listdir(temp_results_dir) if f.endswith(".pkl")]
    
    for temp_file in temp_files:
        temp_filepath = os.path.join(temp_results_dir, temp_file)
        try:
            with open(temp_filepath, 'rb') as f:
                feature_part = pickle.load(f)
                feature_dict.update(feature_part)
            os.remove(temp_filepath) 
        except Exception as e:
            print(f"Error merging temp file '{temp_file}': {e}")
    try:
        os.rmdir(temp_results_dir) 
    except OSError as e:
        print(f"Error removing temporary directory: {e}")

    print(f"Total features calculated for {len(feature_dict)} contig pairs.")

    print("\n--- Printing 3 random samples from the generated dataset ---")
    if feature_dict:
        all_pairs = list(feature_dict.keys())
        num_samples_to_print = min(3, len(all_pairs))
        if num_samples_to_print > 0:
            random_samples = random.sample(all_pairs, num_samples_to_print)
            for i, pair in enumerate(random_samples):
                print(f"\n--- Sample {i+1} ---")
                print(f"Contig Pair: {pair}")
                for direction, cosine_values in feature_dict[pair].items():
                    print(f"    Direction {direction}: [{cosine_values[0]:.4f}, {cosine_values[1]:.4f}]")
        else:
            print("Generated cosine similarity dataset is empty. Cannot print samples.")
    else:
        print("Generated cosine similarity dataset is empty. Cannot print samples.")

    # Save results to pkl file
    output_file = os.path.join(output_dir, "cosine_model_features_softmax.pkl")

    print(f"\nSaving results to {output_file}...")
    try:
        with open(output_file, 'wb') as f:
            pickle.dump(feature_dict, f, protocol=4)
        
        print(f"Total time: {time.time()-start_time:.2f}s")
        total_directional_combinations = len(feature_dict) * 4
        total_scores = total_directional_combinations * 2
        print(f"Saved features for {len(feature_dict)} unique contig pairs.")
        print(f"Total {total_directional_combinations} directional combinations, {total_scores} total cosine scores.")
        print(f"*** Feature Vector Dimension is now NUM_CLASSES (i.e., {NUM_CLASSES}), representing SOFTMAX PROBABILITIES ***")

    except Exception as e:
        print(f"Error saving feature file: {e}")
        sys.exit(1)


if __name__ == "__main__":
    try:
        # Set start method for multiprocessing compatibility
        multiprocessing.set_start_method("spawn", force=True)
    except RuntimeError as e:
        if "context has already been set" not in str(e):
            print(f"Warning: Could not set multiprocessing start method: {e}")
        pass

    if len(sys.argv) != 5:
        print("Usage: python gene_feature_extractor_softmax.py <embedding_dir> <junc_file.tsv> <model_path.pth> <output_dir>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])