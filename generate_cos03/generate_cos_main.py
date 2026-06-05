from .depth_aver import *
from .generate_bins import *
from .contigs_index import *
from .clean_faa import *
from . import generate_cosine_bin1
from .generate_inputpkl import *
from . import pprodigal
from typing import List
import glob
import os
import subprocess
import sys


def require_existing_path(*candidates):
    for candidate in candidates:
        if os.path.exists(candidate):
            return candidate
    raise FileNotFoundError(f"Required file not found. Checked: {', '.join(candidates)}")


def merge_fasta_files(binning_dir, output_fasta):
    """Merge all FASTA files (.fa/.fasta) inside subdirectories."""
    if os.path.exists(output_fasta):
        print(f"Merged FASTA already exists: {output_fasta}")
        return

    print("Merging all FASTA files...")

    # Create an empty output file
    with open(output_fasta, "w"):
        pass

    # Collect all .fa / .fasta files
    fasta_files = []
    for root, _, files in os.walk(binning_dir):
        for name in files:
            if name.endswith(".fa") or name.endswith(".fasta"):
                fasta_files.append(os.path.join(root, name))

    # Append contents
    with open(output_fasta, "a") as outfile:
        for fasta_file in fasta_files:
            print(f"Processing: {fasta_file}")
            with open(fasta_file, "r") as f:
                outfile.write(f.read())

    print(f"FASTA merge finished: {output_fasta}")


def generate_inputpkl_main(args, bam_graph_path, bin_index_path, out_path, final_input_pkl_path):
    """Main pipeline to generate input PKL files for model training."""

    junc_file, index_file, fasta_file = bam_graph_path, bin_index_path, out_path
    output_dir = final_input_pkl_path
    num_gpus = args.num_gpus
    chunk_size = args.chunk_size

    os.makedirs(output_dir, exist_ok=True)

    print("Step 1: Reading JUNC file to collect contig pairs...")
    all_junc_pairs = read_junc_contig_pairs(junc_file)

    all_junc_contigs = set()
    for c1, c2 in all_junc_pairs:
        all_junc_contigs.add(c1)
        all_junc_contigs.add(c2)

    print(f"Found {len(all_junc_pairs)} contig pairs involving {len(all_junc_contigs)} unique contigs.")

    print("Step 2: Reading contig-gene index...")
    contig_dict = read_file_A(index_file)
    print(f"Loaded gene info for {len(contig_dict)} contigs.")

    print("Step 3A: Double filtering: (1) ≥3 genes (2) must appear in JUNC...")
    gene_count_filtered = {cid: genes for cid, genes in contig_dict.items() if len(genes) >= 3}

    filtered_contig_dict = {cid: genes for cid, genes in gene_count_filtered.items() if cid in all_junc_contigs}
    print(f"{len(filtered_contig_dict)} contigs remain after filtering.")

    # Final eligible contig pairs
    target_contig_pairs = set()
    for c1, c2 in all_junc_pairs:
        if c1 in filtered_contig_dict and c2 in filtered_contig_dict:
            target_contig_pairs.add(tuple(sorted((c1, c2))))

    print(f"{len(target_contig_pairs)} valid contig pairs will be processed.")

    print("\nStep 3B: Collecting unique target genes (only first 6 and last 6 genes per contig)...")
    all_target_genes = set()

    for cid, genes in filtered_contig_dict.items():
        first6 = genes[:6]
        last6 = genes[-6:]
        all_target_genes.update(first6)
        all_target_genes.update(last6)

    print(f"Collected {len(all_target_genes)} unique target genes.")

    all_target_genes_list = list(all_target_genes)

    # Split into GPU tasks
    chunk_size_per_gpu = max(1, len(all_target_genes_list) // num_gpus)
    gene_chunks = [all_target_genes_list[i:i + chunk_size_per_gpu]
                   for i in range(0, len(all_target_genes_list), chunk_size_per_gpu)]

    actual_num_processes = min(len(gene_chunks), num_gpus)

    print(f"\nStep 4: Launching multiprocessing with {actual_num_processes} workers...")

    ctx = get_context('spawn')
    all_gene_embeddings = {}

    with ctx.Pool(processes=actual_num_processes) as pool:
        try:
            embedding_dicts = pool.starmap(
                generate_embeddings_worker,
                zip(gene_chunks,
                    repeat(fasta_file),
                    list(range(actual_num_processes)),
                    repeat(output_dir))
            )

            for d in embedding_dicts:
                all_gene_embeddings.update(d)

            print("All embedding workers finished.")
        except Exception as e:
            print(f"Error during multiprocessing: {e}")
            pool.terminate()
            pool.join()
            sys.exit(1)

    print(f"Generated embeddings for {len(all_gene_embeddings)} genes.")

    print("\nStep 5: Building contig-level embeddings and saving in chunks...")

    output_data_chunk = {}
    contig_idx_map = {}
    contig_list = list(filtered_contig_dict.keys())

    for i, cid in tqdm(enumerate(contig_list), total=len(contig_list), desc="Building contig embeddings"):
        genes = filtered_contig_dict[cid]

        # Gene groups:
        first5 = genes[:5]
        last5 = genes[-5:]
        first6 = genes[:6]
        f5_minus_1 = first6[1:]
        last6 = genes[-6:]
        l5_minus_1 = last6[:-1]

        # Safe embedding padding
        f5_emb = pad_embeddings([all_gene_embeddings.get(g) for g in first5 if g in all_gene_embeddings], pad_at_head=False)
        l5_emb = pad_embeddings([all_gene_embeddings.get(g) for g in last5 if g in all_gene_embeddings], pad_at_head=True)
        f5_emb2 = pad_embeddings([all_gene_embeddings.get(g) for g in f5_minus_1 if g in all_gene_embeddings], pad_at_head=False)
        l5_emb2 = pad_embeddings([all_gene_embeddings.get(g) for g in l5_minus_1 if g in all_gene_embeddings], pad_at_head=True)

        output_data_chunk[cid] = {
            "f5_embedding": f5_emb,
            "l5_embedding": l5_emb,
            "f5_embedding2": f5_emb2,
            "l5_embedding2": l5_emb2,
            "f_labels": [],
            "l_labels": []
        }

        # Save chunk
        if (i + 1) % chunk_size == 0 or i == len(contig_list) - 1:
            file_index = i // chunk_size
            file_name = f"contig_embeddings_{file_index:04d}.pkl"
            file_path = os.path.join(output_dir, file_name)

            try:
                with open(file_path, 'wb') as f:
                    pickle.dump(output_data_chunk, f, protocol=4)

                print(f"Saved {len(output_data_chunk)} contigs to {file_path}")

                for c_id in output_data_chunk:
                    contig_idx_map[c_id] = {"file": file_path, "index": file_index}

                output_data_chunk = {}
            except Exception as e:
                print(f"Error saving chunk file: {e}")
                sys.exit(1)

    # Save index mapping
    index_file_path = os.path.join(output_dir, "contig_index_map.json")
    with open(index_file_path, 'w') as f:
        json.dump(contig_idx_map, f)

    print(f"\nIndex file generated: {index_file_path}")
    print("Pipeline completed.")

def main(args):
    # --- 变量定义和路径设置 ---
    # 假设 args 中包含了所有必需的路径信息，包括 bam_file_path 作为 BAM 文件的源目录
    
    GLOBAL_MERGED_GRAPH = os.path.join(args.output_path, 'bam.graph')
    GLOBAL_MERGED_BAM = os.path.join(args.output_path, 'merged.bam')
    GLOBAL_MERGED_DEPTH = os.path.join(args.output_path, 'merged.depth')
    SAMTOOLS_BIN = args.samtools_bin if hasattr(args, 'samtools_bin') else "samtools"
    
    # 新增 BAM 源目录变量，并替换掉不再使用的 SRR_FOLDERS 和 GLOBAL_COMBINE_BAM_DIR
    BAM_SOURCE_DIR = args.bam_file 
    
    # 检查最终图文件是否存在
    print(f"\n# 检查最终图文件是否存在: {GLOBAL_MERGED_GRAPH}")
    if os.path.exists(GLOBAL_MERGED_GRAPH):
        print(f"  最终目标文件 {GLOBAL_MERGED_GRAPH} 已存在，跳过全局合并BAM和生成图步骤")
    else:
        # --- Shell 逻辑移植部分 (简化 BAM 收集) ---
        
        
        
        # 1. 执行samtools merge
        print("\n# 1. 执行samtools merge")
        if os.path.exists(GLOBAL_MERGED_BAM):
            print(f"  目标文件 {GLOBAL_MERGED_BAM} 已存在，跳过samtools merge")
        else:
            print("  执行samtools merge合并所有BAM文件...")

            # 1. 查找所有要合并的 BAM 文件 (使用 glob 简化逻辑)
            print(f"\n  正在 {BAM_SOURCE_DIR} 目录下查找所有 *sorted.bam 文件...")
            
            # 使用 glob 查找所有匹配的文件
            search_pattern = os.path.join(BAM_SOURCE_DIR, '*sorted.bam')
            ALL_SORTED_BAMS: List[str] = glob.glob(search_pattern)

            if not ALL_SORTED_BAMS:
                print(f"错误: 在目录 {BAM_SOURCE_DIR} 中没有找到任何 *sorted.bam 文件可供合并")
                sys.exit(1)
            
            print(f"  找到 {len(ALL_SORTED_BAMS)} 个 BAM 文件准备合并。")
            
            # 使用 subprocess.run 替换 os.system，参数列表包含所有 BAM 文件路径
            command_list = [
                SAMTOOLS_BIN, 'merge', '-@', '8', '-c', '-f', GLOBAL_MERGED_BAM
            ] + ALL_SORTED_BAMS
            
            try:
                # check=True: 如果返回码非零，则抛出异常
                subprocess.run(
                    command_list,
                    check=True,
                    capture_output=True, # 捕获 stdout 和 stderr
                    text=True
                )
                print(f"  所有BAM文件已合并到 {GLOBAL_MERGED_BAM}")
            except subprocess.CalledProcessError as e:
                print(f"错误: samtools merge失败。命令: {' '.join(command_list)}")
                print(f"stderr:\n{e.stderr}")
                sys.exit(1)
        
        # 2. 计算合并后BAM文件的深度
        print("\n# 2. 计算合并后BAM文件的深度")
        if os.path.exists(GLOBAL_MERGED_DEPTH):
            print(f"  目标文件 {GLOBAL_MERGED_DEPTH} 已存在，跳过samtools depth")
        else:
            print("  计算合并后BAM文件的深度...")
            
            # samtools depth "${GLOBAL_MERGED_BAM}" > "${GLOBAL_MERGED_DEPTH}"
            command_list = [SAMTOOLS_BIN, 'depth', GLOBAL_MERGED_BAM]
            
            try:
                # 将标准输出重定向到文件
                with open(GLOBAL_MERGED_DEPTH, 'w') as f:
                    subprocess.run(
                        command_list,
                        check=True,
                        stdout=f, # 将标准输出写入文件对象
                        stderr=subprocess.PIPE,
                        text=True
                    )
                print(f"  深度已计算并保存到 {GLOBAL_MERGED_DEPTH}")
            except subprocess.CalledProcessError as e:
                print(f"错误: samtools depth失败。命令: {' '.join(command_list)}")
                print(f"stderr:\n{e.stderr}")
                sys.exit(1)

        # 3a. Calculate depth average
        # 使用新生成的深度文件路径
        depth_path = GLOBAL_MERGED_DEPTH
        bam_for_graph = GLOBAL_MERGED_BAM
        AVERAGE_CACHE_FILE = os.path.join(args.output_path, 'merged.depth.num')
        average = None
        if os.path.exists(AVERAGE_CACHE_FILE):
            # 如果文件存在，尝试读取
            with open(AVERAGE_CACHE_FILE, 'r') as f:
                # 读取并尝试转换为浮点数
                cached_average_str = f.read().strip()
                average = float(cached_average_str)
        if average is None:
            print(f"⚠️ 缓存未找到或无效。开始执行计算...")
            average = calculate_average_third_column(depth_path)
            with open(AVERAGE_CACHE_FILE, 'w') as f:
                f.write(str(average))        
        # 3b. Set library path and generate G13 graph
        print("\n# 3b. Set library path and generate G13 graph")
        
        # 使用 subprocess.run 结合 env 参数来设置 LD_LIBRARY_PATH，比 os.system 更安全
        env = os.environ.copy()
        env['LD_LIBRARY_PATH'] = args.LD_LIBRARY_PATH

        command_args = [
            os.path.join(args.linking_path, 'save_models', 'generateG13'),
            bam_for_graph, # 使用合并后的 BAM 文件
            GLOBAL_MERGED_GRAPH, # 目标图文件路径
            str(average),
            '-e', '300',
            '--min-count', '1'
        ]
        
        try:
            subprocess.run(
                command_args,
                check=True,
                env=env,
                capture_output=True,
                text=True
            )
            print("  generateG13 完成。")
        except subprocess.CalledProcessError as e:
            print(f"错误: generateG13 失败。命令: {' '.join(command_args)}")
            print(f"stderr:\n{e.stderr}")
            sys.exit(1)

    # --- 后续流程（使用新的全局路径）---
    
    # 3c. generate bin
    bin_path = os.path.join(args.output_path, 'binning', 'bins')
    link_bin_path = os.path.join(args.output_path, 'linking', 'bins_0.9')
    report_path = os.path.join(args.output_path, 'binning', 'checkm2_bins', 'quality_report.tsv')
    filter_and_copy_bins(report_path, bin_path, link_bin_path)

    # 3d. merge fasta
    bining_dir = os.path.join(args.output_path, 'linking')
    fasta_path = os.path.join(bining_dir, 'all_bins_merged.fasta')
    merge_fasta_files(bining_dir, fasta_path)

    #3e. run prodigal and clean faa
    prodigal_log_path = os.path.join(bining_dir, 'prodigal.log') # 定义日志文件路径
    pprodigal.main(args, fasta_path, prodigal_log_path) # 传入日志文件路径
    # pprodigal.main(args, fasta_path)
    out_path = os.path.join(bining_dir, 'binning_clean.faa')
    in_path = os.path.join(bining_dir, 'bining.faa')
    remove_asterisks_from_file(in_path, out_path)

    bin_index_path = os.path.join(bining_dir, 'binning_index.txt')
    process_faa(out_path, bin_index_path)

    # 3f. generate input pkl
    bam_graph_path = GLOBAL_MERGED_GRAPH
    final_input_pkl_path = os.path.join(bining_dir, 'final_input.pkl')
    generate_inputpkl_main(args, bam_graph_path, bin_index_path, out_path, final_input_pkl_path)

    # 3g. generate cosine bin
    cos_path = os.path.join(bining_dir, 'cosine')
    model_path = require_existing_path(
        os.path.join(args.linking_path, 'save_models', 'bacteria_transformer2.pth'),
    )
    generate_cosine_bin1.main(final_input_pkl_path, bam_graph_path, model_path, cos_path)

    return

if __name__ == "__main__":
    main()
