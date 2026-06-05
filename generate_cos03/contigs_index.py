import re
import sys

def process_faa(faa_file, output_file):
    contig_genes = {}
    
    with open(faa_file, 'r') as f:
        current_contig = None
        for line in f:
            if line.startswith('>'):
                # 解析头部信息
                header = line.strip().split()[0]  # 取第一个字段
                contig_gene = header[1:]  # 去掉开头的>
                
                # 分割contig名称和基因编号
                if '_' in contig_gene:
                    contig_id, gene_num = contig_gene.rsplit('_', 1)
                else:
                    contig_id = contig_gene
                    gene_num = "1"
                
                # 解析方向信息
                direction_match = re.search(r'# [^#]+ # [^#]+ # (-?\d+) #', line)
                direction = '-' if direction_match and int(direction_match.group(1)) < 0 else '+'
                
                # 生成基因标识符
                gene_id = f"{direction}{contig_gene}"
                
                # 添加到字典
                if contig_id not in contig_genes:
                    contig_genes[contig_id] = []
                contig_genes[contig_id].append(gene_id)
    
    # 写入输出文件
    with open(output_file, 'w') as f_out:
        for contig_id, genes in contig_genes.items():
            # 按基因编号排序（假设编号代表位置顺序）
            try:
                genes.sort(key=lambda x: int(x.split('_')[-1]))
            except:
                pass
            f_out.write(f"{contig_id}\t{';'.join(genes)}\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <input.faa> <output_index.txt>")
        sys.exit(1)
    
    faa_file = sys.argv[1]
    output_file = sys.argv[2]
    
    process_faa(faa_file, output_file)


# python script/0.4.contigs_index.py ~/project/contigs/gcn/CAMI/CAMI1/low/RL_S001__insert_270/RL_S001__insert_270_1000.faa process/cami1_low/contigs_index.txt