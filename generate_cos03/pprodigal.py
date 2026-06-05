import argparse
import os
import subprocess
import tempfile
import shutil
import sys
import re

def get_fasta_records(input_file):
    """
    获取FASTA文件中每个记录的字节范围
    
    Args:
        input_file (str): 输入文件路径
    
    Returns:
        list: 包含每个FASTA记录的起始和结束字节位置的元组列表
    """
    records = []
    with open(input_file, 'rb') as f:
        start = None
        while True:
            pos = f.tell()
            line = f.readline()
            if not line:
                break
            if line.startswith(b'>'):
                if start is not None:
                    records.append((start, pos))
                start = pos
        # 处理最后一个记录
        if start is not None:
            f.seek(0, 2)
            end = f.tell()
            records.append((start, end))
    return records

def split_input(input_file, num_parts, temp_dir):
    """
    将输入文件分割成多个部分

    Args:
        input_file (str): 输入文件路径
        num_parts (int): 分割的份数
        temp_dir (str): 临时目录路径

    Returns:
        list: 分割后的子文件路径列表
    """
    records = get_fasta_records(input_file)
    total = len(records)
    num_parts = min(num_parts, total)  # 不能超过记录总数
    
    # 计算每个分块包含的记录数
    chunk_size, remainder = divmod(total, num_parts)
    chunks = []
    current = 0
    
    for i in range(num_parts):
        count = chunk_size + (1 if i < remainder else 0)
        if current >= total:
            break
        end = current + count
        chunk = records[current:end]
        chunks.append((chunk[0][0], chunk[-1][1]))
        current = end
    
    # 写入分块文件
    sub_files = []
    for i, (start, end) in enumerate(chunks):
        sub_file = os.path.join(temp_dir, f'part_{i}.fasta')
        with open(input_file, 'rb') as f_in:
            f_in.seek(start)
            data = f_in.read(end - start)
        with open(sub_file, 'wb') as f_out:
            f_out.write(data)
        sub_files.append(sub_file)
    return sub_files

def run_prodigal(sub_files, temp_dir, prodigal_args, meta, output_faa, output_fna, output_gff, log_file_path=None):
    """
    并行运行Prodigal，并将Prodigal的日志输出重定向到指定文件。

    Args:
        sub_files (list): 分割后的子文件列表
        temp_dir (str): 临时目录路径
        prodigal_args (str): 额外的prodigal参数
        meta (str): prodigal模式 ('single' 或 'meta')
        output_faa (bool): 是否输出蛋白序列
        output_fna (bool): 是否输出基因序列
        output_gff (bool): 是否输出GFF文件
        log_file_path (str, optional): Prodigal运行日志的路径。如果提供，
                                       subprocess的stdout/stderr将被追加重定向到此文件。

    Returns:
        tuple: (faa_files, fna_files, gff_files)
    """
    processes = []
    faa_files = []
    fna_files = []
    gff_files = []
    
    # 准备日志文件句柄和重定向目标
    log_file_handle = None
    subprocess_out = None # 默认不重定向 (即输出到终端)

    if log_file_path:
        try:
            # 使用 'a' 模式追加日志
            log_file_handle = open(log_file_path, 'a')
            subprocess_out = log_file_handle
        except IOError as e:
            print(f"Warning: Could not open log file {log_file_path} for writing: {e}", file=sys.stderr)
            # 警告后，继续使用默认输出
            
    try:
        for i, sub_file in enumerate(sub_files):
            cmd = ['prodigal', '-i', sub_file]
            
            # 添加输出参数
            if output_faa:
                faa_output = os.path.join(temp_dir, f'output_{i}.faa')
                cmd.extend(['-a', faa_output])
                faa_files.append(faa_output)
            
            if output_fna:
                fna_output = os.path.join(temp_dir, f'output_{i}.fna')
                cmd.extend(['-d', fna_output])
                fna_files.append(fna_output)
                
            if output_gff:
                gff_output = os.path.join(temp_dir, f'output_{i}.gff')
                cmd.extend(['-o', gff_output, '-f', 'gff'])
                gff_files.append(gff_output)

            # 检查是否至少指定了一种输出格式
            if not (output_faa or output_fna or output_gff):
                error_msg = "Error: No output format specified. Please use at least one of -a, -d, or -g.\n"
                if log_file_handle:
                    log_file_handle.write(error_msg)
                print(error_msg, file=sys.stderr)
                return [], [], []

            if meta == 'meta':
                cmd.extend(['-p', 'meta'])
            if prodigal_args:
                cmd.extend(prodigal_args.split())
                
            command_str = f"Executing Prodigal part {i}: {' '.join(cmd)}\n"
            
            # 将执行命令写入日志文件（如果日志文件已打开）
            if log_file_handle:
                log_file_handle.write(command_str)

            # 保持在交互窗中打印执行命令，便于调试
            print(command_str.strip()) 
            
            # 启动子进程，并将 stdout 和 stderr 重定向到日志文件句柄 (如果存在)
            # Prodigal 通常将进度和非错误信息输出到 stderr
            processes.append(subprocess.Popen(
                cmd, 
                stdout=subprocess_out, 
                stderr=subprocess_out,
                universal_newlines=True # 使用文本模式处理输出
            ))
        
        # 等待所有进程完成
        for p in processes:
            p.wait()
            
    finally:
        # 确保关闭日志文件句柄
        if log_file_handle:
            log_file_handle.close()
    
    return faa_files, fna_files, gff_files

def merge_fasta_outputs(output_files, final_output):
    """
    合并FASTA输出文件
    """
    with open(final_output, 'wb') as f_out:
        for f in output_files:
            with open(f, 'rb') as f_in:
                f_out.write(f_in.read())
            
def merge_gff_outputs(output_files, final_output):
    """
    合并GFF输出文件，移除重复的头部信息
    """
    with open(final_output, 'w') as f_out:
        header_written = False
        for f in output_files:
            with open(f, 'r') as f_in:
                lines = f_in.readlines()
                for line in lines:
                    if line.startswith('##gff-version') or line.startswith('##sequence-region'):
                        if not header_written:
                            f_out.write(line)
                        continue
                    if line.startswith('##FASTA'):
                        header_written = True
                        continue
                    f_out.write(line)
                    
def main(args, fasta_path, prodigal_log_path=None):
    """
    主函数：协调输入分割、Prodigal并行运行和结果合并。
    
    Args:
        args (object): 包含命令行参数的对象。
        fasta_path (str): 输入的FASTA文件路径。
        prodigal_log_path (str, optional): Prodigal运行日志的输出路径。
    """
    # 假设 args 中包含 output_path, gff_output, num_parts, prodigal_args, meta, output_faa, output_fna
    
    # 保持原有逻辑，设置输出路径
    args.prodigal_output = args.output_path+"/linking/bining.faa"
    args.prodigal_input = fasta_path
    
    if not (args.prodigal_output or args.gff_output):
        print("At least one output file must be specified using -o or -gff_output.", file=sys.stderr)
        return # 使用 return 退出函数

    temp_dir = tempfile.mkdtemp()
    
    # 打印运行状态信息
    log_message = f"Using temporary directory: {temp_dir}\n"
    if prodigal_log_path:
        log_message += f"Prodigal logs will be appended to: {prodigal_log_path}\n"
        # 尝试将主程序日志也写入日志文件
        try:
             with open(prodigal_log_path, 'a') as f:
                f.write(log_message)
        except Exception:
            # 写入失败不影响程序运行，只在控制台打印
            pass
            
    print(log_message.strip()) # 在控制台打印状态信息

    try:
        sub_files = split_input(args.prodigal_input, args.num_parts, temp_dir)
        
        # 运行Prodigal，传入日志文件路径
        faa_files, fna_files, gff_files = run_prodigal(
            sub_files, temp_dir, args.prodigal_args, args.meta, 
            args.output_faa, args.output_fna, bool(args.gff_output),
            prodigal_log_path # 传入日志路径
        )
        
        # 合并FASTA文件
        if args.prodigal_output:
            if args.output_faa and args.output_fna:
                print("Error: Cannot specify both -a and -d for a single output file. Please choose one.", file=sys.stderr)
                sys.exit(1)
            elif args.output_faa:
                merge_fasta_outputs(faa_files, args.prodigal_output)
                print(f"Processing complete. Final FAA output is at: {args.prodigal_output}")
            elif args.output_fna:
                merge_fasta_outputs(fna_files, args.prodigal_output)
                print(f"Processing complete. Final FNA output is at: {args.prodigal_output}")
            else:
                print("Error: -o was specified but no output type (-a or -d) was provided.", file=sys.stderr)
                sys.exit(1)

        # 合并GFF文件
        if args.gff_output:
            merge_gff_outputs(gff_files, args.gff_output)
            print(f"Processing complete. Final GFF output is at: {args.gff_output}")

    except Exception as e:
        print(f"An error occurred: {e}", file=sys.stderr)
    finally:
        print(f"Cleaning up temporary directory: {temp_dir}")
        shutil.rmtree(temp_dir)

if __name__ == '__main__':
    main()