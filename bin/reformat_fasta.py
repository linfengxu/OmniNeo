#!/usr/bin/env python3

import sys
import re

def reformat_fasta(input_file, output_file):
    """
    重新格式化FASTA文件，将标题从 ">NAME | NAME" 转换为 ">sp|NAME|NAME"
    """
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith('>'):
                # 提取标题中的名称
                match = re.match(r'^>(\S+)(?:\s+\|\s+(\S+))?', line)
                if match:
                    name1 = match.group(1)
                    name2 = match.group(2) if match.group(2) else name1
                    # 重新格式化为 UniProt 格式
                    new_header = f">sp|{name1}|{name2}\n"
                    outfile.write(new_header)
                else:
                    # 如果无法解析，保持原样
                    outfile.write(line)
            else:
                outfile.write(line)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"用法: {sys.argv[0]} 输入文件.fasta 输出文件.fasta")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    reformat_fasta(input_file, output_file)
    print(f"FASTA文件已重新格式化并保存为 {output_file}") 