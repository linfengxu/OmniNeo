#!/usr/bin/env python
import pandas as pd
import re
import argparse
import os
from Bio.Seq import Seq

class FusionPeptideGenerator:
    def __init__(self, input_file, output_dir, short_peptide_len=(8,11), long_peptide_len=(12,25)):
        """初始化融合肽段生成器"""
        self.input_file = input_file
        self.output_dir = output_dir
        self.short_range = short_peptide_len
        self.long_range = long_peptide_len
        
        # 创建输出目录
        os.makedirs(output_dir, exist_ok=True)
        os.makedirs(os.path.join(output_dir, "fasta"), exist_ok=True)
        os.makedirs(os.path.join(output_dir, "csv"), exist_ok=True)
        
        # 读取融合预测结果
        self.fusions = pd.read_csv(input_file, sep='\t')
        
    def process_fusions(self):
        """处理所有具有蛋白质序列的融合事件"""
        fusion_peptides = {}
        
        # 筛选有蛋白质翻译和有效CDS范围的融合
        valid_fusions = self.fusions[
            (self.fusions['FUSION_TRANSL'].notna()) & 
            (self.fusions['CDS_LEFT_RANGE'].notna()) & 
            (self.fusions['CDS_RIGHT_RANGE'].notna()) &
            (self.fusions['CDS_LEFT_RANGE'] != '.') &
            (self.fusions['CDS_RIGHT_RANGE'] != '.')
        ]
        
        print(f"Processing {len(valid_fusions)} fusions with protein translations...")
        
        for idx, fusion in valid_fusions.iterrows():
            fusion_name = fusion['#FusionName']
            print(f"Processing fusion: {fusion_name}")
                
            # 解析CDS范围
            try:
                left_start, left_end = map(int, fusion['CDS_LEFT_RANGE'].split('-'))
                right_start, right_end = map(int, fusion['CDS_RIGHT_RANGE'].split('-'))
            except (ValueError, AttributeError) as e:
                print(f"  Skipping {fusion_name}: invalid CDS range - {e}")
                continue
            
            # 计算对应的氨基酸位置
            left_aa_len = (left_end - left_start + 1) // 3
            fusion_point = left_aa_len
            
            # 检查是否有移码
            is_frameshift = fusion['PROT_FUSION_TYPE'] == 'FRAMESHIFT'
            frame_offset = (left_end - left_start + 1) % 3
            print(f"  Fusion point at AA position: {fusion_point}")
            print(f"  Frameshift: {is_frameshift}, Frame offset: {frame_offset}")
            
            # 获取融合蛋白序列并处理*终止密码子
            fusion_protein = fusion['FUSION_TRANSL']
            if '*' in fusion_protein:
                # 保留第一个*前的序列
                fusion_protein = fusion_protein.split('*')[0]
            
            # 生成短肽和长肽
            short_peptides = self.generate_junction_peptides(
                fusion_protein, fusion_point, self.short_range[0], self.short_range[1])
            
            long_peptides = self.generate_junction_peptides(
                fusion_protein, fusion_point, self.long_range[0], self.long_range[1])
            
            # 存储结果
            fusion_peptides[fusion_name] = {
                'short_peptides': short_peptides,
                'long_peptides': long_peptides,
                'fusion_point': fusion_point,
                'is_frameshift': is_frameshift
            }
            
            # 保存为文件
            self.save_peptides(fusion_name, short_peptides, 'short')
            self.save_peptides(fusion_name, long_peptides, 'long')
            
        return fusion_peptides
    
    def generate_junction_peptides(self, fusion_protein, fusion_point, min_len, max_len):
        """生成跨越融合点的肽段"""
        junction_peptides = []
        
        # 确保前后至少有1个氨基酸来自各自的基因
        for peptide_len in range(min_len, max_len + 1):
            for aa_before in range(1, peptide_len):
                aa_after = peptide_len - aa_before
                
                # 检查是否有足够的氨基酸
                if fusion_point - aa_before < 0 or fusion_point + aa_after > len(fusion_protein):
                    continue
                
                # 提取肽段
                start_pos = fusion_point - aa_before
                end_pos = fusion_point + aa_after
                peptide = fusion_protein[start_pos:end_pos]
                
                # 肽段元数据
                peptide_info = {
                    'peptide': peptide,
                    'start': start_pos,
                    'end': end_pos,
                    'length': peptide_len,
                    'aa_before': aa_before,
                    'aa_after': aa_after
                }
                
                junction_peptides.append(peptide_info)
        
        return junction_peptides
    
    def save_peptides(self, fusion_name, peptides, peptide_type):
        """保存肽段到FASTA和CSV文件"""
        if not peptides:
            return
            
        # 保存为FASTA
        fasta_path = os.path.join(self.output_dir, "fasta", f"{fusion_name}_{peptide_type}_peptides.fasta")
        with open(fasta_path, 'w') as f:
            for i, pep in enumerate(peptides):
                header = f">{fusion_name}|junction|pos={pep['start']+1}-{pep['end']}|len={pep['length']}|before={pep['aa_before']}|after={pep['aa_after']}"
                f.write(f"{header}\n{pep['peptide']}\n")
        
        # 保存为CSV
        csv_path = os.path.join(self.output_dir, "csv", f"{fusion_name}_{peptide_type}_peptides.csv")
        df = pd.DataFrame(peptides)
        df['fusion_name'] = fusion_name
        df.to_csv(csv_path, index=False)
        
        print(f"  Saved {len(peptides)} {peptide_type} peptides to {fasta_path} and {csv_path}")

def main():
    parser = argparse.ArgumentParser(description='Generate junction peptides from STAR-Fusion results')
    parser.add_argument('-i', '--input', required=True, help='STAR-Fusion coding effect results file')
    parser.add_argument('-o', '--output', required=True, help='Output directory for peptide files')
    parser.add_argument('--short_min', type=int, default=8, help='Minimum length for short peptides')
    parser.add_argument('--short_max', type=int, default=11, help='Maximum length for short peptides')
    parser.add_argument('--long_min', type=int, default=12, help='Minimum length for long peptides')
    parser.add_argument('--long_max', type=int, default=25, help='Maximum length for long peptides')
    
    args = parser.parse_args()
    
    generator = FusionPeptideGenerator(
        args.input, 
        args.output,
        short_peptide_len=(args.short_min, args.short_max),
        long_peptide_len=(args.long_min, args.long_max)
    )
    
    generator.process_fusions()

if __name__ == "__main__":
    main()