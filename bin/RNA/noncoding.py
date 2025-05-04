#!/usr/bin/env python
# -*- coding: utf-8 -*-i


"""
Generation of mutated peptides:
"""

from os import remove
from functools import reduce
import os
import pandas as pd
import sys
import argparse
from Bio import SeqIO
from Bio.Seq import Seq


def rev_comp(s):
    sc = ''
    for i in reversed(s):
        if i == 'A':
            sc += 'T'
        elif i == 'C':
            sc += 'G'
        elif i == 'G':
            sc += 'C'
        elif i == 'T':
            sc += 'A'
    return sc


def get_mutation(f, p_df):
    train_pos_data_path = f 
    alt_value = p_df 
    fh = open(train_pos_data_path, 'r')
    seq = {}
    for line in fh:
        if line.startswith('>'):
            name = line.replace('>', '').split()[0]
            seq[name] = ''
        else:
            seq[name] += line.replace('\n', '')
    fh.close()
    wild_list = []
    for key in seq:
        wild_list.append(seq[key])
    mut_list = list()
    mut_pos_list = alt_value['ALT'].tolist()
    for i in range(0, len(wild_list)):
        re_front = wild_list[i][0:100]
        re_back = wild_list[i][-100:]
        re = re_front + mut_pos_list[i] + re_back
        re = re.upper()
        mut_list.append(re)
    out_dict = {"ID": alt_value["ID"].tolist(), "POS": alt_value["POS"], "X.CHROM": alt_value["X.CHROM"].tolist(),
                "ALT": alt_value["ALT"].tolist(), "REF": alt_value["REF"].tolist(), "Gene": alt_value["Gene"].tolist(),
                "wild": wild_list, "mut": mut_list
                }
    out_df = pd.DataFrame(out_dict)
    return out_df


def sixFrame_translate_fun(strs):
    lst1 = strs.strip("\n")
    lst1 = lst1.upper()
    lst1 = lst1.replace('N', 'A')
    trans_cod1 = "".join([codon_table["".join(lst1[i:i + 3])] for i in range(0, len(lst1) - len(lst1) % 3, 3)])
    trans_cod2 = "".join(
        [codon_table["".join(lst1[i:i + 3])] for i in range(1, (len(lst1) - 1) - (len(lst1) - 1) % 3, 3)])
    trans_cod3 = "".join(
        [codon_table["".join(lst1[i:i + 3])] for i in range(2, (len(lst1) - 2) - (len(lst1) - 2) % 3, 3)])
    return ([trans_cod1, trans_cod2, trans_cod3])


def align_translate_func(srs, order_trans="pl"):
    if order_trans == "pl":
        mut_base_str = srs["mut"]
        wild_base_str = srs["wild"]
    elif order_trans == "nl":
        mut_base_str = rev_comp(srs["mut"])
        wild_base_str = rev_comp(srs["wild"])
    else:
        exit("Either pl or nl is permitted in align_translate_func")

    tr1_list = sixFrame_translate_fun(mut_base_str)
    tr2_list = sixFrame_translate_fun(wild_base_str)
    tmp_list = list(map(lambda x, y: x if (x != y) else "", tr1_list, tr2_list))
    tmp_seq = ["trans" + str(i) for i in range(1, len(tr1_list) + 1)]
    out_dict = {'ID': srs["ID"] + "|" + order_trans, "POS": srs["POS"], "X.CHROM": srs["X.CHROM"],
                "ALT": srs["ALT"], "REF": srs["REF"], "Gene": srs['Gene'],
                'mut_pro': tmp_list, 'mut_base': mut_base_str,
                'wild_pro': tr2_list, 'wild_base': wild_base_str,
                "trans": tmp_seq}
    out_df = pd.DataFrame(out_dict)
    return (out_df)


def mk_fastaStr(pdSeries):
    fastaStr_lst = []
    for strs in ["prot3", "prot2", "prot1"]:
        ant_prot = strs + "-" + 'ant'
        mut_prot = strs + "-" + 'mut'
        mut_prot_pos = mut_prot + "-" + "pos"
        if (pd.isna(pdSeries[ant_prot]) or len(pdSeries[ant_prot]) < 7):
            continue
        else:
            try:
                global header1
                header1 = pdSeries["ID"] + "|" + pdSeries["gene"] + "|" + pdSeries[mut_prot] + "|" + pdSeries[
                    mut_prot_pos] + "\n"
            except:
                print("Mutated protein %s translaste error", pdSeries["ID"])
            header2 = header1.replace(";", "|")
            fastaStr = header2 + pdSeries[ant_prot]
        fastaStr_lst.append(fastaStr)
    return (fastaStr_lst)


class NoncodingPeptideGenerator:
    def __init__(self, anno_file, ref_genome, output_dir, 
                 short_peptide_len=(8,11), long_peptide_len=(12,25)):
        """初始化非编码区突变肽段生成器"""
        self.anno_file = anno_file
        self.ref_genome = ref_genome
        self.output_dir = output_dir
        self.short_range = short_peptide_len
        self.long_range = long_peptide_len
        
        # 创建输出目录
        os.makedirs(output_dir, exist_ok=True)
        os.makedirs(os.path.join(output_dir, "fasta"), exist_ok=True)
        os.makedirs(os.path.join(output_dir, "csv"), exist_ok=True)
        
        # 密码子表
        self.codon_table = {
            "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
            "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
            "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
            "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
            "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
            "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
            "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
            "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
            "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
            "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
            "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
            "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
            "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
            "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
            "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
            "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"
        }
    
    def rev_comp(self, seq):
        """生成反向互补序列，处理非标准碱基字符"""
        complement = {
            'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N',
            'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n',
            # 处理特殊字符
            '-': '-', '.': '.', '*': '*', '?': '?',
            # IUPAC模糊碱基代码
            'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W',
            'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B',
            'D': 'H', 'H': 'D',
            'r': 'y', 'y': 'r', 's': 's', 'w': 'w',
            'k': 'm', 'm': 'k', 'b': 'v', 'v': 'b',
            'd': 'h', 'h': 'd'
        }
        
        result = []
        seq = seq.upper()
        for base in reversed(seq):
            if base in complement:
                result.append(complement[base])
            else:
                # 对未知字符，保持原样
                result.append(base)
        
        return ''.join(result)
    
    def sixFrame_translate(self, seq):
        """六框翻译序列"""
        seq = seq.upper().replace('N', 'A')  # 替换N为A以便翻译
        translations = []
        
        # 正向3个框架
        for i in range(3):
            prot = ""
            for j in range(i, len(seq)-2, 3):
                codon = seq[j:j+3]
                if len(codon) == 3:
                    prot += self.codon_table.get(codon, "X")
            translations.append(prot)
        
        # 反向3个框架
        rev_seq = self.rev_comp(seq)
        for i in range(3):
            prot = ""
            for j in range(i, len(rev_seq)-2, 3):
                codon = rev_seq[j:j+3]
                if len(codon) == 3:
                    prot += self.codon_table.get(codon, "X")
            translations.append(prot)
        
        return translations
    
    def extract_noncoding_mutations(self):
        """从ANNOVAR结果中提取非编码区突变"""
        print(f"Reading ANNOVAR annotation from {self.anno_file}")
        anno_df = pd.read_csv(self.anno_file, sep="\t")
        
        # 筛选非编码区突变
        noncoding_df = anno_df[~anno_df['Func.refGene'].isin(('exonic', 'exonic;splicing'))]
        print(f"Found {len(noncoding_df)} noncoding mutations")
        
        # 重命名列并添加需要的字段
        if 'Otherinfo' in noncoding_df.columns:
            # 扩展otherinfo列（如果采用了不同的ANNOVAR输出格式）
            info_cols = ['chrom', 'pos', 'id', 'ref', 'alt']
            for i, col in enumerate(info_cols):
                colname = f"Otherinfo{i+1}"
                if colname in noncoding_df.columns:
                    noncoding_df[col] = noncoding_df[colname]
        else:
            # 使用标准列名
            noncoding_df['chrom'] = noncoding_df['Chr']
            noncoding_df['pos'] = noncoding_df['Start']
            noncoding_df['ref'] = noncoding_df['Ref']
            noncoding_df['alt'] = noncoding_df['Alt']
        
        # 生成ID和区域范围
        noncoding_df['id'] = '>' + noncoding_df['chrom'] + '|' + noncoding_df['pos'].astype(str) + '|' + noncoding_df['ref'] + '|' + noncoding_df['alt']
        noncoding_df['begin'] = noncoding_df['pos'] - 101  # 提取前100bp
        noncoding_df['end'] = noncoding_df['pos'] + noncoding_df['ref'].str.len() + 100  # 提取后100bp
        
        # 创建BED文件
        bed_path = os.path.join(self.output_dir, "noncoding_regions.bed")
        noncoding_df[['chrom', 'begin', 'end']].to_csv(
            bed_path, sep="\t", header=False, index=False)
        
        return noncoding_df
    
    def extract_sequences(self, noncoding_df):
        """从参考基因组提取区域序列，使用纯Python实现"""
        print(f"Reading reference genome: {self.ref_genome}")
        fasta_path = os.path.join(self.output_dir, "wild_sequences.fasta")
        
        # 尝试使用pyfaidx (如果安装了)
        try:
            import pyfaidx
            print("Using pyfaidx for efficient sequence extraction")
            ref = pyfaidx.Fasta(self.ref_genome)
            
            # 提取序列
            sequences = {}
            with open(fasta_path, 'w') as out_file:
                for idx, row in noncoding_df.iterrows():
                    try:
                        chrom = row['chrom']
                        start = int(row['begin'])
                        end = int(row['end'])
                        region_id = f"{chrom}:{start}-{end}"
                        
                        # 确保坐标有效
                        if start < 0:
                            start = 0
                        
                        # 使用pyfaidx提取序列
                        if chrom in ref:
                            seq = ref[chrom][start:end].seq
                            # 保存到FASTA文件
                            out_file.write(f">{region_id}\n{seq}\n")
                            sequences[region_id] = seq
                        else:
                            print(f"Warning: Chromosome {chrom} not found in reference")
                    except Exception as e:
                        print(f"Error extracting region {chrom}:{start}-{end}: {e}")
            
            return sequences, fasta_path
        
        except ImportError:
            print("pyfaidx not found, using standard BioPython (slower but works)")
        
        # 回退到基本的BioPython实现
        sequences = {}
        
        # 直接从FASTA读取全部序列(慢但可靠)
        ref_dict = {}
        for record in SeqIO.parse(self.ref_genome, "fasta"):
            ref_dict[record.id] = record
        
        # 提取序列
        with open(fasta_path, 'w') as out_file:
            for idx, row in noncoding_df.iterrows():
                chrom = row['chrom']
                start = max(0, int(row['begin']))
                end = int(row['end'])
                region_id = f"{chrom}:{start}-{end}"
                
                if chrom in ref_dict:
                    seq_record = ref_dict[chrom]
                    # 确保坐标不超出序列范围
                    if end > len(seq_record.seq):
                        end = len(seq_record.seq)
                    seq = str(seq_record.seq[start:end])
                    
                    # 保存到FASTA文件
                    out_file.write(f">{region_id}\n{seq}\n")
                    sequences[region_id] = seq
                else:
                    print(f"Warning: Chromosome {chrom} not found in reference")
        
        print(f"Extracted {len(sequences)} sequences and saved to {fasta_path}")
        return sequences, fasta_path
    
    def generate_mutant_sequences(self, noncoding_df, wild_sequences):
        """生成突变序列"""
        print("Generating mutant sequences")
        mutant_data = []
        
        mismatches = 0
        for idx, row in noncoding_df.iterrows():
            chrom = row['chrom']
            pos = row['pos']
            ref = row['ref']
            alt = row['alt']
            region_id = f"{chrom}:{row['begin']}-{row['end']}"
            
            if region_id in wild_sequences:
                wild_seq = wild_sequences[region_id]
                # 计算突变位置在提取序列中的位置（0-based）
                rel_pos = pos - row['begin']
                
                # 生成突变序列
                if 0 <= rel_pos < len(wild_seq):
                    # 确保ref匹配
                    if wild_seq[rel_pos:rel_pos+len(ref)].upper() == ref.upper():
                        mut_seq = wild_seq[:rel_pos] + alt + wild_seq[rel_pos+len(ref):]
                        
                        mutant_data.append({
                            'id': row['id'],
                            'chrom': chrom,
                            'pos': pos,
                            'ref': ref,
                            'alt': alt,
                            'gene': row.get('Gene.refGene', '.'),
                            'wild_seq': wild_seq,
                            'mut_seq': mut_seq,
                            'rel_pos': rel_pos
                        })
                    else:
                        mismatches += 1
                        print(f"Warning: Reference mismatch at {chrom}:{pos} - Expected {ref}, found {wild_seq[rel_pos:rel_pos+len(ref)]}")
                        
                        # 尝试附近位置查找匹配
                        found_match = False
                        for offset in [-1, 1]:
                            check_pos = rel_pos + offset
                            if 0 <= check_pos < len(wild_seq) - len(ref) + 1:
                                if wild_seq[check_pos:check_pos+len(ref)].upper() == ref.upper():
                                    print(f"  Found match with offset {offset}")
                                    rel_pos = check_pos
                                    mut_seq = wild_seq[:rel_pos] + alt + wild_seq[rel_pos+len(ref):]
                                    
                                    mutant_data.append({
                                        'id': row['id'],
                                        'chrom': chrom,
                                        'pos': pos,
                                        'ref': ref,
                                        'alt': alt,
                                        'gene': row.get('Gene.refGene', '.'),
                                        'wild_seq': wild_seq,
                                        'mut_seq': mut_seq,
                                        'rel_pos': rel_pos
                                    })
                                    found_match = True
                                    break
                        
                        if not found_match:
                            print(f"  No match found in nearby positions")
        
        print(f"Generated {len(mutant_data)} mutant sequences ({mismatches} reference mismatches)")
        return pd.DataFrame(mutant_data)
    
    def translate_and_compare(self, mut_df):
        """翻译并比较野生型和突变序列"""
        print(f"Translating and comparing sequences for {len(mut_df)} mutations")
        translation_results = []
        
        # 保存所有翻译结果用于质谱分析
        self.all_translations = []
        
        for idx, row in mut_df.iterrows():
            # 正向翻译
            wild_translations = self.sixFrame_translate(row['wild_seq'])
            mut_translations = self.sixFrame_translate(row['mut_seq'])
            
            # 保存所有突变序列的翻译结果（用于质谱分析）
            for i in range(len(mut_translations)):
                frame = i % 3 + 1
                strand = "+" if i < 3 else "-"
                
                self.all_translations.append({
                    'id': row['id'],
                    'chrom': row['chrom'],
                    'pos': row['pos'],
                    'ref': row['ref'],
                    'alt': row['alt'],
                    'gene': row['gene'],
                    'frame': f"frame{frame}",
                    'strand': strand,
                    'prot_seq': mut_translations[i],
                    'is_mutated': wild_translations[i] != mut_translations[i]
                })
            
            # 比较翻译结果
            mutation_found = False
            for i in range(len(wild_translations)):
                frame = i % 3 + 1
                strand = "+" if i < 3 else "-"
                
                wild_prot = wild_translations[i]
                mut_prot = mut_translations[i]
                
                # 找出差异
                if wild_prot != mut_prot:
                    mutation_found = True
                    translation_results.append({
                        'id': row['id'],
                        'chrom': row['chrom'],
                        'pos': row['pos'],
                        'ref': row['ref'],
                        'alt': row['alt'],
                        'gene': row['gene'],
                        'frame': f"frame{frame}",
                        'strand': strand,
                        'wild_prot': wild_prot,
                        'mut_prot': mut_prot,
                    })
            
            if not mutation_found:
                print(f"No protein differences found for mutation at {row['chrom']}:{row['pos']}")
        
        print(f"Found {len(translation_results)} protein differences in {len(mut_df)} mutations")
        print(f"Saved {len(self.all_translations)} total protein translations for MS analysis")
        return pd.DataFrame(translation_results)

    def generate_all_peptides_for_ms(self):
        """从所有翻译序列中生成可能的肽段，用于质谱分析"""
        print("Generating all possible peptides for mass spectrometry analysis")
        all_ms_peptides = []
        
        # 遍历所有翻译
        for trans in self.all_translations:
            prot_seq = trans['prot_seq']
            
            # 只处理没有终止密码子的序列部分
            if '*' in prot_seq:
                # 分割序列，每个终止密码子前的部分独立处理
                segments = prot_seq.split('*')
                valid_segments = [seg for seg in segments if len(seg) >= self.long_range[0]]
            else:
                valid_segments = [prot_seq]
            
            # 为每个有效段生成所有可能的长肽段
            for segment in valid_segments:
                # 确保序列足够长
                if len(segment) < self.long_range[0]:
                    continue
                    
                # 滑动窗口生成各种长度的肽段
                for peptide_len in range(self.long_range[0], min(self.long_range[1] + 1, len(segment) + 1)):
                    for start in range(0, len(segment) - peptide_len + 1):
                        end = start + peptide_len
                        peptide = segment[start:end]
                        
                        # 确保肽段不含终止密码子
                        if '*' not in peptide:
                            all_ms_peptides.append({
                                'id': trans['id'],
                                'chrom': trans['chrom'], 
                                'pos': trans['pos'],
                                'ref': trans['ref'],
                                'alt': trans['alt'],
                                'gene': trans['gene'],
                                'frame': trans['frame'],
                                'strand': trans['strand'],
                                'peptide': peptide,
                                'length': peptide_len,
                                'start': start,
                                'end': end,
                                'is_mutated': trans['is_mutated']
                            })
        
        print(f"Generated {len(all_ms_peptides)} peptides for mass spectrometry analysis")
        return pd.DataFrame(all_ms_peptides)

    def find_mutation_position(self, trans_df):
        """找到突变位置"""
        print(f"Finding mutation positions in {len(trans_df)} translated sequences")
        mut_pos_results = []
        not_found = 0
        
        for idx, row in trans_df.iterrows():
            wild_prot = row['wild_prot']
            mut_prot = row['mut_prot']
            
            # 找出第一个不同的位置
            mut_pos = -1
            for i in range(min(len(wild_prot), len(mut_prot))):
                if wild_prot[i] != mut_prot[i]:
                    mut_pos = i
                    break
            
            # 如果长度不同且没找到不同，则突变位置在末尾
            if mut_pos == -1 and len(wild_prot) != len(mut_prot):
                mut_pos = min(len(wild_prot), len(mut_prot))
            
            if mut_pos >= 0:
                mut_pos_results.append({
                    'id': row['id'],
                    'chrom': row['chrom'],
                    'pos': row['pos'],
                    'ref': row['ref'],
                    'alt': row['alt'],
                    'gene': row['gene'],
                    'frame': row['frame'],
                    'strand': row['strand'],
                    'wild_prot': wild_prot,
                    'mut_prot': mut_prot,
                    'mut_pos': mut_pos
                })
            else:
                not_found += 1
                print(f"Warning: Could not identify mutation position for {row['chrom']}:{row['pos']}")
        
        print(f"Identified mutation positions for {len(mut_pos_results)} translations ({not_found} not found)")
        return pd.DataFrame(mut_pos_results)
    
    def generate_junction_peptides(self, mut_df):
        """生成跨越突变点的肽段"""
        print(f"Generating peptides spanning mutation sites for {len(mut_df)} mutations")
        all_peptides = []
        empty_mutations = 0
        
        for idx, row in mut_df.iterrows():
            mut_pos = row['mut_pos']
            mut_prot = row['mut_prot']
            
            # 生成短肽
            short_peptides = self._generate_peptides(
                mut_prot, mut_pos, self.short_range[0], self.short_range[1])
            
            # 生成长肽
            long_peptides = self._generate_peptides(
                mut_prot, mut_pos, self.long_range[0], self.long_range[1])
            
            total_peptides = len(short_peptides) + len(long_peptides)
            if total_peptides == 0:
                empty_mutations += 1
                print(f"Warning: No valid peptides found for mutation at {row['chrom']}:{row['pos']} in protein position {mut_pos}")
                
                # 检查是否有终止密码子阻止了肽段生成
                if '*' in mut_prot:
                    stop_pos = mut_prot.find('*')
                    if stop_pos <= mut_pos:
                        print(f"  Stop codon found at position {stop_pos}, before mutation position {mut_pos}")
                    else:
                        print(f"  Stop codon found at position {stop_pos}, after mutation position {mut_pos}")
            
            # 添加标识信息
            for peptide in short_peptides + long_peptides:
                peptide.update({
                    'id': row['id'],
                    'chrom': row['chrom'],
                    'pos': row['pos'],
                    'ref': row['ref'],
                    'alt': row['alt'],
                    'gene': row['gene'],
                    'frame': row['frame'],
                    'strand': row['strand'],
                    'mut_pos': mut_pos  # 确保把突变位置添加到每个肽段
                })
                all_peptides.append(peptide)
        
        short_count = len([p for p in all_peptides if p['peptide_type'] == 'short'])
        long_count = len([p for p in all_peptides if p['peptide_type'] == 'long'])
        print(f"Generated {len(all_peptides)} peptides spanning mutation sites:")
        print(f"  Short peptides (8-11 aa): {short_count}")
        print(f"  Long peptides (12-25 aa): {long_count}")
        print(f"  {empty_mutations} mutations did not produce any valid peptides")
        
        if long_count == 0:
            print("WARNING: No long peptides were generated. Possible reasons:")
            print("  1. All potential long peptides contain stop codons ('*')")
            print("  2. Mutations are close to protein termini or stop codons")
            print("  3. Sequence regions are too short for valid peptide generation")
        
        return pd.DataFrame(all_peptides)
    
    def _generate_peptides(self, sequence, mut_pos, min_len, max_len):
        """生成给定长度范围内跨越突变点的肽段"""
        peptides = []
        rejected = 0
        
        for peptide_len in range(min_len, max_len + 1):
            # 确保至少包含1个突变前和1个突变后的氨基酸
            for aa_before in range(1, peptide_len):
                aa_after = peptide_len - aa_before
                
                # 计算肽段的起始和结束位置
                start_pos = mut_pos - aa_before
                end_pos = mut_pos + aa_after
                
                # 检查位置是否有效
                if start_pos < 0 or end_pos > len(sequence):
                    continue
                
                # 提取肽段
                peptide = sequence[start_pos:end_pos]
                
                # 过滤含有终止密码子的肽段
                if '*' in peptide:
                    rejected += 1
                    continue
                
                # 检查是否包含突变点（确保peptide跨越突变位置）
                if not (start_pos <= mut_pos < end_pos):
                    continue
                
                # 添加元数据
                peptides.append({
                    'peptide': peptide,
                    'start': start_pos,
                    'end': end_pos,
                    'length': peptide_len,
                    'aa_before': aa_before,
                    'aa_after': aa_after,
                    'peptide_type': 'short' if peptide_len <= 11 else 'long'
                })
        
        return peptides
    
    def save_results(self, peptide_df):
        """保存结果到文件"""
        print("Saving results to files")
        
        # 保存完整结果到CSV
        csv_path = os.path.join(self.output_dir, "noncoding_mutation_peptides.csv")
        peptide_df.to_csv(csv_path, index=False)
        print(f"Saved all peptides to {csv_path}")
        
        # 按照长度类型分组
        for peptide_type in ['short', 'long']:
            type_df = peptide_df[peptide_df['peptide_type'] == peptide_type]
            
            if not type_df.empty:
                # 生成FASTA文件
                fasta_path = os.path.join(
                    self.output_dir, "fasta", f"noncoding_{peptide_type}_peptides.fasta")
                
                with open(fasta_path, 'w') as f:
                    for idx, row in type_df.iterrows():
                        gene_info = row['gene'] if row['gene'] != '.' else 'unknown'
                        header = f">{row['chrom']}|{row['pos']}|{row['ref']}>{row['alt']}|{gene_info}|{row['frame']}|{row['strand']}|pos={row['start']+1}-{row['end']}|len={row['length']}"
                        f.write(f"{header}\n{row['peptide']}\n")
                
                print(f"Saved {len(type_df)} {peptide_type} peptides to {fasta_path}")
        
        # 为质谱分析生成所有可能的长肽段
        ms_peptides_df = self.generate_all_peptides_for_ms()
        long_ms_peptides = ms_peptides_df[ms_peptides_df['length'] >= self.long_range[0]]
        
        if not long_ms_peptides.empty:
            # 生成用于质谱分析的FASTA文件
            ms_fasta_path = os.path.join(self.output_dir, "fasta", "noncoding_long_peptides_for_ms.fasta")
            
            with open(ms_fasta_path, 'w') as f:
                for idx, row in long_ms_peptides.iterrows():
                    gene_info = row['gene'] if row['gene'] != '.' else 'unknown'
                    mut_status = "mutated" if row['is_mutated'] else "unmutated"
                    header = f">{row['chrom']}|{row['pos']}|{row['ref']}>{row['alt']}|{gene_info}|{row['frame']}|{row['strand']}|{mut_status}|len={row['length']}"
                    f.write(f"{header}\n{row['peptide']}\n")
            
            print(f"Saved {len(long_ms_peptides)} long peptides for mass spectrometry analysis to {ms_fasta_path}")
            
            # 按基因分组统计突变肽段
            gene_counts = long_ms_peptides.groupby('gene').size()
            print("Long peptide counts by gene:")
            for gene, count in gene_counts.items():
                if gene != '.':
                    print(f"  {gene}: {count} peptides")
            
            # 保存统计信息到文本文件
            stats_path = os.path.join(self.output_dir, "long_peptide_stats.txt")
            with open(stats_path, 'w') as f:
                f.write(f"Total long peptides for MS analysis: {len(long_ms_peptides)}\n")
                f.write("Peptides by length:\n")
                length_counts = long_ms_peptides.groupby('length').size()
                for length, count in length_counts.items():
                    f.write(f"  Length {length}: {count} peptides\n")
                
                f.write("\nPeptides by gene:\n")
                for gene, count in gene_counts.items():
                    if gene != '.':
                        f.write(f"  {gene}: {count} peptides\n")
                
                # 添加突变/非突变统计
                mutated_count = len(long_ms_peptides[long_ms_peptides['is_mutated']])
                unmutated_count = len(long_ms_peptides[~long_ms_peptides['is_mutated']])
                f.write(f"\nMutated peptides: {mutated_count}\n")
                f.write(f"Unmutated peptides: {unmutated_count}\n")
            
            print(f"Saved long peptide statistics to {stats_path}")
    
    def run(self):
        """执行完整的处理流程"""
        # 1. 提取非编码区突变
        noncoding_df = self.extract_noncoding_mutations()
        
        # 2. 从参考基因组提取序列
        wild_sequences, fasta_path = self.extract_sequences(noncoding_df)
        
        # 3. 生成突变序列
        mut_df = self.generate_mutant_sequences(noncoding_df, wild_sequences)
        
        # 4. 翻译并比较序列
        trans_df = self.translate_and_compare(mut_df)
        
        # 5. 找到突变位置
        mut_pos_df = self.find_mutation_position(trans_df)
        
        # 6. 生成肽段
        peptide_df = self.generate_junction_peptides(mut_pos_df)
        
        # 输出统计信息
        total_short_peptides = len(peptide_df[peptide_df['peptide_type'] == 'short'])
        total_long_peptides = len(peptide_df[peptide_df['peptide_type'] == 'long'])
        print(f"\nGenerated peptide statistics:")
        print(f"  Short peptides (8-11 aa): {total_short_peptides}")
        print(f"  Long peptides (12-25 aa): {total_long_peptides}")
        print(f"  Total peptides: {len(peptide_df)}")
        
        # 7. 保存结果
        self.save_results(peptide_df)
        
        print("Noncoding mutation peptide generation completed successfully")

def main():
    parser = argparse.ArgumentParser(description='Generate peptides spanning noncoding mutation sites')
    parser.add_argument('-i', '--input', required=True, help='ANNOVAR annotation file')
    parser.add_argument('-r', '--reference', required=True, help='Reference genome FASTA file')
    parser.add_argument('-o', '--output', required=True, help='Output directory')
    parser.add_argument('--short_min', type=int, default=8, help='Minimum length for short peptides')
    parser.add_argument('--short_max', type=int, default=11, help='Maximum length for short peptides')
    parser.add_argument('--long_min', type=int, default=12, help='Minimum length for long peptides')
    parser.add_argument('--long_max', type=int, default=25, help='Maximum length for long peptides')
    parser.add_argument('--no_stop_codons', action='store_true', help='Filter out peptides containing stop codons')
    parser.add_argument('--try_fix_mismatches', action='store_true', 
                      help='Try to fix reference mismatches by checking nearby positions')
    parser.add_argument('--debug', action='store_true', help='Enable detailed debug logging')
    
    args = parser.parse_args()
    
    generator = NoncodingPeptideGenerator(
        args.input,
        args.reference,
        args.output,
        short_peptide_len=(args.short_min, args.short_max),
        long_peptide_len=(args.long_min, args.long_max)
    )
    
    generator.run()

if __name__ == "__main__":
    main()