#!/usr/bin/env nextflow

// RNA-Seq分析工作流

// 导入模块
// MERGE_READS
include { MERGE_READS } from '../modules/process_merge'
// FASTP
include { FASTP } from '../modules/fastp'
// bwa
include { BWA_INDEX } from '../modules/bwa'
include { BWA } from '../modules/bwa'
//samtools
include { SAMTOOLS_FIXMATE } from '../modules/samtools'
include { SAMTOOLS_SORT } from '../modules/samtools'
include { SAMTOOLS_INDEX } from '../modules/samtools'
// GATK
include { GATK_MARKDUPLICATES } from '../modules/gatk'
include { GATK_ADD_OR_REPLACE_READ_GROUPS } from '../modules/gatk'
include { GATK_BASE_RECALIBRATOR } from '../modules/gatk'
include { GATK_APPLY_BQSR } from '../modules/gatk'
include { GATK_MUTECT2_NORMAL } from '../modules/gatk'
include { GATK_MUTECT2_PAIRED } from '../modules/gatk'
include { GATK_FILTER_MUTECT_CALLS } from '../modules/gatk'
include { BCFTOOLS_FILTER_PASS } from '../modules/gatk'
// ANNOVAR
include { ANNOVAR } from '../modules/annovar'
// Picard
include { PICARD_ADD_OR_REPLACE_READ_GROUPS } from '../modules/picard'
// STAR-Fusion
include { STAR_FUSION } from '../modules/star_fusion'
// Kallisto
include { KALLISTO_INDEX } from '../modules/kallisto'
include { KALLISTO_QUANT } from '../modules/kallisto'
// HLA-I Typing
include { HLA_I_WORKFLOW } from '../modules/HLA1OptiType'
// HLA-II Typing
include { HLAMINER_WORKFLOW } from '../modules/HLAminer'

workflow RNA_WORKFLOW {
    take:
    raw_samples_ch  // 原始样本通道 [sample_id, dna_normal_r1, dna_normal_r2, dna_tumor_r1, dna_tumor_r2, rna_normal_r1, rna_normal_r2, rna_tumor_r1, rna_tumor_r2]
    reference_file  // 参考基因组文件

    main:
    // 准备仅包含RNA数据的输入通道
    rna_samples_ch = raw_samples_ch.map { id, dn_r1, dn_r2, dt_r1, dt_r2, rn_r1, rn_r2, rt_r1, rt_r2 ->
        // 将空字符串替换为null，以便process_merge正确处理
        def nullify = { str -> str == '' ? null : str }
        tuple(id, 
              null, null, null, null,  // DNA数据设置为null
              nullify(rn_r1), nullify(rn_r2), nullify(rt_r1), nullify(rt_r2)) // RNA数据
    }
    
    // 首先使用 MERGE_READS 模块合并所有分割文件，但仅处理RNA数据
    merged_reads = MERGE_READS(rna_samples_ch)
    
    // 从合并后的读数中提取 RNA 数据
    rna_normal_reads = merged_reads.rna_normal_paired
    rna_tumor_reads = merged_reads.rna_tumor_paired
    
    // 准备用于 FASTP 的拆分样本通道
    rna_split_samples = rna_normal_reads
        .map { id, r1, r2 -> tuple("${id}_rna_normal", r1, r2) }
        .mix(
            rna_tumor_reads.map { id, r1, r2 -> tuple("${id}_rna_tumor", r1, r2) }
        )
    
    // 准备用于配对分析的样本通道
    rna_paired_samples = rna_normal_reads
        .map { id, r1, r2 -> tuple(id, ['normal': [r1, r2]]) }
        .join(
            rna_tumor_reads.map { id, r1, r2 -> tuple(id, ['tumor': [r1, r2]]) }
        )
        .map { id, normal, tumor -> tuple(id, normal + tumor) }
    
    // 调用 FASTP 模块处理原始读数
    fastp_result = FASTP(rna_split_samples)
    
    // 提取 FASTP 模块中经过质控后的配对读数通道
    cleaned_reads_ch = fastp_result.cleaned_reads
    
    // 为RNA肿瘤样本运行HLA-I分型 - 仅处理肿瘤样本
    rna_tumor_for_hla = cleaned_reads_ch
        .filter { sample_id, r1, r2 -> sample_id.toString().contains('_rna_tumor') }
    
    // 运行HLA-I分型工作流
    hla_typing_results = HLA_I_WORKFLOW(rna_tumor_for_hla)
    
    // 运行HLA-II分型工作流 (HLAminer) - 同样只处理肿瘤样本
    hlaminer_typing_results = HLAMINER_WORKFLOW(rna_tumor_for_hla)
    
    // 为RNA肿瘤样本运行STAR-Fusion
    // 仅选择肿瘤样本进行融合检测
    rna_tumor_cleaned_reads = cleaned_reads_ch
        .filter { sample_id, r1, r2 -> sample_id.toString().contains('_rna_tumor') }
    
    // 运行STAR-Fusion
    star_fusion_results = STAR_FUSION(rna_tumor_cleaned_reads)
    
    // 运行Kallisto表达量定量
    // 1. 构建Kallisto索引
    transcriptome_fasta = file(params.kallisto.transcriptome_fasta)
    kallisto_index_result = KALLISTO_INDEX(transcriptome_fasta)
    
    // 2. 对RNA样本进行定量
    kallisto_quant_results = KALLISTO_QUANT(
        cleaned_reads_ch,
        kallisto_index_result.kallisto_index
    )
    
    // 创建 BWA 索引
    bwa_index_result = BWA_INDEX(reference_file)

    // 将 FASTP 的输出和参考基因组及索引文件传递给 BWA 模块做比对
    bwa_result = BWA(
        cleaned_reads_ch,        // 清洗后的reads
        reference_file,          // 参考基因组文件
        bwa_index_result.index   // 索引文件
    )

    // 将 BWA 输出的对齐结果传递给 SAMTOOLS_FIXMATE 进程
    fixmate_result = SAMTOOLS_FIXMATE(bwa_result)

    // 将 SAMTOOLS_FIXMATE 输出的结果传递给 SAMTOOLS_SORT 进程
    sort_result = SAMTOOLS_SORT(fixmate_result)
    
    // 将排序后的 BAM 文件传递给 GATK MarkDuplicates 进程
    markdup_result = GATK_MARKDUPLICATES(sort_result)
    
    // 将标记重复后的 BAM 文件传递给 GATK AddOrReplaceReadGroups 进程
    readgroups_result = GATK_ADD_OR_REPLACE_READ_GROUPS(markdup_result.marked_bam)
    
    // 将添加读组后的 BAM 文件传递给 GATK BaseRecalibrator 进程
    recal_result = GATK_BASE_RECALIBRATOR(
        readgroups_result.bam_with_read_groups,
        reference_file
    )
    
    // 将重校准数据和BAM文件传递给 GATK ApplyBQSR 进程
    bqsr_result = GATK_APPLY_BQSR(
        readgroups_result.bam_with_read_groups.join(recal_result.recal_table),
        reference_file
    )
    
    // 将重校准后的BAM文件传递给 Picard AddOrReplaceReadGroups 进程
    picard_result = PICARD_ADD_OR_REPLACE_READ_GROUPS(bqsr_result.recalibrated_bam)
    
    // 为最终的BAM文件创建索引
    index_result = SAMTOOLS_INDEX(picard_result.final_bam)
    
    // 将最终的BAM文件和索引文件组合在一起
    final_bams_with_index = picard_result.final_bam
        .map { sample_id, bam -> tuple(sample_id, bam) }
        .join(index_result.indexed_bam
            .map { sample_id, bam, bai -> tuple(sample_id, bai) }
        )
        .map { sample_id, bam, bai -> tuple(sample_id, bam, bai) }
    
    // 将最终的BAM文件分为normal和tumor通道
    final_bams_normal = final_bams_with_index
        .filter { it[0].toString().endsWith('_rna_normal') }
    
    final_bams_tumor = final_bams_with_index
        .filter { it[0].toString().endsWith('_rna_tumor') }
    
    // 对normal样本运行Mutect2
    normal_vcf = GATK_MUTECT2_NORMAL(
        final_bams_normal,
        reference_file
    )
    
    // 准备配对分析的输入数据
    paired_input = final_bams_tumor
        .map { tumor_id, tumor_bam, tumor_bai -> 
            def base_id = tumor_id.replaceAll('_rna_tumor$', '')
            tuple(base_id, tumor_id, tumor_bam, tumor_bai)
        }
        .join(
            final_bams_normal.map { normal_id, normal_bam, normal_bai ->
                def base_id = normal_id.replaceAll('_rna_normal$', '')
                tuple(base_id, normal_id, normal_bam, normal_bai)
            }
        )
        .join(
            normal_vcf.normal_vcf.map { normal_id, vcf, idx ->
                def base_id = normal_id.replaceAll('_rna_normal$', '')
                tuple(base_id, vcf, idx)
            }
        )
        .map { base_id, 
               tumor_id, tumor_bam, tumor_bai,
               normal_id, normal_bam, normal_bai,
               vcf, idx ->
            tuple(tumor_id, tumor_bam, tumor_bai,
                  normal_id, normal_bam, normal_bai,
                  vcf, idx)
        }
    
    // 运行配对样本分析
    paired_vcf = GATK_MUTECT2_PAIRED(
        paired_input,
        reference_file
    )
    
    // 过滤突变结果
    filtered_vcf = GATK_FILTER_MUTECT_CALLS(
        paired_vcf.paired_vcf,
        reference_file
    )
    
    // 使用bcftools进一步过滤PASS的变异
    pass_filtered_vcf = BCFTOOLS_FILTER_PASS(
        filtered_vcf.filtered_vcf
    )
    
    // 使用 ANOVAR 注释 PASS 过滤后的变异
    anovar_results = ANNOVAR(
        pass_filtered_vcf
    )
    
 
    // 输出结果
    emit:
    fusion_results = star_fusion_results.fusion_results  // 确保通道名称正确
    kallisto_abundance = kallisto_quant_results.quant_results
    hla_results = hla_typing_results.hla_results
    hla2_results = hlaminer_typing_results.hla_results  // HLA-II分型结果
    annotation_results = anovar_results.anovar_txt
    annotation_vcf = anovar_results.anovar_vcf
}
