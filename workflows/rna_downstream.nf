#!/usr/bin/env nextflow

// RNA分析的下游分析工作流，处理注释结果和kallisto定量结果，进行后续表达过滤和肽段生成

// 导入模块
include { RNA_TPM_FILTER } from '../modules/downstream/rnatpmfilter'
include { RNA_NONCODING_ANALYSIS } from '../modules/downstream/rnanoncoding'
include { RNA_STAR_FUSION_PEPTIDES } from '../modules/downstream/rnastarfusion'

workflow RNA_DOWNSTREAM {
    take:
    annovar_results         // ANNOVAR的RNA突变注释结果 [sample_id, annovar_txt]
    kallisto_abundance      // Kallisto定量结果 [sample_id, abundance_tsv]
    star_fusion_results     // STAR-Fusion结果 [sample_id, fusion_results]
    
    main:
    // 步骤1: 基于基因表达值进行RNA变异过滤
    // 需要转录本-基因映射文件 - 使用配置文件中已定义的路径
    gene_ref_file = file(params.kallisto.gene_transcript_map)
    
    // 应用TPM过滤
    tpm_filtered_results = RNA_TPM_FILTER(
        kallisto_abundance,
        gene_ref_file,
        annovar_results
    )
    
    // 步骤2: 对过滤后的变异进行非编码区突变分析
    noncoding_results = RNA_NONCODING_ANALYSIS(
        tpm_filtered_results.filtered_variants
    )
    
    // 步骤3: 生成融合肽段
    // 检查star_fusion_results是否为空
    star_fusion_results_exists = star_fusion_results
        .map { it -> return it }
        .ifEmpty { log.warn("No STAR-Fusion results provided, skipping fusion peptide generation"); return Channel.empty() }
    
    // 仅当存在融合结果时运行肽段生成过程
    fusion_results = RNA_STAR_FUSION_PEPTIDES(
        star_fusion_results_exists
    )
    
    // 准备空通道，以确保即使没有融合结果，输出通道也能正常工作
    empty_fusion_short = Channel.empty()
    empty_fusion_long = Channel.empty()
    empty_fusion_short_csv = Channel.empty()
    empty_fusion_long_csv = Channel.empty()
    
    // 合并实际结果与空通道
    fusion_short_peptides = fusion_results.fusion_short_peptides.ifEmpty { empty_fusion_short }
    fusion_long_peptides = fusion_results.fusion_long_peptides.ifEmpty { empty_fusion_long }
    fusion_short_csv = fusion_results.fusion_short_csv.ifEmpty { empty_fusion_short_csv }
    fusion_long_csv = fusion_results.fusion_long_csv.ifEmpty { empty_fusion_long_csv }
    
    // 整合所有结果
    
    emit:
    // TPM过滤结果
    filtered_variants = tpm_filtered_results.filtered_variants
    gene_expression = tpm_filtered_results.gene_expression
    
    // 非编码区肽段
    noncoding_peptides_csv = noncoding_results.noncoding_peptides_csv
    noncoding_short_peptides = noncoding_results.noncoding_short_peptides
    noncoding_long_peptides = noncoding_results.noncoding_long_peptides
    noncoding_ms_peptides = noncoding_results.noncoding_ms_peptides
    noncoding_stats = noncoding_results.noncoding_stats
    
    // 融合肽段
    fusion_short_peptides = fusion_short_peptides
    fusion_long_peptides = fusion_long_peptides
    fusion_short_csv = fusion_short_csv
    fusion_long_csv = fusion_long_csv
    
    // MS分析需要的肽段通道
    noncoding_peptides = noncoding_results.noncoding_ms_peptides  // 非编码区肽段用于质谱分析
    fusion_peptides = fusion_long_peptides                        // 融合长肽段用于质谱分析
    
    // 所有肽段的集合（可用于后续整合分析）
    all_short_peptides = noncoding_results.noncoding_short_peptides
        .mix(fusion_short_peptides)
    all_long_peptides = noncoding_results.noncoding_long_peptides
        .mix(fusion_long_peptides)
}
