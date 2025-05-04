#!/usr/bin/env nextflow

// RNA分析的下游分析工作流，处理注释结果和kallisto定量结果，进行后续表达过滤和肽段生成

// 导入模块
include { RNA_TPM_FILTER } from '../modules/downstream/rnatpmfilter'
include { RNA_NONCODING_ANALYSIS } from '../modules/downstream/rnanoncoding'
include { RNA_STAR_FUSION_PEPTIDES } from '../modules/downstream/rnastarfusion'

// 处理缺失文件的辅助进程
process CREATE_EMPTY_FILES {
    tag "$sample_id"
    
    input:
    val(sample_id)
    
    output:
    tuple val(sample_id), path("${sample_id}_empty_noncoding_ms_peptides.fasta"), emit: empty_noncoding_ms
    tuple val(sample_id), path("${sample_id}_empty_fusion_long_peptides.fasta"), emit: empty_fusion_long
    tuple val(sample_id), path("${sample_id}_empty_noncoding_short.fasta"), emit: empty_noncoding_short
    tuple val(sample_id), path("${sample_id}_empty_noncoding_long.fasta"), emit: empty_noncoding_long
    tuple val(sample_id), path("${sample_id}_empty_noncoding_csv.csv"), emit: empty_noncoding_csv
    tuple val(sample_id), path("${sample_id}_empty_noncoding_stats.txt"), emit: empty_noncoding_stats
    tuple val(sample_id), path("${sample_id}_empty_fusion_short.fasta"), emit: empty_fusion_short
    tuple val(sample_id), path("${sample_id}_empty_fusion_short.csv"), emit: empty_fusion_short_csv
    tuple val(sample_id), path("${sample_id}_empty_fusion_long.csv"), emit: empty_fusion_long_csv
    
    script:
    """
    # 创建空的非编码区肽段文件
    echo ">Empty_noncoding_${sample_id}" > ${sample_id}_empty_noncoding_ms_peptides.fasta
    echo "EMPTYPEPTIDE" >> ${sample_id}_empty_noncoding_ms_peptides.fasta
    
    # 创建空的融合肽段文件
    echo ">Empty_fusion_${sample_id}" > ${sample_id}_empty_fusion_long_peptides.fasta
    echo "EMPTYPEPTIDE" >> ${sample_id}_empty_fusion_long_peptides.fasta
    
    # 创建其他空文件
    echo ">Empty_noncoding_short_${sample_id}" > ${sample_id}_empty_noncoding_short.fasta
    echo "EMPTY" >> ${sample_id}_empty_noncoding_short.fasta
    
    echo ">Empty_noncoding_long_${sample_id}" > ${sample_id}_empty_noncoding_long.fasta
    echo "EMPTYPEPTIDE" >> ${sample_id}_empty_noncoding_long.fasta
    
    echo "SampleID,GeneID,GeneName,MutationType,Count" > ${sample_id}_empty_noncoding_csv.csv
    echo "${sample_id},NA,NA,NA,0" >> ${sample_id}_empty_noncoding_csv.csv
    
    echo "No noncoding mutations found for ${sample_id}" > ${sample_id}_empty_noncoding_stats.txt
    
    echo ">Empty_fusion_short_${sample_id}" > ${sample_id}_empty_fusion_short.fasta
    echo "EMPTY" >> ${sample_id}_empty_fusion_short.fasta
    
    echo "SampleID,Gene1,Gene2,Count,Peptides" > ${sample_id}_empty_fusion_short.csv
    echo "${sample_id},NA,NA,0,0" >> ${sample_id}_empty_fusion_short.csv
    
    echo "SampleID,Gene1,Gene2,Count,Peptides" > ${sample_id}_empty_fusion_long.csv
    echo "${sample_id},NA,NA,0,0" >> ${sample_id}_empty_fusion_long.csv
    """
}

// 创建一个合并肽段的进程
process MERGE_RNA_PEPTIDES {
    tag "${sample_id}"
    
    publishDir "${params.outdir}/${sample_id}/rna/12_integrated", mode: 'copy'
    
    input:
    tuple val(sample_id), path(noncoding_short), path(fusion_short)
    tuple val(sample_id), path(noncoding_long), path(fusion_long)
    
    output:
    tuple val(sample_id), path("${sample_id}_all_short_peptides.fasta"), emit: all_short_peptides
    tuple val(sample_id), path("${sample_id}_all_long_peptides.fasta"), emit: all_long_peptides
    
    script:
    """
    # 合并短肽段文件
    if [ -s "${noncoding_short}" ] || [ -s "${fusion_short}" ]; then
        cat "${noncoding_short}" "${fusion_short}" > "${sample_id}_all_short_peptides.fasta"
    else
        echo ">Empty_all_short_peptides_${sample_id}" > "${sample_id}_all_short_peptides.fasta"
        echo "EMPTY" >> "${sample_id}_all_short_peptides.fasta"
    fi
    
    # 合并长肽段文件
    if [ -s "${noncoding_long}" ] || [ -s "${fusion_long}" ]; then
        cat "${noncoding_long}" "${fusion_long}" > "${sample_id}_all_long_peptides.fasta"
    else
        echo ">Empty_all_long_peptides_${sample_id}" > "${sample_id}_all_long_peptides.fasta"
        echo "EMPTYPEPTIDE" >> "${sample_id}_all_long_peptides.fasta"
    fi
    
    # 输出统计信息
    echo "Merged RNA peptides for sample ${sample_id}:"
    echo "Total short peptides: \$(grep -c '>' ${sample_id}_all_short_peptides.fasta)"
    echo "Total long peptides: \$(grep -c '>' ${sample_id}_all_long_peptides.fasta)"
    """
}

workflow RNA_DOWNSTREAM {
    take:
    annovar_results         // ANNOVAR的RNA突变注释结果 [sample_id, annovar_txt]
    kallisto_abundance      // Kallisto定量结果 [sample_id, abundance_tsv]
    star_fusion_results     // STAR-Fusion结果 [sample_id, fusion_results]
    
    main:
    // 提取样本ID列表，以便处理缺失文件
    sample_ids = annovar_results.map { it[0] }
    
    // 创建空文件通道
    empty_files = CREATE_EMPTY_FILES(sample_ids)
    
    // 步骤1: 基于基因表达值进行RNA变异过滤
    // 需要转录本-基因映射文件
    gene_ref_file = file(params.kallisto.gene_transcript_map)
    
    // 应用TPM过滤，确保所有输入都存在
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
    fusion_results = RNA_STAR_FUSION_PEPTIDES(
        star_fusion_results
    )
    
    // 处理可能缺失的结果
    // 非编码区肽段文件
    noncoding_ms_peptides_ch = noncoding_results.noncoding_ms_peptides
        .mix(empty_files.empty_noncoding_ms)
        .groupTuple(size: 1, remainder: true)
        .map { sample_id, files -> 
            if (files.size() > 1) {
                // 如果有多个文件，保留非空文件
                def nonEmpty = files.find { it.size() > 0 }
                return tuple(sample_id, nonEmpty ?: files[0])
            } else {
                return tuple(sample_id, files[0])
            }
        }
    
    // 融合肽段文件
    fusion_long_peptides_ch = fusion_results.fusion_long_peptides
        .mix(empty_files.empty_fusion_long)
        .groupTuple(size: 1, remainder: true)
        .map { sample_id, files -> 
            if (files.size() > 1) {
                // 如果有多个文件，保留非空文件
                def nonEmpty = files.find { it.size() > 0 }
                return tuple(sample_id, nonEmpty ?: files[0])
            } else {
                return tuple(sample_id, files[0])
            }
        }
    
    // 同样处理其他可能缺失的结果
    noncoding_short_ch = noncoding_results.noncoding_short_peptides
        .mix(empty_files.empty_noncoding_short)
        .groupTuple(size: 1, remainder: true)
        .map { sample_id, files -> 
            if (files.size() > 1) {
                def nonEmpty = files.find { it.size() > 0 }
                return tuple(sample_id, nonEmpty ?: files[0])
            } else {
                return tuple(sample_id, files[0])
            }
        }
    
    noncoding_long_ch = noncoding_results.noncoding_long_peptides
        .mix(empty_files.empty_noncoding_long)
        .groupTuple(size: 1, remainder: true)
        .map { sample_id, files -> 
            if (files.size() > 1) {
                def nonEmpty = files.find { it.size() > 0 }
                return tuple(sample_id, nonEmpty ?: files[0])
            } else {
                return tuple(sample_id, files[0])
            }
        }
    
    fusion_short_ch = fusion_results.fusion_short_peptides
        .mix(empty_files.empty_fusion_short)
        .groupTuple(size: 1, remainder: true)
        .map { sample_id, files -> 
            if (files.size() > 1) {
                def nonEmpty = files.find { it.size() > 0 }
                return tuple(sample_id, nonEmpty ?: files[0])
            } else {
                return tuple(sample_id, files[0])
            }
        }
    
    // 整合所有结果: 合并RNA肽段文件
    merged_peptides = MERGE_RNA_PEPTIDES(
        noncoding_short_ch.join(fusion_short_ch),
        noncoding_long_ch.join(fusion_long_peptides_ch)
    )
    
    emit:
    // TPM过滤结果
    filtered_variants = tpm_filtered_results.filtered_variants
    gene_expression = tpm_filtered_results.gene_expression
    
    // 非编码区肽段
    noncoding_peptides_csv = noncoding_results.noncoding_peptides_csv
    noncoding_short_peptides = noncoding_short_ch
    noncoding_long_peptides = noncoding_long_ch
    noncoding_ms_peptides = noncoding_results.noncoding_ms_peptides
    noncoding_stats = noncoding_results.noncoding_stats
    
    // 融合肽段
    fusion_short_peptides = fusion_short_ch
    fusion_long_peptides = fusion_results.fusion_long_peptides
    fusion_short_csv = fusion_results.fusion_short_csv
    fusion_long_csv = fusion_results.fusion_long_csv
    
    // MS分析需要的肽段通道
    noncoding_peptides = noncoding_ms_peptides_ch  // 非编码区肽段用于质谱分析
    fusion_peptides = fusion_long_peptides_ch      // 融合长肽段用于质谱分析
    
    // 所有肽段的集合（可用于后续整合分析）
    all_short_peptides = merged_peptides.all_short_peptides
    all_long_peptides = merged_peptides.all_long_peptides
}
