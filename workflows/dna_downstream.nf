#!/usr/bin/env nextflow

// DNA分析的下游分析工作流，处理注释结果并进行后续分析

// 导入模块
include { SNV_ANALYSIS } from '../modules/downstream/SNV'
include { FRAMESHIFT_DELETION_ANALYSIS } from '../modules/downstream/frameshift_deletion'
include { NONFRAMESHIFT_DELETION_ANALYSIS } from '../modules/downstream/nonframeshift_deletion'
include { NONFRAMESHIFT_INSERTION_ANALYSIS } from '../modules/downstream/nonframeshift_insertion'
include { FRAMESHIFT_INSERTION_ANALYSIS as FRAMESHIFT_INSERTION } from '../modules/downstream/frameshift_insertion'
include { NONFRAMESHIFT_SUBSTITUTION_ANALYSIS as NONFRAMESHIFT_SUBSTITUTION } from '../modules/downstream/nonframeshift_substitution'
include { STOPLOSS_MUTATION_ANALYSIS as STOPLOSS_MUTATION }
// 未来可能添加的其他下游分析模块:
// include { INDEL_ANALYSIS } from '../modules/downstream/INDEL'
// include { CNV_ANALYSIS } from '../modules/downstream/CNV'
// include { PATHWAY_ANALYSIS } from '../modules/downstream/PATHWAY'

// 创建用于整合结果的进程
process INTEGRATE_RESULTS {
    tag "$sample_id"
    label 'process_low'
    
    publishDir "${params.outdir}/${sample_id}/dna/11_integrated", mode: 'copy'
    
    input:
    tuple val(sample_id), path(gene_summaries)
    tuple val(sample_id), path(short_peptides)
    tuple val(sample_id), path(long_peptides)
    
    output:
    tuple val(sample_id), path("${sample_id}_integrated_gene_summary.csv"), emit: integrated_summary
    tuple val(sample_id), path("${sample_id}_integrated_short_peptides.fasta"), emit: integrated_short_peptides
    tuple val(sample_id), path("${sample_id}_integrated_long_peptides.fasta"), emit: integrated_long_peptides
    
    script:
    """
    # 检查是否有任何基因摘要文件
    if [ \$(find . -name "*.csv" -type f | wc -l) -eq 0 ]; then
        # 创建空的基因摘要，如果没有有效的输入文件
        echo "SampleID,GeneID,GeneName,MutationType,VariantCount,AnchorCount,PeptideCount" > ${sample_id}_integrated_gene_summary.csv
        echo "${sample_id},NA,NA,NA,0,0,0" >> ${sample_id}_integrated_gene_summary.csv
    else
        # 合并基因摘要CSV文件
        echo "SampleID,GeneID,GeneName,MutationType,VariantCount,AnchorCount,PeptideCount" > ${sample_id}_integrated_gene_summary.csv
        for file in ${gene_summaries}; do
            if [ -s "\$file" ]; then
                tail -n +2 \$file >> ${sample_id}_integrated_gene_summary.csv
            fi
        done
    fi
    
    # 检查是否有短肽段文件
    if [ \$(find . -name "*_short_peptides.fasta" -type f | wc -l) -eq 0 ]; then
        # 创建空的短肽段文件
        echo ">Empty_short_peptide_${sample_id}" > ${sample_id}_integrated_short_peptides.fasta
        echo "EMPTY" >> ${sample_id}_integrated_short_peptides.fasta
    else
        # 合并有效的短肽段FASTA文件
        touch ${sample_id}_integrated_short_peptides.fasta
        for file in ${short_peptides}; do
            if [ -s "\$file" ]; then
                cat \$file >> ${sample_id}_integrated_short_peptides.fasta
            fi
        done
        # 如果合并后依然为空，添加一个空记录
        if [ ! -s "${sample_id}_integrated_short_peptides.fasta" ]; then
            echo ">Empty_short_peptide_${sample_id}" > ${sample_id}_integrated_short_peptides.fasta
            echo "EMPTY" >> ${sample_id}_integrated_short_peptides.fasta
        fi
    fi
    
    # 检查是否有长肽段文件
    if [ \$(find . -name "*_long_peptides.fasta" -type f | wc -l) -eq 0 ]; then
        # 创建空的长肽段文件
        echo ">Empty_long_peptide_${sample_id}" > ${sample_id}_integrated_long_peptides.fasta
        echo "EMPTYPEPTIDE" >> ${sample_id}_integrated_long_peptides.fasta
    else
        # 合并有效的长肽段FASTA文件
        touch ${sample_id}_integrated_long_peptides.fasta
        for file in ${long_peptides}; do
            if [ -s "\$file" ]; then
                cat \$file >> ${sample_id}_integrated_long_peptides.fasta
            fi
        done
        # 如果合并后依然为空，添加一个空记录
        if [ ! -s "${sample_id}_integrated_long_peptides.fasta" ]; then
            echo ">Empty_long_peptide_${sample_id}" > ${sample_id}_integrated_long_peptides.fasta
            echo "EMPTYPEPTIDE" >> ${sample_id}_integrated_long_peptides.fasta
        fi
    fi
    
    # 输出统计信息
    echo "Integrated results for sample ${sample_id}:"
    echo "Total gene entries: \$(tail -n +2 ${sample_id}_integrated_gene_summary.csv | wc -l)"
    echo "Total short peptides: \$(grep -c '>' ${sample_id}_integrated_short_peptides.fasta)"
    echo "Total long peptides: \$(grep -c '>' ${sample_id}_integrated_long_peptides.fasta)"
    """
}

// 处理frameshift deletion肽段文件以提取junction肽段
process EXTRACT_FSDEL_JUNCTION_PEPTIDES {
    tag "$sample_id"
    label 'process_low'
    
    publishDir "${params.outdir}/${sample_id}/dna/11_integrated/ms_peptides", mode: 'copy'
    
    input:
    tuple val(sample_id), path(fsdel_files)
    
    output:
    tuple val(sample_id), 
          path("${sample_id}_fs_del_junction_peptides_1.fasta"),
          path("${sample_id}_fs_del_junction_peptides_2.fasta"),
          emit: fsdel_junction_peptides
    
    script:
    """
    # 查找junction肽段文件
    fs_del_junction_1=\$(find . -name "*fs_del_junction_peptides_1.fasta" || echo "")
    fs_del_junction_2=\$(find . -name "*fs_del_junction_peptides_2.fasta" || echo "")
    
    # 如果找到junction肽段文件1，则复制它
    if [ -n "\$fs_del_junction_1" ] && [ -s "\$fs_del_junction_1" ]; then
        cp \$fs_del_junction_1 ${sample_id}_fs_del_junction_peptides_1.fasta
    else
        # 否则创建空文件
        echo ">Empty_fs_del_junction_1_${sample_id}" > ${sample_id}_fs_del_junction_peptides_1.fasta
        echo "EMPTYPEPTIDE" >> ${sample_id}_fs_del_junction_peptides_1.fasta
    fi
    
    # 如果找到junction肽段文件2，则复制它
    if [ -n "\$fs_del_junction_2" ] && [ -s "\$fs_del_junction_2" ]; then
        cp \$fs_del_junction_2 ${sample_id}_fs_del_junction_peptides_2.fasta
    else
        # 否则创建空文件
        echo ">Empty_fs_del_junction_2_${sample_id}" > ${sample_id}_fs_del_junction_peptides_2.fasta
        echo "EMPTYPEPTIDE" >> ${sample_id}_fs_del_junction_peptides_2.fasta
    fi
    """
}

// 处理其他类型的肽段文件，提取MS肽段
process EXTRACT_MS_PEPTIDES {
    tag "$sample_id"
    label 'process_low'
    
    publishDir "${params.outdir}/${sample_id}/dna/11_integrated/ms_peptides", mode: 'copy'
    
    input:
    tuple val(sample_id), path(var_sequence_files), val(mutation_type)
    
    output:
    tuple val(sample_id), path("${sample_id}_${mutation_type}_ms_peptides.fasta"), emit: ms_peptides
    
    script:
    """
    # 查找变异序列文件
    var_seq_file=\$(find . -name "*${mutation_type}*Varsequence.fasta" || echo "")
    
    # 如果找到变异序列文件，则复制它
    if [ -n "\$var_seq_file" ] && [ -s "\$var_seq_file" ]; then
        cp \$var_seq_file ${sample_id}_${mutation_type}_ms_peptides.fasta
    else
        # 否则创建空文件
        echo ">Empty_${mutation_type}_${sample_id}" > ${sample_id}_${mutation_type}_ms_peptides.fasta
        echo "EMPTYPEPTIDE" >> ${sample_id}_${mutation_type}_ms_peptides.fasta
    fi
    """
}

workflow DNA_DOWNSTREAM {
    take:
    annotation_results_ch
    annotation_vcf_ch
    
    main:
    // 对注释结果进行SNV分析
    snv_results = SNV_ANALYSIS(
        annotation_results_ch
    )
    
    // 对注释结果进行Frameshift Deletion分析
    frameshift_del_results = FRAMESHIFT_DELETION_ANALYSIS(
        annotation_results_ch
    )
    
    // 对注释结果进行Nonframeshift Deletion分析
    nonframeshift_del_results = NONFRAMESHIFT_DELETION_ANALYSIS(
        annotation_results_ch
    )
    
    // 对注释结果进行Nonframeshift Insertion分析
    nonframeshift_ins_results = NONFRAMESHIFT_INSERTION_ANALYSIS(
        annotation_results_ch
    )
    
    // 对注释结果进行Frameshift Insertion分析
    frameshift_insertion_results = FRAMESHIFT_INSERTION(
        annotation_results_ch
    )
    
    // 对注释结果进行Nonframeshift Substitution分析
    nonframeshift_substitution_results = NONFRAMESHIFT_SUBSTITUTION(
        annotation_results_ch
    )
    
    // 对注释结果进行Stoploss Mutation分析
    stoploss_mutation_results = STOPLOSS_MUTATION(
        annotation_results_ch
    )
    
    // 收集所有基因摘要文件
    all_gene_summaries = frameshift_insertion_results.fsins_gene_summary
        .mix(nonframeshift_substitution_results.nfssub_gene_summary)
        .mix(stoploss_mutation_results.stoploss_gene_summary)
        .mix(snv_results.snv_gene_summary)
        .mix(frameshift_del_results.fsdel_gene_summary)
        .mix(nonframeshift_del_results.nfsdel_gene_summary)
        .mix(nonframeshift_ins_results.nfsins_gene_summary)
        .groupTuple()
    
    // 收集所有短肽段文件
    all_short_peptides = frameshift_insertion_results.fsins_fasta_files
        .map { sample_id, files -> 
            def short_peptides = files ? files.findAll { it.name.contains('_short_peptides.fasta') } : []
            tuple(sample_id, short_peptides)
        }
        .mix(
            nonframeshift_substitution_results.nfssub_fasta_files
                .map { sample_id, files -> 
                    def short_peptides = files ? files.findAll { it.name.contains('_short_peptides.fasta') } : []
                    tuple(sample_id, short_peptides)
                }
        )
        .mix(
            stoploss_mutation_results.stoploss_fasta_files
                .map { sample_id, files -> 
                    def short_peptides = files ? files.findAll { it.name.contains('_short_peptides.fasta') } : []
                    tuple(sample_id, short_peptides)
                }
        )
        .mix(
            snv_results.snv_fasta_files
                .map { sample_id, files -> 
                    def short_peptides = files ? files.findAll { it.name.contains('_short_peptides.fasta') } : []
                    tuple(sample_id, short_peptides)
                }
        )
        .mix(
            frameshift_del_results.fsdel_fasta_files
                .map { sample_id, files -> 
                    def short_peptides = files ? files.findAll { it.name.contains('_short_peptides.fasta') } : []
                    tuple(sample_id, short_peptides)
                }
        )
        .mix(
            nonframeshift_del_results.nfsdel_fasta_files
                .map { sample_id, files -> 
                    def short_peptides = files ? files.findAll { it.name.contains('_short_peptides.fasta') } : []
                    tuple(sample_id, short_peptides)
                }
        )
        .mix(
            nonframeshift_ins_results.nfsins_fasta_files
                .map { sample_id, files -> 
                    def short_peptides = files ? files.findAll { it.name.contains('_short_peptides.fasta') } : []
                    tuple(sample_id, short_peptides)
                }
        )
        .groupTuple()
    
    // 收集所有长肽段文件
    all_long_peptides = frameshift_insertion_results.fsins_fasta_files
        .map { sample_id, files -> 
            def long_peptides = files ? files.findAll { it.name.contains('_long_peptides.fasta') } : []
            tuple(sample_id, long_peptides)
        }
        .mix(
            nonframeshift_substitution_results.nfssub_fasta_files
                .map { sample_id, files -> 
                    def long_peptides = files ? files.findAll { it.name.contains('_long_peptides.fasta') } : []
                    tuple(sample_id, long_peptides)
                }
        )
        .mix(
            stoploss_mutation_results.stoploss_fasta_files
                .map { sample_id, files -> 
                    def long_peptides = files ? files.findAll { it.name.contains('_long_peptides.fasta') } : []
                    tuple(sample_id, long_peptides)
                }
        )
        .mix(
            snv_results.snv_fasta_files
                .map { sample_id, files -> 
                    def long_peptides = files ? files.findAll { it.name.contains('_long_peptides.fasta') } : []
                    tuple(sample_id, long_peptides)
                }
        )
        .mix(
            frameshift_del_results.fsdel_fasta_files
                .map { sample_id, files -> 
                    def long_peptides = files ? files.findAll { it.name.contains('_long_peptides.fasta') } : []
                    tuple(sample_id, long_peptides)
                }
        )
        .mix(
            nonframeshift_del_results.nfsdel_fasta_files
                .map { sample_id, files -> 
                    def long_peptides = files ? files.findAll { it.name.contains('_long_peptides.fasta') } : []
                    tuple(sample_id, long_peptides)
                }
        )
        .mix(
            nonframeshift_ins_results.nfsins_fasta_files
                .map { sample_id, files -> 
                    def long_peptides = files ? files.findAll { it.name.contains('_long_peptides.fasta') } : []
                    tuple(sample_id, long_peptides)
                }
        )
        .groupTuple()
    
    // 提取frameshift deletion junction肽段
    fsdel_junction_peptides = EXTRACT_FSDEL_JUNCTION_PEPTIDES(
        frameshift_del_results.fsdel_var_sequence
    )
    
    // 提取SNV肽段用于质谱分析
    snv_ms_peptides = EXTRACT_MS_PEPTIDES(
        snv_results.snv_var_sequence.map { sample_id, files -> tuple(sample_id, files, "snv") }
    )
    
    // 提取Nonframeshift Deletion肽段用于质谱分析
    nfsdel_ms_peptides = EXTRACT_MS_PEPTIDES(
        nonframeshift_del_results.nfsdel_var_sequence.map { sample_id, files -> tuple(sample_id, files, "nfsdel") }
    )
    
    // 提取Nonframeshift Insertion肽段用于质谱分析
    nfsins_ms_peptides = EXTRACT_MS_PEPTIDES(
        nonframeshift_ins_results.nfsins_var_sequence.map { sample_id, files -> tuple(sample_id, files, "nfsins") }
    )
    
    // 提取Frameshift Insertion肽段用于质谱分析
    fsins_ms_peptides = EXTRACT_MS_PEPTIDES(
        frameshift_insertion_results.fsins_var_sequence.map { sample_id, files -> tuple(sample_id, files, "fsins") }
    )
    
    // 提取Nonframeshift Substitution肽段用于质谱分析
    nfssub_ms_peptides = EXTRACT_MS_PEPTIDES(
        nonframeshift_substitution_results.nfssub_var_sequence.map { sample_id, files -> tuple(sample_id, files, "nfssub") }
    )
    
    // 提取Stoploss Mutation肽段用于质谱分析
    stoploss_ms_peptides = EXTRACT_MS_PEPTIDES(
        stoploss_mutation_results.stoploss_var_sequence.map { sample_id, files -> tuple(sample_id, files, "stoploss") }
    )

    // 整合所有基因摘要
    all_integrated_results = INTEGRATE_RESULTS(
        all_gene_summaries,
        all_short_peptides,
        all_long_peptides
    )

    emit:
    all_integrated_results = all_integrated_results.integrated_summary
    all_short_peptides = all_integrated_results.integrated_short_peptides
    all_long_peptides = all_integrated_results.integrated_long_peptides
    
    // 用于质谱分析的肽段文件
    fs_del_junction_1 = fsdel_junction_peptides.fsdel_junction_peptides.map { sample_id, file1, file2 -> tuple(sample_id, file1) }
    fs_del_junction_2 = fsdel_junction_peptides.fsdel_junction_peptides.map { sample_id, file1, file2 -> tuple(sample_id, file2) }
    snv_peptides = snv_ms_peptides.ms_peptides
    nfsdel_peptides = nfsdel_ms_peptides.ms_peptides
    nfsins_peptides = nfsins_ms_peptides.ms_peptides
    fsins_peptides = fsins_ms_peptides.ms_peptides
    nfssub_peptides = nfssub_ms_peptides.ms_peptides
    stoploss_peptides = stoploss_ms_peptides.ms_peptides
    
    // 保留原始输出以保持兼容性
    variant_sequence_files = all_variant_sequences
    
    // SNV分析结果
    snv_summary = snv_results.snv_gene_summary       // SNV基因摘要统计
    snv_var_sequence = snv_results.snv_var_sequence  // 变异蛋白序列FASTA文件
    snv_fasta_files = snv_results.snv_fasta_files    // 短肽段FASTA文件(8-11aa)
    snv_csv_files = snv_results.snv_csv_files        // SNV详细结果CSV文件
    snv_txt_files = snv_results.snv_txt_files        // SNV文本结果文件
    
    // Frameshift Deletion分析结果
    fsdel_summary = frameshift_del_results.fsdel_gene_summary       // Frameshift删除基因摘要统计
    fsdel_var_sequence = frameshift_del_results.fsdel_var_sequence  // Frameshift删除变异蛋白序列FASTA文件
    fsdel_fasta_files = frameshift_del_results.fsdel_fasta_files    // Frameshift删除FASTA文件
    fsdel_csv_files = frameshift_del_results.fsdel_csv_files        // Frameshift删除详细结果CSV文件
    fsdel_txt_files = frameshift_del_results.fsdel_txt_files        // Frameshift删除文本结果文件
    
    // Nonframeshift Deletion分析结果
    nfsdel_summary = nonframeshift_del_results.nfsdel_gene_summary       // Nonframeshift删除基因摘要统计
    nfsdel_var_sequence = nonframeshift_del_results.nfsdel_var_sequence  // Nonframeshift删除变异蛋白序列FASTA文件
    nfsdel_fasta_files = nonframeshift_del_results.nfsdel_fasta_files    // Nonframeshift删除FASTA文件
    nfsdel_csv_files = nonframeshift_del_results.nfsdel_csv_files        // Nonframeshift删除详细结果CSV文件
    nfsdel_txt_files = nonframeshift_del_results.nfsdel_txt_files        // Nonframeshift删除文本结果文件
    
    // Nonframeshift Insertion分析结果
    nfsins_summary = nonframeshift_ins_results.nfsins_gene_summary       // Nonframeshift插入基因摘要统计
    nfsins_var_sequence = nonframeshift_ins_results.nfsins_var_sequence  // Nonframeshift插入变异蛋白序列FASTA文件
    nfsins_fasta_files = nonframeshift_ins_results.nfsins_fasta_files    // Nonframeshift插入FASTA文件
    nfsins_csv_files = nonframeshift_ins_results.nfsins_csv_files        // Nonframeshift插入详细结果CSV文件
    nfsins_txt_files = nonframeshift_ins_results.nfsins_txt_files        // Nonframeshift插入文本结果文件
    
    // Frameshift Insertion分析结果
    fsins_summary = frameshift_insertion_results.fsins_gene_summary       // Frameshift插入基因摘要统计
    fsins_var_sequence = frameshift_insertion_results.fsins_var_sequence  // Frameshift插入变异蛋白序列FASTA文件
    fsins_fasta_files = frameshift_insertion_results.fsins_fasta_files    // Frameshift插入FASTA文件
    fsins_csv_files = frameshift_insertion_results.fsins_csv_files        // Frameshift插入详细结果CSV文件
    fsins_txt_files = frameshift_insertion_results.fsins_txt_files        // Frameshift插入文本结果文件
    
    // Nonframeshift Substitution分析结果
    nfssub_summary = nonframeshift_substitution_results.nfssub_gene_summary       // Nonframeshift替换基因摘要统计
    nfssub_var_sequence = nonframeshift_substitution_results.nfssub_var_sequence  // Nonframeshift替换变异蛋白序列FASTA文件
    nfssub_fasta_files = nonframeshift_substitution_results.nfssub_fasta_files    // Nonframeshift替换FASTA文件
    nfssub_csv_files = nonframeshift_substitution_results.nfssub_csv_files        // Nonframeshift替换详细结果CSV文件
    nfssub_txt_files = nonframeshift_substitution_results.nfssub_txt_files        // Nonframeshift替换文本结果文件
    
    // Stoploss Mutation分析结果
    stoploss_summary = stoploss_mutation_results.stoploss_gene_summary       // Stoploss突变基因摘要统计
    stoploss_var_sequence = stoploss_mutation_results.stoploss_var_sequence  // Stoploss突变变异蛋白序列FASTA文件
    stoploss_fasta_files = stoploss_mutation_results.stoploss_fasta_files    // Stoploss突变FASTA文件
    stoploss_csv_files = stoploss_mutation_results.stoploss_csv_files        // Stoploss突变详细结果CSV文件
    stoploss_txt_files = stoploss_mutation_results.stoploss_txt_files        // Stoploss突变文本结果文件
}
