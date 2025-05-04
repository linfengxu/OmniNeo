#!/usr/bin/env nextflow

// Workflow to combine DNA and RNA peptide sequences for mass spectrometry analysis

// Import processes
process COMBINE_MS_DATABASE {
    tag "$sample_id"
    label 'process_medium'
    container "${params.combine_ms.container}"
    cpus params.combine_ms.threads
    
    publishDir "${params.combine_ms.ms_outdir}/${sample_id}", mode: 'copy'
    
    input:
    tuple val(sample_id), 
          path(dna_fs_del_junction_peptides_1), 
          path(dna_fs_del_junction_peptides_2),
          path(dna_snv_peptides),
          path(dna_nfsdel_peptides),
          path(dna_nfsins_peptides),
          path(dna_fsins_peptides),
          path(dna_nfssub_peptides),
          path(dna_stoploss_peptides),
          path(rna_noncoding_long_peptides),
          path(rna_noncoding_short_peptides),
          path(rna_fusion_peptides)
    path(crap_fasta)
    path(uniprot_fasta)
    
    output:
    tuple val(sample_id), path("${sample_id}_combined_ms_database.fasta"), emit: combined_database
    tuple val(sample_id), path("${sample_id}_ms_database_stats.txt"), emit: database_stats
    
    script:
    """
    # 创建文件列表
    echo "构建质谱数据库文件列表:" > file_list.txt
    
    # 添加DNA变异序列文件（如果存在并且非空）
    echo "DNA变异序列文件:" >> file_list.txt
    
    # 检查frameshift deletion肽段文件1
    if [ -s "${dna_fs_del_junction_peptides_1}" ]; then
        echo "${dna_fs_del_junction_peptides_1}" >> file_list.txt
        echo "  - \$(basename ${dna_fs_del_junction_peptides_1}) [\$(grep -c '>' ${dna_fs_del_junction_peptides_1})个序列]" >> file_list.txt
    fi
    
    # 检查frameshift deletion肽段文件2
    if [ -s "${dna_fs_del_junction_peptides_2}" ]; then
        echo "${dna_fs_del_junction_peptides_2}" >> file_list.txt
        echo "  - \$(basename ${dna_fs_del_junction_peptides_2}) [\$(grep -c '>' ${dna_fs_del_junction_peptides_2})个序列]" >> file_list.txt
    fi
    
    # 检查SNV肽段文件
    if [ -s "${dna_snv_peptides}" ]; then
        echo "${dna_snv_peptides}" >> file_list.txt
        echo "  - \$(basename ${dna_snv_peptides}) [\$(grep -c '>' ${dna_snv_peptides})个序列]" >> file_list.txt
    fi
    
    # 检查Non-frameshift deletion肽段文件
    if [ -s "${dna_nfsdel_peptides}" ]; then
        echo "${dna_nfsdel_peptides}" >> file_list.txt
        echo "  - \$(basename ${dna_nfsdel_peptides}) [\$(grep -c '>' ${dna_nfsdel_peptides})个序列]" >> file_list.txt
    fi
    
    # 检查Non-frameshift insertion肽段文件
    if [ -s "${dna_nfsins_peptides}" ]; then
        echo "${dna_nfsins_peptides}" >> file_list.txt
        echo "  - \$(basename ${dna_nfsins_peptides}) [\$(grep -c '>' ${dna_nfsins_peptides})个序列]" >> file_list.txt
    fi
    
    # 检查Frameshift insertion肽段文件
    if [ -s "${dna_fsins_peptides}" ]; then
        echo "${dna_fsins_peptides}" >> file_list.txt
        echo "  - \$(basename ${dna_fsins_peptides}) [\$(grep -c '>' ${dna_fsins_peptides})个序列]" >> file_list.txt
    fi
    
    # 检查Non-frameshift substitution肽段文件
    if [ -s "${dna_nfssub_peptides}" ]; then
        echo "${dna_nfssub_peptides}" >> file_list.txt
        echo "  - \$(basename ${dna_nfssub_peptides}) [\$(grep -c '>' ${dna_nfssub_peptides})个序列]" >> file_list.txt
    fi
    
    # 检查Stoploss mutation肽段文件
    if [ -s "${dna_stoploss_peptides}" ]; then
        echo "${dna_stoploss_peptides}" >> file_list.txt
        echo "  - \$(basename ${dna_stoploss_peptides}) [\$(grep -c '>' ${dna_stoploss_peptides})个序列]" >> file_list.txt
    fi
    
    # 添加RNA非编码长肽段（如果存在并且非空）
    if [ -s "${rna_noncoding_long_peptides}" ]; then
        echo "RNA非编码长肽段文件:" >> file_list.txt
        echo "${rna_noncoding_long_peptides}" >> file_list.txt
        echo "  - \$(basename ${rna_noncoding_long_peptides}) [\$(grep -c '>' ${rna_noncoding_long_peptides})个序列]" >> file_list.txt
    fi
    
    # 添加RNA非编码短肽段（如果存在并且非空）
    if [ -s "${rna_noncoding_short_peptides}" ]; then
        echo "RNA非编码短肽段文件:" >> file_list.txt
        echo "${rna_noncoding_short_peptides}" >> file_list.txt
        echo "  - \$(basename ${rna_noncoding_short_peptides}) [\$(grep -c '>' ${rna_noncoding_short_peptides})个序列]" >> file_list.txt
    fi
    
    # 添加RNA融合肽段（如果存在并且非空）
    if [ -s "${rna_fusion_peptides}" ]; then
        echo "RNA融合肽段文件:" >> file_list.txt
        echo "${rna_fusion_peptides}" >> file_list.txt
        echo "  - \$(basename ${rna_fusion_peptides}) [\$(grep -c '>' ${rna_fusion_peptides})个序列]" >> file_list.txt
    fi
    
    # 添加参考数据库
    echo "参考数据库文件:" >> file_list.txt
    echo "${crap_fasta}" >> file_list.txt
    echo "  - \$(basename ${crap_fasta}) [\$(grep -c '>' ${crap_fasta})个序列]" >> file_list.txt
    echo "${uniprot_fasta}" >> file_list.txt
    echo "  - \$(basename ${uniprot_fasta}) [\$(grep -c '>' ${uniprot_fasta})个序列]" >> file_list.txt
    
    # 创建一个临时输入文件列表，只包含实际文件路径
    grep -v "^  - " file_list.txt | grep -v ":" | grep -v "^$" > input_files.txt
    
    # 合并所有文件到一个数据库
    cat \$(cat input_files.txt) > ${sample_id}_combined_ms_database.fasta
    
    # 生成统计信息
    echo "质谱数据库统计 - 样本: ${sample_id}" > ${sample_id}_ms_database_stats.txt
    echo "======================================================" >> ${sample_id}_ms_database_stats.txt
    echo "总序列数: \$(grep -c '>' ${sample_id}_combined_ms_database.fasta)" >> ${sample_id}_ms_database_stats.txt
    echo "" >> ${sample_id}_ms_database_stats.txt
    
    # 添加组成文件信息
    echo "组成文件:" >> ${sample_id}_ms_database_stats.txt
    grep "^  - " file_list.txt >> ${sample_id}_ms_database_stats.txt
    
    # 添加附加信息
    echo "" >> ${sample_id}_ms_database_stats.txt
    echo "数据库文件大小: \$(du -h ${sample_id}_combined_ms_database.fasta | cut -f1)" >> ${sample_id}_ms_database_stats.txt
    echo "创建时间: \$(date)" >> ${sample_id}_ms_database_stats.txt
    """
}

// 创建默认空文件的进程
process CREATE_EMPTY_FILES {
    tag "$sample_id"
    container "${params.combine_ms.container}"
    
    input:
    val(sample_id)
    
    output:
    tuple val(sample_id), 
          path("${sample_id}_empty_fs_del_junction_1.fasta"), 
          path("${sample_id}_empty_fs_del_junction_2.fasta"),
          path("${sample_id}_empty_snv.fasta"),
          path("${sample_id}_empty_nfsdel.fasta"),
          path("${sample_id}_empty_nfsins.fasta"),
          path("${sample_id}_empty_fsins.fasta"),
          path("${sample_id}_empty_nfssub.fasta"),
          path("${sample_id}_empty_stoploss.fasta"),
          emit: empty_dna_files
    tuple val(sample_id), 
          path("${sample_id}_empty_noncoding_long.fasta"),
          path("${sample_id}_empty_noncoding_short.fasta"),
          emit: empty_rna_noncoding_files
    tuple val(sample_id), 
          path("${sample_id}_empty_fusion.fasta"),
          emit: empty_rna_fusion
    
    script:
    """
    # 创建空的DNA frameshift deletion肽段文件1
    echo ">Empty_dna_fs_del_1_${sample_id}" > ${sample_id}_empty_fs_del_junction_1.fasta
    echo "EMPTYPEPTIDE" >> ${sample_id}_empty_fs_del_junction_1.fasta
    
    # 创建空的DNA frameshift deletion肽段文件2
    echo ">Empty_dna_fs_del_2_${sample_id}" > ${sample_id}_empty_fs_del_junction_2.fasta
    echo "EMPTYPEPTIDE" >> ${sample_id}_empty_fs_del_junction_2.fasta
    
    # 创建空的DNA SNV肽段文件
    echo ">Empty_dna_snv_${sample_id}" > ${sample_id}_empty_snv.fasta
    echo "EMPTYPEPTIDE" >> ${sample_id}_empty_snv.fasta
    
    # 创建空的DNA Non-frameshift deletion肽段文件
    echo ">Empty_dna_nfsdel_${sample_id}" > ${sample_id}_empty_nfsdel.fasta
    echo "EMPTYPEPTIDE" >> ${sample_id}_empty_nfsdel.fasta
    
    # 创建空的DNA Non-frameshift insertion肽段文件
    echo ">Empty_dna_nfsins_${sample_id}" > ${sample_id}_empty_nfsins.fasta
    echo "EMPTYPEPTIDE" >> ${sample_id}_empty_nfsins.fasta
    
    # 创建空的DNA Frameshift insertion肽段文件
    echo ">Empty_dna_fsins_${sample_id}" > ${sample_id}_empty_fsins.fasta
    echo "EMPTYPEPTIDE" >> ${sample_id}_empty_fsins.fasta
    
    # 创建空的DNA Non-frameshift substitution肽段文件
    echo ">Empty_dna_nfssub_${sample_id}" > ${sample_id}_empty_nfssub.fasta
    echo "EMPTYPEPTIDE" >> ${sample_id}_empty_nfssub.fasta
    
    # 创建空的DNA Stoploss mutation肽段文件
    echo ">Empty_dna_stoploss_${sample_id}" > ${sample_id}_empty_stoploss.fasta
    echo "EMPTYPEPTIDE" >> ${sample_id}_empty_stoploss.fasta
    
    # 创建空的RNA非编码长肽段文件
    echo ">Empty_noncoding_long_${sample_id}" > ${sample_id}_empty_noncoding_long.fasta
    echo "EMPTYPEPTIDE" >> ${sample_id}_empty_noncoding_long.fasta
    
    # 创建空的RNA非编码短肽段文件
    echo ">Empty_noncoding_short_${sample_id}" > ${sample_id}_empty_noncoding_short.fasta
    echo "EMPTYPEPTIDE" >> ${sample_id}_empty_noncoding_short.fasta
    
    # 创建空的RNA融合肽段文件
    echo ">Empty_fusion_${sample_id}" > ${sample_id}_empty_fusion.fasta
    echo "EMPTYPEPTIDE" >> ${sample_id}_empty_fusion.fasta
    """
}

// Define the main workflow
workflow COMBINE_MS {
    take:
    dna_fs_del_1_ch      // Channel with [sample_id, fs_del_junction_peptides_1.fasta]
    dna_fs_del_2_ch      // Channel with [sample_id, fs_del_junction_peptides_2.fasta]
    dna_snv_ch           // Channel with [sample_id, snv_var_sequence.fasta]
    dna_nfsdel_ch        // Channel with [sample_id, nfsdel_var_sequence.fasta]
    dna_nfsins_ch        // Channel with [sample_id, nfsins_var_sequence.fasta]
    dna_fsins_ch         // Channel with [sample_id, fsins_var_sequence.fasta]
    dna_nfssub_ch        // Channel with [sample_id, nfssub_var_sequence.fasta]
    dna_stoploss_ch      // Channel with [sample_id, stoploss_var_sequence.fasta]
    rna_noncoding_long_ch   // Channel with [sample_id, noncoding_long_peptides.fasta]
    rna_noncoding_short_ch  // Channel with [sample_id, noncoding_short_peptides.fasta]
    rna_fusion_ch        // Channel with [sample_id, fusion_long_peptides.fasta]
    
    main:
    // Reference databases
    crap_fasta = Channel.fromPath(params.combine_ms.crap_fasta, checkIfExists: true)
        .ifEmpty { error "CRAP database file not found at: ${params.combine_ms.crap_fasta}" }
    
    uniprot_fasta = Channel.fromPath(params.combine_ms.uniprot_fasta, checkIfExists: true)
        .ifEmpty { error "UniProt database file not found at: ${params.combine_ms.uniprot_fasta}" }
    
    // 收集所有样本ID
    all_sample_ids = dna_fs_del_1_ch.map { it[0] }
        .mix(dna_fs_del_2_ch.map { it[0] })
        .mix(dna_snv_ch.map { it[0] })
        .mix(dna_nfsdel_ch.map { it[0] })
        .mix(dna_nfsins_ch.map { it[0] })
        .mix(dna_fsins_ch.map { it[0] })
        .mix(dna_nfssub_ch.map { it[0] })
        .mix(dna_stoploss_ch.map { it[0] })
        .mix(rna_noncoding_long_ch.map { it[0] })
        .mix(rna_noncoding_short_ch.map { it[0] })
        .mix(rna_fusion_ch.map { it[0] })
        .unique()
    
    // 为每个样本创建空文件以防缺失数据
    empty_files = CREATE_EMPTY_FILES(all_sample_ids)
    
    // 处理DNA fs_del_junction_peptides_1文件
    dna_fs_del_1_safe = dna_fs_del_1_ch
        .map { sample_id, file -> tuple(sample_id, file) }
        .mix(empty_files.empty_dna_files.map { sample_id, file1, file2, file3, file4, file5, file6, file7, file8 -> tuple(sample_id, file1) })
        .groupTuple()
        .map { sample_id, files -> 
            if (files.size() > 1) {
                def nonEmptyFiles = files.findAll { !it.name.contains("empty") }
                return tuple(sample_id, nonEmptyFiles ? nonEmptyFiles[0] : files[0])
            } else {
                return tuple(sample_id, files[0])
            }
        }
    
    // 处理DNA fs_del_junction_peptides_2文件
    dna_fs_del_2_safe = dna_fs_del_2_ch
        .map { sample_id, file -> tuple(sample_id, file) }
        .mix(empty_files.empty_dna_files.map { sample_id, file1, file2, file3, file4, file5, file6, file7, file8 -> tuple(sample_id, file2) })
        .groupTuple()
        .map { sample_id, files -> 
            if (files.size() > 1) {
                def nonEmptyFiles = files.findAll { !it.name.contains("empty") }
                return tuple(sample_id, nonEmptyFiles ? nonEmptyFiles[0] : files[0])
            } else {
                return tuple(sample_id, files[0])
            }
        }
    
    // 处理DNA snv肽段文件
    dna_snv_safe = dna_snv_ch
        .map { sample_id, file -> tuple(sample_id, file) }
        .mix(empty_files.empty_dna_files.map { sample_id, file1, file2, file3, file4, file5, file6, file7, file8 -> tuple(sample_id, file3) })
        .groupTuple()
        .map { sample_id, files -> 
            if (files.size() > 1) {
                def nonEmptyFiles = files.findAll { !it.name.contains("empty") }
                return tuple(sample_id, nonEmptyFiles ? nonEmptyFiles[0] : files[0])
            } else {
                return tuple(sample_id, files[0])
            }
        }
    
    // 处理DNA nfsdel肽段文件
    dna_nfsdel_safe = dna_nfsdel_ch
        .map { sample_id, file -> tuple(sample_id, file) }
        .mix(empty_files.empty_dna_files.map { sample_id, file1, file2, file3, file4, file5, file6, file7, file8 -> tuple(sample_id, file4) })
        .groupTuple()
        .map { sample_id, files -> 
            if (files.size() > 1) {
                def nonEmptyFiles = files.findAll { !it.name.contains("empty") }
                return tuple(sample_id, nonEmptyFiles ? nonEmptyFiles[0] : files[0])
            } else {
                return tuple(sample_id, files[0])
            }
        }
    
    // 处理DNA nfsins肽段文件
    dna_nfsins_safe = dna_nfsins_ch
        .map { sample_id, file -> tuple(sample_id, file) }
        .mix(empty_files.empty_dna_files.map { sample_id, file1, file2, file3, file4, file5, file6, file7, file8 -> tuple(sample_id, file5) })
        .groupTuple()
        .map { sample_id, files -> 
            if (files.size() > 1) {
                def nonEmptyFiles = files.findAll { !it.name.contains("empty") }
                return tuple(sample_id, nonEmptyFiles ? nonEmptyFiles[0] : files[0])
            } else {
                return tuple(sample_id, files[0])
            }
        }
    
    // 处理DNA fsins肽段文件
    dna_fsins_safe = dna_fsins_ch
        .map { sample_id, file -> tuple(sample_id, file) }
        .mix(empty_files.empty_dna_files.map { sample_id, file1, file2, file3, file4, file5, file6, file7, file8 -> tuple(sample_id, file6) })
        .groupTuple()
        .map { sample_id, files -> 
            if (files.size() > 1) {
                def nonEmptyFiles = files.findAll { !it.name.contains("empty") }
                return tuple(sample_id, nonEmptyFiles ? nonEmptyFiles[0] : files[0])
            } else {
                return tuple(sample_id, files[0])
            }
        }
    
    // 处理DNA nfssub肽段文件
    dna_nfssub_safe = dna_nfssub_ch
        .map { sample_id, file -> tuple(sample_id, file) }
        .mix(empty_files.empty_dna_files.map { sample_id, file1, file2, file3, file4, file5, file6, file7, file8 -> tuple(sample_id, file7) })
        .groupTuple()
        .map { sample_id, files -> 
            if (files.size() > 1) {
                def nonEmptyFiles = files.findAll { !it.name.contains("empty") }
                return tuple(sample_id, nonEmptyFiles ? nonEmptyFiles[0] : files[0])
            } else {
                return tuple(sample_id, files[0])
            }
        }
    
    // 处理DNA stoploss肽段文件
    dna_stoploss_safe = dna_stoploss_ch
        .map { sample_id, file -> tuple(sample_id, file) }
        .mix(empty_files.empty_dna_files.map { sample_id, file1, file2, file3, file4, file5, file6, file7, file8 -> tuple(sample_id, file8) })
        .groupTuple()
        .map { sample_id, files -> 
            if (files.size() > 1) {
                def nonEmptyFiles = files.findAll { !it.name.contains("empty") }
                return tuple(sample_id, nonEmptyFiles ? nonEmptyFiles[0] : files[0])
            } else {
                return tuple(sample_id, files[0])
            }
        }
    
    // 处理RNA noncoding_long_peptides文件
    rna_noncoding_long_safe = rna_noncoding_long_ch
        .map { sample_id, file -> tuple(sample_id, file) }
        .mix(empty_files.empty_rna_noncoding_files.map { sample_id, file1, file2 -> tuple(sample_id, file1) })
        .groupTuple()
        .map { sample_id, files -> 
            if (files.size() > 1) {
                def nonEmptyFiles = files.findAll { !it.name.contains("empty") }
                return tuple(sample_id, nonEmptyFiles ? nonEmptyFiles[0] : files[0])
            } else {
                return tuple(sample_id, files[0])
            }
        }
    
    // 处理RNA noncoding_short_peptides文件
    rna_noncoding_short_safe = rna_noncoding_short_ch
        .map { sample_id, file -> tuple(sample_id, file) }
        .mix(empty_files.empty_rna_noncoding_files.map { sample_id, file1, file2 -> tuple(sample_id, file2) })
        .groupTuple()
        .map { sample_id, files -> 
            if (files.size() > 1) {
                def nonEmptyFiles = files.findAll { !it.name.contains("empty") }
                return tuple(sample_id, nonEmptyFiles ? nonEmptyFiles[0] : files[0])
            } else {
                return tuple(sample_id, files[0])
            }
        }
    
    // 处理RNA融合肽段文件
    rna_fusion_safe = rna_fusion_ch
        .map { sample_id, file -> tuple(sample_id, file) }
        .mix(empty_files.empty_rna_fusion)
        .groupTuple()
        .map { sample_id, files -> 
            if (files.size() > 1) {
                def nonEmptyFiles = files.findAll { !it.name.contains("empty") }
                return tuple(sample_id, nonEmptyFiles ? nonEmptyFiles[0] : files[0])
            } else {
                return tuple(sample_id, files[0])
            }
        }
    
    // 结合安全通道
    combined_input = dna_fs_del_1_safe
        .join(dna_fs_del_2_safe, by: 0)
        .join(dna_snv_safe, by: 0)
        .join(dna_nfsdel_safe, by: 0)
        .join(dna_nfsins_safe, by: 0)
        .join(dna_fsins_safe, by: 0)
        .join(dna_nfssub_safe, by: 0)
        .join(dna_stoploss_safe, by: 0)
        .join(rna_noncoding_long_safe, by: 0)
        .join(rna_noncoding_short_safe, by: 0)
        .join(rna_fusion_safe, by: 0)
    
    // 运行合并进程
    ms_database = COMBINE_MS_DATABASE(
        combined_input,
        crap_fasta,
        uniprot_fasta
    )
    
    // 输出处理完成信息
    ms_database.combined_database.subscribe { sample_id, ms_db ->
        log.info "创建质谱数据库完成，样本: ${sample_id}"
    }
    
    emit:
    ms_database = ms_database.combined_database
    database_stats = ms_database.database_stats
}
