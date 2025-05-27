// GATK MarkDuplicates process
process GATK_MARKDUPLICATES {
    tag "${sample_id}"
    label 'process_medium'
    
    container { "${params.gatk.container}" }
    
    publishDir {
        def base_sample_id = sample_id.replaceAll('_(dna|rna)_(normal|tumor)$', '')
        
        // Determine data type (DNA or RNA)
        def data_type = "dna"
        if (sample_id.contains("_rna_")) {
            data_type = "rna"
        }
        
        // Determine sample type (normal or tumor)
        def sample_type = "normal"
        if (sample_id.contains("_tumor")) {
            sample_type = "tumor"
        }
        
        return "${params.outdir}/${base_sample_id}/${data_type}/04_gatk/markduplicates"
    },
    mode: 'copy',
    saveAs: { filename ->
        if (filename.endsWith('.bam')) "bam/${filename}"
        else if (filename.endsWith('.txt')) "metrics/${filename}"
        else null
    }

    input:
    tuple val(sample_id), path(sorted_bam)

    output:
    tuple val(sample_id), 
          path("${sample_id}_marked_duplicates.bam"),
          emit: marked_bam
    path "${sample_id}_marked_dup_metrics.txt",
          emit: metrics

    script:
    def base_sample_id = sample_id.replaceAll('_(dna|rna)_(normal|tumor)$', '')
    
    // Determine data type (DNA or RNA)
    def data_type = "dna"
    if (sample_id.contains("_rna_")) {
        data_type = "rna"
    }
    
    // Determine sample type (normal or tumor)
    def sample_type = "normal"
    if (sample_id.contains("_tumor")) {
        sample_type = "tumor"
    }
    
    def output_dir = "${params.outdir}/${base_sample_id}/${data_type}/04_gatk/markduplicates"
    def marked_bam = "${output_dir}/bam/${sample_id}_marked_duplicates.bam"
    def metrics_file = "${output_dir}/metrics/${sample_id}_marked_dup_metrics.txt"
    def tmp_dir = "/tmp/${workflow.sessionId}/${sample_id}"

    """
    if [[ -f "${marked_bam}" && -f "${metrics_file}" ]]; then
        echo "MarkDuplicates output files already exist for ${sample_id}, creating symlinks..."
        ln -s "${marked_bam}" "${sample_id}_marked_duplicates.bam"
        ln -s "${metrics_file}" "${sample_id}_marked_dup_metrics.txt"
    else
        echo "Processing ${sample_id} with GATK MarkDuplicates..."
        mkdir -p ${tmp_dir}
        
        gatk MarkDuplicates \\
            -I ${sorted_bam} \\
            -O ${sample_id}_marked_duplicates.bam \\
            -M ${sample_id}_marked_dup_metrics.txt \\
            --TMP_DIR ${tmp_dir} \\
            --CREATE_INDEX true \\
            --VALIDATION_STRINGENCY SILENT
            
        rm -rf ${tmp_dir}
    fi
    """
}

// GATK AddOrReplaceReadGroups process
process GATK_ADD_OR_REPLACE_READ_GROUPS {
    tag "${sample_id}"
    label 'process_medium'
    
    container { "${params.gatk.container}" }
    
    publishDir {
        def base_sample_id = sample_id.replaceAll('_(dna|rna)_(normal|tumor)$', '')
        
        // Determine data type (DNA or RNA)
        def data_type = "dna"
        if (sample_id.contains("_rna_")) {
            data_type = "rna"
        }
        
        // Determine sample type (normal or tumor)
        def sample_type = "normal"
        if (sample_id.contains("_tumor")) {
            sample_type = "tumor"
        }
        
        return "${params.outdir}/${base_sample_id}/${data_type}/04_gatk/read_groups"
    },
    mode: 'copy',
    saveAs: { filename ->
        if (filename.endsWith('.bam')) "bam/${filename}"
        else null
    }

    input:
    tuple val(sample_id), path(marked_bam)

    output:
    tuple val(sample_id), 
          path("${sample_id}_with_read_groups.bam"),
          emit: bam_with_read_groups

    script:
    def base_sample_id = sample_id.replaceAll('_(dna|rna)_(normal|tumor)$', '')
    
    //  Determine data type (DNA or RNA)
    def data_type = "dna"
    if (sample_id.contains("_rna_")) {
        data_type = "rna"
    }
    
    // Determine sample type (normal or tumor)
    def sample_type = "normal"
    if (sample_id.contains("_tumor")) {
        sample_type = "tumor"
    }
    
    def output_dir = "${params.outdir}/${base_sample_id}/${data_type}/04_gatk/read_groups"
    def output_bam = "${output_dir}/bam/${sample_id}_with_read_groups.bam"
    def tmp_dir = "/tmp/${workflow.sessionId}/${sample_id}"

    """
    if [[ -f "${output_bam}" ]]; then
        echo "AddOrReplaceReadGroups output file already exists for ${sample_id}, creating symlink..."
        ln -s "${output_bam}" "${sample_id}_with_read_groups.bam"
    else
        echo "Processing ${sample_id} with GATK AddOrReplaceReadGroups..."
        mkdir -p ${tmp_dir}
        
        gatk AddOrReplaceReadGroups \\
            -I ${marked_bam} \\
            -O ${sample_id}_with_read_groups.bam \\
            -ID "${params.gatk.read_groups.id}" \\
            -LB "${params.gatk.read_groups.library}" \\
            -PL "${params.gatk.read_groups.platform}" \\
            -PU "${params.gatk.read_groups.platform_unit}" \\
            -SM ${sample_id} \\
            --TMP_DIR ${tmp_dir}
            
        rm -rf ${tmp_dir}
    fi
    """
}

// GATK BaseRecalibrator process
process GATK_BASE_RECALIBRATOR {
    tag "${sample_id}"
    label 'process_medium'
    
    container { "${params.gatk.container}" }
    
    publishDir {
        def base_sample_id = sample_id.replaceAll('_(dna|rna)_(normal|tumor)$', '')
        
        //  Determine data type (DNA or RNA)
        def data_type = "dna"
        if (sample_id.contains("_rna_")) {
            data_type = "rna"
        }
        
        // Determine sample type (normal or tumor)
        def sample_type = "normal"
        if (sample_id.contains("_tumor")) {
            sample_type = "tumor"
        }
        
        return "${params.outdir}/${base_sample_id}/${data_type}/04_gatk/base_recalibrator"
    },
    mode: 'copy',
    saveAs: { filename ->
        if (filename.endsWith('.table')) "recal/${filename}"
        else null
    }

    input:
    tuple val(sample_id), path(bam_with_read_groups)
    path reference_file
    
    output:
    tuple val(sample_id), 
          path("${sample_id}_recal_data.table"),
          emit: recal_table

    script:
    def base_sample_id = sample_id.replaceAll('_(dna|rna)_(normal|tumor)$', '')
    
    // Determine data type (DNA or RNA)
    def data_type = "dna"
    if (sample_id.contains("_rna_")) {
        data_type = "rna"
    }
    
    // Determine sample type (normal or tumor)
    def sample_type = "normal"
    if (sample_id.contains("_tumor")) {
        sample_type = "tumor"
    }
    
    def output_dir = "${params.outdir}/${base_sample_id}/${data_type}/04_gatk/base_recalibrator"
    def recal_table = "${output_dir}/recal/${sample_id}_recal_data.table"
    def tmp_dir = "/tmp/${workflow.sessionId}/${sample_id}"

    """
    # 获取参考基因组文件的路径      
    ref_path=\$(readlink -f ${reference_file})
    ref_dir=\$(dirname \$ref_path)
    ref_base=\$(basename \$ref_path)
    dict_base=\${ref_base%.*}  # 移除.fa或.fasta扩展名
    
    # 为参考基因组及其索引文件创建软链接
    for f in \${ref_dir}/\${ref_base}*; do
        ln -sf \$f ./\$(basename \$f)
    done
    
    # 检查并创建必要的索引文件
    if [[ ! -f "\${ref_base}.fai" ]]; then
        echo "Creating FASTA index file..."
        samtools faidx \$ref_base
    fi
    
    if [[ ! -f "\${dict_base}.dict" ]]; then
        echo "Creating sequence dictionary..."
        gatk CreateSequenceDictionary -R \$ref_base -O \${dict_base}.dict
    fi
    
    if [[ -f "${recal_table}" ]]; then
        echo "BaseRecalibrator output file already exists for ${sample_id}, creating symlink..."
        ln -s "${recal_table}" "${sample_id}_recal_data.table"
    else
        echo "Processing ${sample_id} with GATK BaseRecalibrator..."
        mkdir -p ${tmp_dir}
        
        gatk BaseRecalibrator \\
            -I ${bam_with_read_groups} \\
            -R \$ref_base \\
            --known-sites ${params.gatk.known_sites.dbsnp} \\
            --known-sites ${params.gatk.known_sites.mills_1000g} \\
            -O ${sample_id}_recal_data.table \\
            --tmp-dir ${tmp_dir} \\
            --java-options "-XX:ParallelGCThreads=40 -Xmx64g"
            
        rm -rf ${tmp_dir}
    fi
    """
}

// GATK ApplyBQSR process
process GATK_APPLY_BQSR {
    tag "${sample_id}"
    label 'process_medium'
    
    container { "${params.gatk.container}" }
    
    publishDir {
        def base_sample_id = sample_id.replaceAll('_(dna|rna)_(normal|tumor)$', '')
        
        // 判断数据类型（DNA或RNA）
        def data_type = "dna"
        if (sample_id.contains("_rna_")) {
            data_type = "rna"
        }
        
        // 判断样本类型（normal或tumor）
        def sample_type = "normal"
        if (sample_id.contains("_tumor")) {
            sample_type = "tumor"
        }
        
        return "${params.outdir}/${base_sample_id}/${data_type}/04_gatk/apply_bqsr"
    },
    mode: 'copy',
    saveAs: { filename ->
        if (filename.endsWith('.bam')) "bam/${filename}"
        else null
    }

    input:
    tuple val(sample_id), 
          path(bam_with_read_groups),
          path(recal_table)
    path reference_file

    output:
    tuple val(sample_id), 
          path("${sample_id}_recal.bam"),
          emit: recalibrated_bam

    script:
    def base_sample_id = sample_id.replaceAll('_(dna|rna)_(normal|tumor)$', '')
    
    // 判断数据类型（DNA或RNA）
    def data_type = "dna"
    if (sample_id.contains("_rna_")) {
        data_type = "rna"
    }
    
    // 判断样本类型（normal或tumor）
    def sample_type = "normal"
    if (sample_id.contains("_tumor")) {
        sample_type = "tumor"
    }
    
    def output_dir = "${params.outdir}/${base_sample_id}/${data_type}/04_gatk/apply_bqsr"
    def output_bam = "${output_dir}/bam/${sample_id}_recal.bam"
    def tmp_dir = "/tmp/${workflow.sessionId}/${sample_id}"

    """
    # 获取参考基因组文件的路径
    ref_path=\$(readlink -f ${reference_file})
    ref_dir=\$(dirname \$ref_path)
    ref_base=\$(basename \$ref_path)
    dict_base=\${ref_base%.*}  # 移除.fa或.fasta扩展名
    
    # 为参考基因组及其索引文件创建软链接
    for f in \${ref_dir}/\${ref_base}*; do
        ln -sf \$f ./\$(basename \$f)
    done
    
    # 检查并创建必要的索引文件
    if [[ ! -f "\${ref_base}.fai" ]]; then
        echo "Creating FASTA index file..."
        samtools faidx \$ref_base
    fi
    
    if [[ ! -f "\${dict_base}.dict" ]]; then
        echo "Creating sequence dictionary..."
        gatk CreateSequenceDictionary -R \$ref_base -O \${dict_base}.dict
    fi
    
    # 确保.dict文件存在于正确位置
    if [[ ! -f "\${ref_base}.dict" && -f "\${dict_base}.dict" ]]; then
        ln -sf "\${dict_base}.dict" "\${ref_base}.dict"
    fi
    
    if [[ -f "${output_bam}" ]]; then
        echo "ApplyBQSR output file already exists for ${sample_id}, creating symlink..."
        ln -s "${output_bam}" "${sample_id}_recal.bam"
    else
        echo "Processing ${sample_id} with GATK ApplyBQSR..."
        mkdir -p ${tmp_dir}
        
        gatk ApplyBQSR \\
            -R \$ref_base \\
            -I ${bam_with_read_groups} \\
            --bqsr-recal-file ${recal_table} \\
            -O ${sample_id}_recal.bam \\
            --tmp-dir ${tmp_dir} \\
            --java-options "-XX:ParallelGCThreads=40 -Xmx64g"
            
        rm -rf ${tmp_dir}
    fi
    """
}

// GATK Mutect2 Normal
process GATK_MUTECT2_NORMAL {
    tag "$sample_id"
    label 'process_high'
    
    container "$params.gatk.container"
    
    publishDir(
        path: {
            def base_sample_id = sample_id.replaceAll('_(dna|rna)_(normal|tumor)$', '')
            
            // 判断数据类型（DNA或RNA）
            def data_type = "dna"
            if (sample_id.contains("_rna_")) {
                data_type = "rna"
            }
            
            return "${params.outdir}/${base_sample_id}/${data_type}/06_mutect2/normal"
        },
        mode: 'copy',
        saveAs: { filename ->
            if (filename.endsWith('.vcf')) return "vcf/${filename}"
            else if (filename.endsWith('.vcf.idx')) return "vcf/${filename}"
            else return null
        }
    )
    
    input:
    tuple val(sample_id), path(normal_bam), path(normal_bam_index)
    path reference_file
    
    output:
    tuple val(sample_id),
          path {
              def base_sample_id = sample_id.replaceAll('_(dna|rna)_(normal|tumor)$', '')
              
              // 判断数据类型（DNA或RNA）
              def data_type = "dna"
              if (sample_id.contains("_rna_")) {
                  data_type = "rna"
              }
              
              return "${base_sample_id}_${data_type}_normal.vcf"
          },
          path {
              def base_sample_id = sample_id.replaceAll('_(dna|rna)_(normal|tumor)$', '')
              
              // 判断数据类型（DNA或RNA）
              def data_type = "dna"
              if (sample_id.contains("_rna_")) {
                  data_type = "rna"
              }
              
              return "${base_sample_id}_${data_type}_normal.vcf.idx"
          },
          emit: normal_vcf
    
    script:
    def base_sample_id = sample_id.replaceAll('_(dna|rna)_(normal|tumor)$', '')
    
    // 判断数据类型（DNA或RNA）
    def data_type = "dna"
    if (sample_id.contains("_rna_")) {
        data_type = "rna"
    }
    
    def output_vcf = "${base_sample_id}_${data_type}_normal.vcf"
    def tmp_dir = "/tmp/${workflow.sessionId}/${sample_id}"
    
    """
    # 获取参考基因组文件的路径
    ref_path=\$(readlink -f "${reference_file}")
    ref_dir=\$(dirname "\$ref_path")
    ref_base=\$(basename "\$ref_path")
    dict_base=\${ref_base%.*}  # 移除.fa或.fasta扩展名
    
    # 为参考基因组及其索引文件创建软链接
    for f in "\${ref_dir}/\${ref_base}"*; do
        ln -sf "\$f" ./\$(basename "\$f")
    done
    
    # 检查并创建必要的索引文件
    if [[ ! -f "\${ref_base}.fai" ]]; then
        echo "Creating FASTA index file..."
        samtools faidx "\$ref_base"
    fi
    
    if [[ ! -f "\${dict_base}.dict" ]]; then
        echo "Creating sequence dictionary..."
        gatk CreateSequenceDictionary -R "\$ref_base" -O "\${dict_base}.dict"
    fi
    
    # 确保.dict文件存在于正确位置
    if [[ ! -f "\${ref_base}.dict" && -f "\${dict_base}.dict" ]]; then
        ln -sf "\${dict_base}.dict" "\${ref_base}.dict"
    fi
    
    echo "Step 1: Calling mutations on normal sample..."
    mkdir -p "${tmp_dir}"
    
    gatk Mutect2 \\
        -R "\$ref_base" \\
        -I "${normal_bam}" \\
        -O "${output_vcf}" \\
        --tmp-dir "${tmp_dir}"
    
    rm -rf "${tmp_dir}"
    """
}

// GATK Mutect2 Paired
process GATK_MUTECT2_PAIRED {
    tag "${tumor_id}"
    label 'process_high'
    
    container "$params.gatk.container"
    
    // 失败时不要删除工作目录，保留所有中间文件
    errorStrategy { task.exitStatus == 143 ? 'retry' : 'finish' }
    maxRetries 3
    
    publishDir(
        path: {
            def base_sample_id = tumor_id.replaceAll('_(dna|rna)_(normal|tumor)$', '')
            
            // 判断数据类型（DNA或RNA）
            def data_type = "dna"
            if (tumor_id.contains("_rna_")) {
                data_type = "rna"
            }
            
            return "${params.outdir}/${base_sample_id}/${data_type}/06_mutect2/paired"
        },
        mode: 'copy',
        saveAs: { filename ->
            if (filename.endsWith('.vcf')) return "vcf/${filename}"
            else if (filename.endsWith('.vcf.idx')) return "vcf/${filename}"
            else if (filename.endsWith('.vcf.stats')) return "vcf/${filename}"
            else return null
        }
    )
    
    // 添加额外的发布目录，用于保存所有中间文件（在调试模式下）
    publishDir(
        path: {
            def base_sample_id = tumor_id.replaceAll('_(dna|rna)_(normal|tumor)$', '')
            def data_type = tumor_id.contains("_rna_") ? "rna" : "dna"
            String timestamp = new java.text.SimpleDateFormat("yyyyMMdd_HHmmss").format(new java.util.Date())
            return "${params.outdir}/${base_sample_id}/${data_type}/06_mutect2/paired/debug/${timestamp}"
        },
        mode: 'copy',
        pattern: '*',
        enabled: true  // 始终启用，确保任务失败时也保存文件
    )
    
    input:
    tuple val(tumor_id), path(tumor_bam), path(tumor_bam_index),
          val(normal_id), path(normal_bam), path(normal_bam_index),
          path(normal_vcf), path(normal_vcf_idx)
    path reference_file
    
    output:
    tuple val(tumor_id),
         path { 
            def base_sample_id = tumor_id.replaceAll('_(dna|rna)_(normal|tumor)$', '')
            
            // 判断数据类型（DNA或RNA）
            def data_type = "dna"
            if (tumor_id.contains("_rna_")) {
                data_type = "rna"
            }
            
            return "${base_sample_id}_${data_type}_paired.vcf"
         },
         path { 
            def base_sample_id = tumor_id.replaceAll('_(dna|rna)_(normal|tumor)$', '')
            
            // 判断数据类型（DNA或RNA）
            def data_type = "dna"
            if (tumor_id.contains("_rna_")) {
                data_type = "rna"
            }
            
            return "${base_sample_id}_${data_type}_paired.vcf.idx"
         },
         path { 
            def base_sample_id = tumor_id.replaceAll('_(dna|rna)_(normal|tumor)$', '')
            
            // 判断数据类型（DNA或RNA）
            def data_type = "dna"
            if (tumor_id.contains("_rna_")) {
                data_type = "rna"
            }
            
            return "${base_sample_id}_${data_type}_paired.vcf.stats"
         },
         emit: paired_vcf
    // 添加额外的调试输出通道，便于检查各种文件
    path "*", optional: true, emit: debug_files
    
    script:
    def base_sample_id = tumor_id.replaceAll('_(dna|rna)_(normal|tumor)$', '')
    
    // 判断数据类型（DNA或RNA）
    def data_type = "dna"
    if (tumor_id.contains("_rna_")) {
        data_type = "rna"
    }
    
    def output_vcf = "${base_sample_id}_${data_type}_paired.vcf"
    def tmp_dir = "${params.outdir}/tmp/${tumor_id}"
    
    """
    # 输出所有环境变量，便于调试
    echo "===== ENVIRONMENT VARIABLES ====="
    env | sort
    echo "================================="
    
    # 输出所有输入文件信息
    echo "===== INPUT FILES ====="
    echo "Tumor BAM: ${tumor_bam}"
    ls -la ${tumor_bam}*
    echo "Normal BAM: ${normal_bam}"
    ls -la ${normal_bam}*
    echo "Normal VCF: ${normal_vcf}"
    ls -la ${normal_vcf}*
    echo "Reference: ${reference_file}"
    ls -la ${reference_file}
    echo "=========================="
    
    # 获取参考基因组文件的路径
    ref_path=\$(readlink -f "${reference_file}")
    ref_dir=\$(dirname "\$ref_path")
    ref_base=\$(basename "\$ref_path")
    dict_base=\${ref_base%.*}  # 移除.fa或.fasta扩展名
    
    # 为参考基因组及其索引文件创建软链接
    for f in "\${ref_dir}/\${ref_base}"*; do
        ln -sf "\$f" ./\$(basename "\$f")
    done
    
    # 检查并创建必要的索引文件
    if [[ ! -f "\${ref_base}.fai" ]]; then
        echo "Creating FASTA index file..."
        samtools faidx "\$ref_base"
    fi
    
    if [[ ! -f "\${dict_base}.dict" ]]; then
        echo "Creating sequence dictionary..."
        gatk CreateSequenceDictionary -R "\$ref_base" -O "\${dict_base}.dict"
    fi
    
    # 确保.dict文件存在于正确位置
    if [[ ! -f "\${ref_base}.dict" && -f "\${dict_base}.dict" ]]; then
        ln -sf "\${dict_base}.dict" "\${ref_base}.dict"
    fi
    
    echo "Step 2: Calling mutations on tumor sample (paired with control)..."
    mkdir -p "${tmp_dir}"
    
    # 展示所有准备好的文件
    echo "===== PREPARED FILES ====="
    ls -la
    echo "=========================="
    
    # 保存完整的命令到文件中，便于复现
    echo "gatk Mutect2 \\
        -R \$ref_base \\
        -I ${tumor_bam} \\
        -tumor ${tumor_id} \\
        -I ${normal_bam} \\
        -normal ${normal_id} \\
        -pon ${normal_vcf} \\
        -O ${output_vcf} \\
        --tmp-dir ${tmp_dir}" > mutect2_command.txt
    
    gatk Mutect2 \\
        -R "\$ref_base" \\
        -I "${tumor_bam}" \\
        -tumor "${tumor_id}" \\
        -I "${normal_bam}" \\
        -normal "${normal_id}" \\
        -pon "${normal_vcf}" \\
        -O "${output_vcf}" \\
        --tmp-dir "${tmp_dir}"
    
    # 检查命令执行结果
    MUTECT2_EXIT=\$?
    echo "Mutect2 exit code: \$MUTECT2_EXIT"
    
    # 输出生成的文件
    echo "===== OUTPUT FILES ====="
    ls -la
    echo "========================"
    
    # 如果没有生成stats文件，创建一个空的
    if [[ ! -f "${output_vcf}.stats" ]]; then
        echo "Stats file not found, creating empty one"
        touch "${output_vcf}.stats"
    fi
    
    rm -rf "${tmp_dir}"
    
    # 如果Mutect2执行失败，返回错误码
    if [[ \$MUTECT2_EXIT -ne 0 ]]; then
        echo "ERROR: Mutect2 failed with exit code \$MUTECT2_EXIT"
        exit \$MUTECT2_EXIT
    fi
    """
}

// GATK FilterMutectCalls
process GATK_FILTER_MUTECT_CALLS {
    tag "$sample_id"
    label 'process_medium'
    
    container "$params.gatk.container"
    
    publishDir(
        path: {
            def base_sample_id = sample_id.replaceAll('_(dna|rna)_(normal|tumor)$', '')
            
            // 判断数据类型（DNA或RNA）
            def data_type = "dna"
            if (sample_id.contains("_rna_")) {
                data_type = "rna"
            }
            
            return "${params.outdir}/${base_sample_id}/${data_type}/06_mutect2/filtered"
        },
        mode: 'copy',
        saveAs: { filename ->
            if (filename.endsWith('.vcf')) return "vcf/${filename}"
            else if (filename.endsWith('.vcf.idx')) return "vcf/${filename}"
            else if (filename.endsWith('.vcf.stats')) return "vcf/${filename}"
            else return null
        }
    )
    
    input:
    tuple val(sample_id), path(vcf), path(vcf_idx), path(vcf_stats)
    path reference_file
    
    output:
    tuple val(sample_id),
          path {
              def base_sample_id = sample_id.replaceAll('_(dna|rna)_(normal|tumor)$', '')
              
              // 判断数据类型（DNA或RNA）
              def data_type = "dna"
              if (sample_id.contains("_rna_")) {
                  data_type = "rna"
              }
              
              return "${base_sample_id}_${data_type}_filtered.vcf"
          },
          path {
              def base_sample_id = sample_id.replaceAll('_(dna|rna)_(normal|tumor)$', '')
              
              // 判断数据类型（DNA或RNA）
              def data_type = "dna"
              if (sample_id.contains("_rna_")) {
                  data_type = "rna"
              }
              
              return "${base_sample_id}_${data_type}_filtered.vcf.idx"
          },
          emit: filtered_vcf
    
    script:
    def base_sample_id = sample_id.replaceAll('_(dna|rna)_(normal|tumor)$', '')
    
    // 判断数据类型（DNA或RNA）
    def data_type = "dna"
    if (sample_id.contains("_rna_")) {
        data_type = "rna"
    }
    
    def output_vcf = "${base_sample_id}_${data_type}_filtered.vcf"
    def tmp_dir = "/tmp/${workflow.sessionId}/${sample_id}"
    
    """
    mkdir -p ${tmp_dir}
    
    # 获取参考基因组文件的路径
    ref_path=\$(readlink -f "${reference_file}")
    ref_dir=\$(dirname "\$ref_path")
    ref_base=\$(basename "\$ref_path")
    dict_base=\${ref_base%.*}  # 移除.fa或.fasta扩展名
    
    # 为参考基因组及其索引文件创建软链接
    for f in "\${ref_dir}/\${ref_base}"*; do
        ln -sf "\$f" ./\$(basename "\$f")
    done
    
    # 检查并创建必要的索引文件
    if [[ ! -f "\${ref_base}.fai" ]]; then
        echo "Creating FASTA index file..."
        samtools faidx "\$ref_base"
    fi
    
    if [[ ! -f "\${dict_base}.dict" ]]; then
        echo "Creating sequence dictionary..."
        gatk CreateSequenceDictionary -R "\$ref_base" -O "\${dict_base}.dict"
    fi
    
    # 确保.dict文件存在于正确位置
    if [[ ! -f "\${ref_base}.dict" && -f "\${dict_base}.dict" ]]; then
        ln -sf "\${dict_base}.dict" "\${ref_base}.dict"
    fi
    
    # 检查stats文件是否为空
    stats_size=\$(stat -c%s "${vcf_stats}" 2>/dev/null || echo "0")
    stats_param=""
    if [[ -f "${vcf_stats}" && "\$stats_size" -gt 0 ]]; then
        stats_param="--stats ${vcf_stats}"
        echo "Using stats file: ${vcf_stats}"
    else
        stats_param="--dummy-stats"
        echo "Stats file not found or empty, using --dummy-stats"
    fi
    
    gatk FilterMutectCalls \\
        -R "\$ref_base" \\
        -V ${vcf} \\
        -O ${output_vcf} \\
        \$stats_param \\
        --tmp-dir ${tmp_dir} \\
        --java-options "-XX:ParallelGCThreads=40 -Xmx64g"

    rm -rf ${tmp_dir}
    """
}

// BCFtools Filter PASS
process BCFTOOLS_FILTER_PASS {
    tag "$sample_id"
    label 'process_low'
    
    container "$params.bcftools.container"
    
    publishDir(
        path: {
            def base_sample_id = sample_id.replaceAll('_(dna|rna)_(normal|tumor)$', '')
            
            // 判断数据类型（DNA或RNA）
            def data_type = "dna"
            if (sample_id.contains("_rna_")) {
                data_type = "rna"
            }
            
            return "${params.outdir}/${base_sample_id}/${data_type}/06_mutect2/filtered_pass"
        },
        mode: 'copy',
        saveAs: { filename ->
            if (filename.endsWith('.vcf')) return "vcf/${filename}"
            else return null
        }
    )
    
    input:
    tuple val(sample_id), path(vcf), path(vcf_idx)
    
    output:
    tuple val(sample_id),
          path {
              def base_sample_id = sample_id.replaceAll('_(dna|rna)_(normal|tumor)$', '')
              
              // 判断数据类型（DNA或RNA）
              def data_type = "dna"
              if (sample_id.contains("_rna_")) {
                  data_type = "rna"
              }
              
              return "${base_sample_id}_${data_type}.PASS.vcf"
          },
          emit: pass_filtered_vcf
    
    script:
    def base_sample_id = sample_id.replaceAll('_(dna|rna)_(normal|tumor)$', '')
    
    // 判断数据类型（DNA或RNA）
    def data_type = "dna"
    if (sample_id.contains("_rna_")) {
        data_type = "rna"
    }
    
    def output_vcf = "${base_sample_id}_${data_type}.PASS.vcf"
    
    """
    bcftools view -f PASS -O v \\
        -o ${output_vcf} \\
        "${vcf}"
    """
}