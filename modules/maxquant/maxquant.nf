#!/usr/bin/env nextflow

// Module for running MaxQuant analysis using MS files specified in the samplesheet

process GENERATE_MQPAR {
    tag "$sample_id"
    label 'process_low'
    container "${params.maxquant.python_container}"
    
    publishDir "${params.maxquant.mqpar_dir}", mode: 'copy', pattern: "mqpar_${sample_id}.xml"
    
    input:
    tuple val(sample_id), path(ms_database), val(ms_files)
    path(template_xml)
    
    output:
    tuple val(sample_id), path("mqpar_${sample_id}.xml"), path(ms_database), path("ms"), emit: mqpar_config
    
    script:
    // 创建拷贝MS文件的命令
    def ms_files_copy = ms_files.collect { ms_file -> "cp -L ${ms_file} ms/" }.join("\n")
    
    """
    # 创建MS文件夹
    mkdir -p ms
    
    # 复制MS文件到ms文件夹
    ${ms_files_copy}
    
    # 列出MS文件夹内容进行确认
    echo "MS files copied to ms folder:"
    ls -la ms/
    
    # 获取MS文件夹的绝对路径
    MS_DIR=\$(readlink -f ./ms)
    echo "MS folder absolute path: \$MS_DIR"
    
    # 生成MaxQuant配置
    python ${projectDir}/bin/gen_mqpar.py \
        ${template_xml} \
        \$MS_DIR \
        -o mqpar_${sample_id}.xml \
        -t ${params.maxquant.threads} \
        -f ${ms_database}
    """
}

process RUN_MAXQUANT {
    tag "$sample_id"
    label 'process_high'
    container "${params.maxquant.container}"
    time { params.maxquant.max_time }
    memory { params.maxquant.max_memory }
    cpus { params.maxquant.threads }
    
    publishDir "${params.maxquant.results_dir}/${sample_id}", mode: 'copy'
    
    input:
    tuple val(sample_id), path(mqpar_file), path(ms_database), path(ms_folder)
    
    output:
    tuple val(sample_id), path("maxquant_results"), emit: results
    
    script:
    """
    # 创建输出目录
    mkdir -p maxquant_results
    
    # 确保MS文件夹是绝对路径
    MS_DIR=\$(readlink -f ${ms_folder})
    
    # 显示MS文件夹内容
    echo "MS folder absolute path: \$MS_DIR"
    echo "MS folder contents:"
    ls -la \$MS_DIR
    
    # 获取mqpar文件的绝对路径
    MQPAR_PATH=\$(readlink -f ${mqpar_file})
    echo "Using mqpar file absolute path: \$MQPAR_PATH"
    
    # 查看mqpar文件内容中的文件路径
    echo "Checking paths in mqpar file:"
    grep -A 3 "<filePaths>" \$MQPAR_PATH || echo "No filePaths found in mqpar file"
    grep -A 1 "<fastaFilePath>" \$MQPAR_PATH || echo "No fastaFilePath found in mqpar file"
    
    # 运行MaxQuant，使用绝对路径
    echo "Running MaxQuant with configuration: \$MQPAR_PATH"
    mono /usr/local/bin/MaxQuantCmd.exe \$MQPAR_PATH
    
    # 复制结果到输出目录
    if [ -d "combined/txt" ]; then
        cp -r combined/txt maxquant_results/
    fi
    
    if [ -d "combined/proc" ]; then
        cp -r combined/proc maxquant_results/
    fi
    
    if [ -d "combined/andromeda" ]; then
        cp -r combined/andromeda maxquant_results/
    fi
    
    # 添加完成文件
    echo "MaxQuant analysis completed for sample ${sample_id} at \$(date)" > maxquant_results/completed.txt
    echo "MS database absolute path: \$(readlink -f ${ms_database})" >> maxquant_results/completed.txt
    echo "mqpar file absolute path: \$MQPAR_PATH" >> maxquant_results/completed.txt
    echo "MS folder absolute path: \$MS_DIR" >> maxquant_results/completed.txt
    echo "MS files used:" >> maxquant_results/completed.txt
    find \$MS_DIR -type f -name "*.raw" | while read file; do
        echo "  - \$(basename \$file) [\$file]" >> maxquant_results/completed.txt
    done
    """
}

workflow MAXQUANT_ANALYSIS {
    take:
    ms_database_ch        // [sample_id, ms_database.fasta]
    samplesheet_ch        // Path to samplesheet file
    
    main:
    // Get template XML
    template_xml = Channel.fromPath(params.maxquant.template_xml, checkIfExists: true)
        .ifEmpty { error "MaxQuant template XML file not found at: ${params.maxquant.template_xml}" }
    
    // Extract MS file paths from samplesheet for each sample
    Channel.fromPath(samplesheet_ch)
        .splitCsv(header:true, sep:'\t')
        .map { row -> 
            def sample_id = row.sampleID
            def ms_files = row.MSTumor ? row.MSTumor.split(',').collect { file(it.trim()) } : []
            if (ms_files.isEmpty()) {
                log.warn "No MS files found for sample ${sample_id}!"
            }
            return tuple(sample_id, ms_files)
        }
        .set { ms_files_ch }
    
    // Combine MS database with MS file paths
    ms_database_ch
        .join(ms_files_ch)
        .set { combined_input_ch }
    
    // Generate MaxQuant configuration
    mqpar_configs = GENERATE_MQPAR(
        combined_input_ch,
        template_xml
    )
    
    // Run MaxQuant if not skipped
    if (!params.maxquant.skip_run) {
        maxquant_results = RUN_MAXQUANT(
            mqpar_configs.mqpar_config
        )
    }
    
    emit:
    mqpar_files = mqpar_configs.mqpar_config.map { it -> tuple(it[0], it[1]) }
    results = params.maxquant.skip_run ? Channel.empty() : maxquant_results.results
}
