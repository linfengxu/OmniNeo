#!/usr/bin/env nextflow

// Module for running MaxQuant analysis using MS files specified in the samplesheet

process GENERATE_MQPAR {
    tag "$sample_id"
    label 'process_high'
   // containerOptions "--bind ${params.tmp_dir}:/tmp --bind /mnt"
    publishDir "${params.maxquant.mqpar_dir}", mode: 'copy', pattern: "mqpar_${sample_id}.xml"
    
    input:
    tuple val(sample_id), path(ms_database), val(ms_files)
    
    output:
    tuple val(sample_id), path("ms_folder_path.txt"), emit: mqpar_config
    
    script:
    // Create command to copy MS files
    def ms_files_copy = ms_files.collect { ms_file -> "cp -L ${ms_file} ${params.maxquant.mqpar_tmp_dir}/ms" }.join("\n")
    
    """
    mkdir -p ${params.maxquant.mqpar_tmp_dir}/ms
    
    # Copy MS files to ms folder
    ${ms_files_copy}
    
    # Get absolute path of ms_database
    MS_DATABASE_ABS=\$(readlink -f ${ms_database})
    echo "MS Database absolute path: \${MS_DATABASE_ABS}"
    
    # Create MaxQuant configuration file
    singularity exec --bind /mnt:/mnt ${params.maxquant.container} \\
        dotnet /opt/MaxQuant_v2.6.7.0/bin/MaxQuantCmd.dll dummy.xml \\
        --create "${params.maxquant.mqpar_dir}/mqpar_${sample_id}.xml" \\
        --instrumentType TO \\
        --LCMSType ST \\
        --useLFQ \\
        --useMBR \\
        --numThreads 24 \\
        --pathFasta "\${MS_DATABASE_ABS}" \\
        --pathRawFileFolder "${params.maxquant.mqpar_tmp_dir}/ms"
    
    # Create output file
    echo "${params.maxquant.mqpar_tmp_dir}/ms" > ms_folder_path.txt
    
    # Run MaxQuant analysis
    singularity exec --bind /mnt:/mnt ${params.maxquant.container} \\
        dotnet /opt/MaxQuant_v2.6.7.0/bin/MaxQuantCmd.dll \\
        "${params.maxquant.mqpar_dir}/mqpar_${sample_id}.xml"
    """
}

workflow MAXQUANT_ANALYSIS {
    take:
    ms_database_ch        // [sample_id, ms_database.fasta]
    samplesheet_ch        // Path to samplesheet file
    
    main:
    // Extract MS file paths from samplesheet for each sample
    samplesheet_ch
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
    mqpar_configs = GENERATE_MQPAR(combined_input_ch)
    emit:
    mqpar_files = mqpar_configs.mqpar_config.map { it -> tuple(it[0], it[1]) }
    results = mqpar_configs.mqpar_config
}
