// HLA-II typing module, using HLAminer for HLA typing of RNA-seq tumor samples

// Step 1: Use BWA mem to align reads to HLA reference sequences
process BWA_ALIGN {
    tag "${sample_id}"
    label 'process_medium'
    
    container "${params.bwa.container}"
    
    containerOptions "--bind ${params.hlaminer.database_bwamem}:/data"
    
    input:
    tuple val(sample_id), path(read1), path(read2)
    
    output:
    tuple val(sample_id), 
          path("${sample_id}_vs_HLA.sam"),
          emit: align_result
    
    script:
    def hla_ref_path = "/data/HLA-I_II_CDS.fasta"
    def threads = params.bwa.threads
    
    """
    # Use BWA mem to align reads to HLA reference sequences
    bwa mem -a -t ${threads} ${hla_ref_path} ${read1} ${read2} > ${sample_id}_vs_HLA.sam
    """
}

// Step 2: Run HLAminer for HLA typing
process HLAMINER_PREDICT {
    tag "${sample_id}"
    label 'process_medium'
        
    publishDir(
        path: {
            def base_sample_id = sample_id.replaceAll('_(dna|rna)_(normal|tumor)$', '')
            return "${params.outdir}/${base_sample_id}/rna/09_hla_typing"
        },
        mode: 'copy',
        saveAs: { filename ->
            return "${sample_id}/${filename}"
        }
    )
    
    input:
    tuple val(sample_id), path(sam_file)
    
    output:
    tuple val(sample_id), 
          path("${sample_id}_hla_predictions.csv"),
          path("${sample_id}_HLAminer_HPRA.csv"),
          path("${sample_id}_HLAminer.log"),
          emit: hla_results
    
    script:
    def base_sample_id = sample_id.replaceAll('_(dna|rna)_(normal|tumor)$', '')
    def output_dir = "${params.outdir}/${base_sample_id}/rna/09_hla_typingii/${sample_id}"
    def hla_db_path = "${params.hlaminer.database}/HLA-I_II_CDS.fasta"
    def hla_p_path = "${params.hlaminer.database}/hla_nom_p.txt"
    def log_file = "${sample_id}_HLAminer.log"
    
    """
    # Run HLAminer for HLA typing
    perl ${params.hlaminer.path}/HLAminer.pl \\
        -a ${sam_file} \\
        -h ${hla_db_path} \\
        -p ${hla_p_path} \\
        -l ${log_file}
    
    # Rename output files for better identification
    cp HLAminer_HPRA_${log_file}.csv ${sample_id}_HLAminer_HPRA.csv
    
    # Check if .csv.pred file exists, if it exists then copy it, otherwise use .csv file as a substitute
    if [ -f HLAminer_HPRA_${log_file}.csv.pred ]; then
        cp HLAminer_HPRA_${log_file}.csv.pred ${sample_id}_hla_predictions.csv
    else
        # If .csv.pred file doesn't exist, use .csv file as a substitute
        cp HLAminer_HPRA_${log_file}.csv ${sample_id}_hla_predictions.csv
    fi
    
    # Save log file
    cp HLAminer_HPRA_${log_file}.log ${log_file}
    """
}

// Complete HLAminer HLA typing workflow
workflow HLAMINER_WORKFLOW {
    take:
    reads_ch  // Input channel containing sample ID and paired reads
    
    main:
    // 1. Use BWA to align reads to HLA reference sequences
    bwa_result = BWA_ALIGN(reads_ch)
    
    // 2. Run HLAminer for HLA typing
    hlaminer_result = HLAMINER_PREDICT(bwa_result.align_result)
    
    emit:
    hla_results = hlaminer_result.hla_results
}
