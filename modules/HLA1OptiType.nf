// HLA-I typing module, using OptiType for HLA typing of RNA-seq samples

// Step 1: Use Razers3 to extract reads matching HLA reference sequences
process RAZERS3_EXTRACT {
    tag "${sample_id}"
    label 'process_medium'
    
    container "${params.optitype.container}"
    
    containerOptions "--bind  ${params.optitype.hla_reference_rna}:/data"
    
    input:
    tuple val(sample_id), path(read1), path(read2)
    
    output:
    tuple val(sample_id), 
          path("${sample_id}_fished_1.bam"),
          path("${sample_id}_fished_2.bam"),
          emit: extract_result
    
    script:
    def optitype_tmp = "${params.tmp_dir}/${sample_id}_optitype"
    def hla_ref_path = "/data/hla_reference_rna.fasta"
    
    """
    # Step 1: Use Razers3 to extract reads matching HLA reference sequences - forward reads
    razers3 -i 95 -m 1 -dr 0 \\
        -o ${sample_id}_fished_1.bam \\
        ${hla_ref_path} \\
        ${read1}
        
    # Step 2: Use Razers3 to extract reads matching HLA reference sequences - reverse reads
    razers3 -i 95 -m 1 -dr 0 \\
        -o ${sample_id}_fished_2.bam \\
        ${hla_ref_path} \\
        ${read2}
    """
}

// Step 2: Convert BAM files to FASTQ format
process SAMTOOLS_BAM_TO_FASTQ {
    tag "${sample_id}"
    label 'process_low'
    
    container "${params.optitype.samtools_container}"
    
    input:
    tuple val(sample_id), path(fished_bam1), path(fished_bam2)
    
    output:
    tuple val(sample_id), 
          path("${sample_id}_fished_1.fastq"),
          path("${sample_id}_fished_2.fastq"),
          emit: fastq_result
    
    script:
    """
    # Convert BAM file to FASTQ format - forward reads
    samtools bam2fq ${fished_bam1} > ${sample_id}_fished_1.fastq
    
    # Convert BAM file to FASTQ format - reverse reads
    samtools bam2fq ${fished_bam2} > ${sample_id}_fished_2.fastq
    """
}

// Run OptiType HLA typing
process OPTITYPE_TYPING {
    tag "${sample_id}"
    label 'process_medium'
    
    container "${params.optitype.container}"
    
    containerOptions "--bind ${params.optitype.hla_reference_rna}:/data"
    
    publishDir(
        path: {
            def base_sample_id = sample_id.replaceAll('_(dna|rna)_(normal|tumor)$', '')
            return "${params.outdir}/${base_sample_id}/rna/11_hla_typing"
        },
        mode: 'copy',
        saveAs: { filename ->
            return "${sample_id}/${filename}"
        }
    )
    
    input:
    tuple val(sample_id), path(fished_fastq1), path(fished_fastq2)
    
    output:
    tuple val(sample_id), 
          path("${sample_id}_result.tsv"),
          path("${sample_id}_coverage_plot.pdf"),
          emit: hla_results
    
    script:
    def base_sample_id = sample_id.replaceAll('_(dna|rna)_(normal|tumor)$', '')
    def output_dir = "${params.outdir}/${base_sample_id}/rna/11_hla_typing/${sample_id}"
    def result_tsv = "${output_dir}/${sample_id}_result.tsv"
    def coverage_pdf = "${output_dir}/${sample_id}_coverage_plot.pdf"
    
    """
    if [[ -f "${result_tsv}" && -f "${coverage_pdf}" ]]; then
        echo "HLA typing results already exist for ${sample_id}, creating symlinks..."
        ln -s "${result_tsv}" "${sample_id}_result.tsv"
        ln -s "${coverage_pdf}" "${sample_id}_coverage_plot.pdf"
    else
        mkdir -p optitype_out
        
        OptiTypePipeline.py \\
            -i ${fished_fastq1} ${fished_fastq2} \\
            -r \\
            -o ./ \\
            -p ${sample_id} \\
            -v
    fi
    """
}

// HLA-I workflow
workflow HLA_I_WORKFLOW {
    take:
    reads_ch  
    
    main:
  
    razers3_result = RAZERS3_EXTRACT(reads_ch)
    
    samtools_result = SAMTOOLS_BAM_TO_FASTQ(razers3_result.extract_result)
    
    optitype_result = OPTITYPE_TYPING(samtools_result.fastq_result)
    
    emit:
    hla_results = optitype_result.hla_results
}
