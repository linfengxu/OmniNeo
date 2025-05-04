process NONFRAMESHIFT_DELETION_ANALYSIS {
    tag "$sample_id"
    label 'process_medium'
    container "${params.snv_analysis.container}"  // 复用同样的Python容器
    
    publishDir(
        path: {
            def base_sample_id = sample_id.replaceAll('_(dna|rna)_(normal|tumor)$', '')
            def data_type = "dna"
            if (sample_id.contains("_rna_")) {
                data_type = "rna"
            }
            return "${params.outdir}/${base_sample_id}/${data_type}/10_peptides/nonframeshift_deletion"
        },
        mode: 'copy',
        saveAs: { filename ->
            if (filename.endsWith('.csv')) return "results/${filename}"
            else if (filename.endsWith('.fasta')) return "fasta/${filename}"
            else if (filename.endsWith('.txt')) return "txt/${filename}"
            else null
        }
    )
    
    input:
    tuple val(sample_id), path(anno_result)
    
    output:
    tuple val(sample_id), 
          path("csv_files/*.csv"), 
          optional: true,
          emit: nfsdel_csv_files
    tuple val(sample_id), 
          path("fasta_files/*.fasta"), 
          optional: true,
          emit: nfsdel_fasta_files
    tuple val(sample_id), 
          path("txt_files/*.txt"), 
          optional: true,
          emit: nfsdel_txt_files
    tuple val(sample_id), 
          path("NFSDEL_gene_summary.csv"), 
          emit: nfsdel_gene_summary
    tuple val(sample_id), 
          path("NFSDEL_Varsequence.fasta"), 
          emit: nfsdel_var_sequence
    
    script:
    """
    # Create output directories
    mkdir -p csv_files fasta_files txt_files
    
    # Run nonframeshift deletion analysis
    python ${projectDir}/bin/DNA/nonframeshift_deletion.py \\
        -r ${params.uniport_ref} \\
        -i ${anno_result} \\
        -o ./  \\
        -t ${params.kallisto.transcriptome_fasta} \\
        -p -v
        
    # Print analysis completion message for logging
    echo "Nonframeshift deletion analysis completed for sample ${sample_id}"
    echo "Generated output files:"
    ls -la csv_files/ fasta_files/ txt_files/ NFSDEL_gene_summary.csv NFSDEL_Varsequence.fasta 2>/dev/null || true
    """
}
