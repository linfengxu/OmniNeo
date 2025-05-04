process NONFRAMESHIFT_SUBSTITUTION_ANALYSIS {
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
            return "${params.outdir}/${base_sample_id}/${data_type}/10_peptides/nonframeshift_substitution"
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
          emit: nfssub_csv_files
    tuple val(sample_id), 
          path("fasta_files/*.fasta"), 
          optional: true,
          emit: nfssub_fasta_files
    tuple val(sample_id), 
          path("txt_files/*.txt"), 
          optional: true,
          emit: nfssub_txt_files
    tuple val(sample_id), 
          path("AASUB_gene_summary.csv"), 
          emit: nfssub_gene_summary
    tuple val(sample_id), 
          path("AASUB_Varsequence.fasta"), 
          emit: nfssub_var_sequence
    
    script:
    """
    # Create output directories
    mkdir -p csv_files fasta_files txt_files
    
    # Run nonframeshift substitution analysis
    python ${projectDir}/bin/DNA/nonframeshift_substitution.py \\
        -r ${params.uniport_ref} \\
        -i ${anno_result} \\
        -o ./  \\
        -t ${params.kallisto.transcriptome_fasta} \\
        --debug -p -v
        
    # Print analysis completion message for logging
    echo "Nonframeshift substitution analysis completed for sample ${sample_id}"
    echo "Generated output files:"
    ls -la csv_files/ fasta_files/ txt_files/ AASUB_gene_summary.csv AASUB_Varsequence.fasta 2>/dev/null || true
    """
}
