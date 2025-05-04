process NONFRAMESHIFT_INSERTION_ANALYSIS {
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
            return "${params.outdir}/${base_sample_id}/${data_type}/10_peptides/nonframeshift_insertion"
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
          emit: nfsins_csv_files
    tuple val(sample_id), 
          path("fasta_files/*.fasta"), 
          optional: true,
          emit: nfsins_fasta_files
    tuple val(sample_id), 
          path("txt_files/*.txt"), 
          optional: true,
          emit: nfsins_txt_files
    tuple val(sample_id), 
          path("NFSINS_gene_summary.csv"), 
          emit: nfsins_gene_summary
    tuple val(sample_id), 
          path("NFSINS_Varsequence.fasta"), 
          emit: nfsins_var_sequence
    
    script:
    """
    # Create output directories
    mkdir -p csv_files fasta_files txt_files
    
    # Run nonframeshift insertion analysis
    python ${projectDir}/bin/DNA/nonframeshift_insertion.py \\
        -r ${params.uniport_ref} \\
        -i ${anno_result} \\
        -o ./  \\
        -t ${params.kallisto.transcriptome_fasta} \\
        --debug -p -v
        
    # Print analysis completion message for logging
    echo "Nonframeshift insertion analysis completed for sample ${sample_id}"
    echo "Generated output files:"
    ls -la csv_files/ fasta_files/ txt_files/ NFSINS_gene_summary.csv NFSINS_Varsequence.fasta 2>/dev/null || true
    """
}
