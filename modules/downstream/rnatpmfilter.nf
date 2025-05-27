process RNA_TPM_FILTER {
    tag "$sample_id"
    label 'process_low'
    container "${params.snv_analysis.container}"
    
    publishDir(
        path: {
            def base_sample_id = sample_id.replaceAll('_(dna|rna)_(normal|tumor)$', '')
            def data_type = "rna"
            return "${params.outdir}/${base_sample_id}/${data_type}/09_kallisto/filtered_tpm"
        },
        mode: 'copy',
        saveAs: { filename ->
            if (filename.endsWith('.txt')) return "gene_expression/${filename}"
            else if (filename.endsWith('.tsv')) return "variants/${filename}"
            else null
        }
    )
    
    input:
    tuple val(sample_id), path(kallisto_abundance), path(abundance_h5), path(run_info)
    path(gene_reference)
    tuple val(anno_sample_id), path(annovar_results)
    
    output:
    tuple val(sample_id), 
          path("anovar_filtered_variants.tsv"), 
          emit: filtered_variants
    tuple val(sample_id), 
          path("kallisto_gene_expression.txt"), 
          emit: gene_expression
    
    script:
    def tpm_threshold = params.rna_tpm_filter ? params.rna_tpm_filter.threshold : 0
    """
    # run TPM filter analyses
    python ${projectDir}/bin/RNA/TpmFilter.py \\
        -k ${kallisto_abundance} \\
        -r ${gene_reference} \\
        -a ${annovar_results} \\
        -o anovar_filtered_variants.tsv \\
        -t ${tpm_threshold}
        
    # log messgages
    echo "RNA TPM filtering completed for sample ${sample_id}"
    echo "TPM threshold used: ${tpm_threshold}"
    echo "Generated output files:"
    ls -la filtered_variants.tsv kallisto_gene_expression.txt 2>/dev/null || true
    """
}
