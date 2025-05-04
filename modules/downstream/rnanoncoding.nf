process RNA_NONCODING_ANALYSIS {
    tag "$sample_id"
    label 'process_medium'
    container "${params.snv_analysis.container}"  // 复用同样的Python容器
    
    publishDir(
        path: {
            def base_sample_id = sample_id.replaceAll('_(dna|rna)_(normal|tumor)$', '')
            def data_type = "rna"
            return "${params.outdir}/${base_sample_id}/${data_type}/10_peptides/noncoding"
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
          path("noncoding_mutation_peptides.csv"), 
          emit: noncoding_peptides_csv
    tuple val(sample_id), 
          path("fasta/noncoding_short_peptides.fasta"), 
          emit: noncoding_short_peptides
    tuple val(sample_id), 
          path("fasta/noncoding_long_peptides.fasta"), 
          optional: true,
          emit: noncoding_long_peptides
    tuple val(sample_id), 
          path("fasta/noncoding_long_peptides_for_ms.fasta"), 
          optional: true,
          emit: noncoding_ms_peptides
    tuple val(sample_id), 
          path("long_peptide_stats.txt"), 
          emit: noncoding_stats
    
    script:
    """
    # 运行RNA非编码区突变分析
    python ${projectDir}/bin/RNA/noncoding.py \\
        -i ${anno_result} \\
        -r ${params.ref_genome_fasta} \\
        -o ./ \\
        --short_min 8 \\
        --short_max 11 \\
        --long_min 15 \\
        --long_max 30
        
    # 打印分析完成信息
    echo "RNA noncoding mutation analysis completed for sample ${sample_id}"
    echo "Generated output files:"
    ls -la noncoding_mutation_peptides.csv fasta/ long_peptide_stats.txt 2>/dev/null || true
    """
}
