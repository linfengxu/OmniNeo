process RNA_STAR_FUSION_PEPTIDES {
    tag "$sample_id"
    label 'process_medium'
    container "${params.snv_analysis.container}"  // 复用同样的Python容器
    
    publishDir(
        path: {
            def base_sample_id = sample_id.replaceAll('_(dna|rna)_(normal|tumor)$', '')
            def data_type = "rna"
            return "${params.outdir}/${base_sample_id}/${data_type}/10_fusion/peptides"
        },
        mode: 'copy',
        saveAs: { filename ->
            if (filename.endsWith('.fasta')) return "fasta/${filename}"
            else if (filename.endsWith('.csv')) return "csv/${filename}"
            else null
        }
    )
    
    input:
    tuple val(sample_id), path(fusion_results)
    
    output:
    tuple val(sample_id), 
          path("fasta/*_short_peptides.fasta"), 
          optional: true,
          emit: fusion_short_peptides
    tuple val(sample_id), 
          path("fasta/*_long_peptides.fasta"), 
          optional: true,
          emit: fusion_long_peptides
    tuple val(sample_id), 
          path("csv/*_short_peptides.csv"), 
          optional: true,
          emit: fusion_short_csv
    tuple val(sample_id), 
          path("csv/*_long_peptides.csv"), 
          optional: true,
          emit: fusion_long_csv
    
    script:
    def short_min = params.fusion_peptides?.short_min ?: 8
    def short_max = params.fusion_peptides?.short_max ?: 11
    def long_min = params.fusion_peptides?.long_min ?: 15
    def long_max = params.fusion_peptides?.long_max ?: 30
    """
    mkdir -p fasta csv
    
    # 检查融合结果文件是否存在且不为空
    if [ -s "${fusion_results}" ]; then
        # 运行融合肽段分析
        python ${projectDir}/bin/RNA/star_fusion.py \\
            -i ${fusion_results} \\
            -o ./ \\
            --short_min ${short_min} \\
            --short_max ${short_max} \\
            --long_min ${long_min} \\
            --long_max ${long_max}
            
        # 打印分析完成信息
        echo "RNA fusion peptide generation completed for sample ${sample_id}"
        echo "Generated fusion peptides:"
        ls -la fasta/ csv/ 2>/dev/null || true
    else
        echo "WARNING: No fusion events found or empty results file for sample ${sample_id}"
        # 创建空文件以确保输出通道正常工作
        touch fasta/${sample_id}_short_peptides.fasta
        touch fasta/${sample_id}_long_peptides.fasta
        touch csv/${sample_id}_short_peptides.csv
        touch csv/${sample_id}_long_peptides.csv
    fi
    """
}
