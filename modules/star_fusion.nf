// STAR-Fusion 进程用于RNA肿瘤样本融合检测
process STAR_FUSION {
    tag "${sample_id}"
    label 'process_high'
    
    container "${params.star_fusion.container}"
    
    // 添加容器绑定选项
    containerOptions "--bind /mnt/data/lumanman/db/star-fusion:/db"
    
    publishDir(
        path: {
            def base_sample_id = sample_id.replaceAll('_(dna|rna)_(normal|tumor)$', '')
            return "${params.outdir}/${base_sample_id}/rna/08_star_fusion"
        },
        mode: 'copy',
        saveAs: { filename ->
            if (filename.endsWith('.tsv')) return filename
            else null
        }
    )
    
    input:
    tuple val(sample_id), path(r1), path(r2)
    
    output:
    tuple val(sample_id),
          path("star_fusion_results/star-fusion.fusion_predictions.tsv"),
          path("star_fusion_results/star-fusion.fusion_predictions.abridged.tsv"),
          path("star_fusion_results/star-fusion.fusion_predictions.abridged.coding_effect.tsv"),
          emit: fusion_results
    
    script:
    def base_sample_id = sample_id.replaceAll('_(dna|rna)_(normal|tumor)$', '')
    def genome_lib_path = "/db/GRCh38_gencode_v22_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir"
    
    """
    # 创建输出目录
    mkdir -p star_fusion_results
    
    # 设置ulimit
    ulimit -n 10000
    
    # 运行STAR-Fusion
    STAR-Fusion \\
        --genome_lib_dir ${genome_lib_path} \\
        --left_fq ${r1} \\
        --right_fq ${r2} \\
        --FusionInspector inspect \\
        --examine_coding_effect \\
        --extract_fusion_reads \\
        --CPU ${task.cpus} \\
        --output_dir star_fusion_results
    """
}
