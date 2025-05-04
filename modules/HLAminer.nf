// HLA-II分型模块，使用HLAminer对RNA-seq样本进行肿瘤样本HLA分型

// 第一步：使用BWA mem将读数比对到HLA参考序列
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
    def hla_ref_path = "/data/HLA_ABC_CDS.fasta"
    def threads = params.bwa.threads
    
    """
    # 使用BWA mem将读数比对到HLA参考序列
    bwa mem -a -t ${threads} ${hla_ref_path} ${read1} ${read2} > ${sample_id}_vs_HLA.sam
    """
}

// 第二步：运行HLAminer进行HLA分型
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
    def hla_db_path = "${params.hlaminer.database}/HLA_ABC_CDS.fasta"
    def hla_p_path = "${params.hlaminer.database}/hla_nom_p.txt"
    
    """
    # 运行HLAminer进行HLA分型
    perl ${params.hlaminer.path}/HLAminer.pl \\
        -a ${sam_file} \\
        -h ${hla_db_path} \\
        -p ${hla_p_path} \\
        -l ${sample_id}_HLAminer.log
    
    # 重命名输出文件以便更好地识别
    cp HLAminer_HPRA.csv ${sample_id}_HLAminer_HPRA.csv
    cp HLAminer_predictions.csv ${sample_id}_hla_predictions.csv
    """
}

// 完整的HLAminer HLA分型工作流
workflow HLAMINER_WORKFLOW {
    take:
    reads_ch  // 输入通道，含样本ID和配对读数
    
    main:
    // 1. 使用BWA比对读数到HLA参考序列
    bwa_result = BWA_ALIGN(reads_ch)
    
    // 2. 运行HLAminer进行HLA分型
    hlaminer_result = HLAMINER_PREDICT(bwa_result.align_result)
    
    emit:
    hla_results = hlaminer_result.hla_results
}
