// modules/fastp.nf
process FASTP {
    tag "${sample_id}"
    label 'process_medium'
    
    // 使用Singularity容器
    container { "${params.fastp.container}" }
    
    // 发布输出文件到结果目录，根据样本ID区分DNA/RNA和normal/tumor
    publishDir {
        def base_sample_id = sample_id.replaceAll('_(dna|rna)_(normal|tumor)$', '')
        
        // 判断数据类型（DNA或RNA）
        def data_type = "dna"
        if (sample_id.contains("_rna_")) {
            data_type = "rna"
        }
        
        // 判断样本类型（normal或tumor）
        def sample_type = "normal"
        if (sample_id.contains("_tumor")) {
            sample_type = "tumor"
        }
        
        return "${params.outdir}/${base_sample_id}/${data_type}/01_fastp/${sample_type}"
    },
        mode: 'copy',
        saveAs: { filename ->
            if (filename.endsWith('.fastq.gz')) "cleaned/$filename"
            else if (filename.endsWith('.json')) "reports/$filename"
            else if (filename.endsWith('.html')) "reports/$filename"
            else null
        }
    
    input:
    tuple val(sample_id), path(read1), path(read2)
    
    output:
    tuple val(sample_id),
          path("${sample_id}_clean_1.fastq.gz"),
          path("${sample_id}_clean_2.fastq.gz"),
          emit: cleaned_reads
    tuple val(sample_id),
          path("${sample_id}_fastp.json"),
          path("${sample_id}_fastp.html"),
          emit: reports
    
    script:
    def adapter_auto = params.fastp.adapters_auto ? "--detect_adapter_for_pe" : ""
    def cut_front = params.fastp.cut_front ? "--cut_front" : ""
    def cut_tail = params.fastp.cut_tail ? "--cut_tail" : ""
    
    """
    fastp \\
        -i ${read1} \\
        -I ${read2} \\
        -o ${sample_id}_clean_1.fastq.gz \\
        -O ${sample_id}_clean_2.fastq.gz \\
        -j ${sample_id}_fastp.json \\
        -h ${sample_id}_fastp.html \\
        -w ${params.fastp.threads} \\
        ${adapter_auto} \\
        ${cut_front} \\
        ${cut_tail} \\
        --cut_window_size ${params.fastp.cut_window_size} \\
        --cut_mean_quality ${params.fastp.cut_mean_quality} \\
        --qualified_quality_phred ${params.fastp.qualified_quality} \\
        --unqualified_percent_limit ${params.fastp.unqualified_percent_limit} \\
        --n_base_limit ${params.fastp.n_base_limit} \\
        --length_required ${params.fastp.min_length}
    """
}