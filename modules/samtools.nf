// 添加一个用于 fixmate 的进程
process SAMTOOLS_FIXMATE {
    tag "${sample_id}"
    label 'process_medium'
    
    container { "${params.samtools.container}" }
    
    publishDir path: {
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
        
        return "${params.outdir}/${base_sample_id}/${data_type}/03_samtools/fixmate/${sample_type}"
    },
    mode: 'copy',
    saveAs: { filename -> 
        if (filename.endsWith('.bam')) "fixmate/${filename}" 
        else null
    }

    input:
    tuple val(sample_id), path(align_sam)

    output:
    tuple val(sample_id), path("${sample_id}_align_fixmate.bam")
    
    script:
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
    
    def output_dir = "${params.outdir}/${base_sample_id}/${data_type}/03_samtools/fixmate/${sample_type}"
    def fixmate_bam = "${output_dir}/fixmate/${sample_id}_align_fixmate.bam"
    """
    if [[ -f "${fixmate_bam}" ]]; then
        echo "Fixmate BAM file already exists for ${sample_id}, creating symlink..."
        ln -s "${fixmate_bam}" "${sample_id}_align_fixmate.bam"
    else
        echo "Processing ${sample_id} with samtools fixmate..."
        samtools fixmate -O bam ${align_sam} ${sample_id}_align_fixmate.bam
    fi
    """
}

// 添加一个用于 sort 的进程
process SAMTOOLS_SORT {
    tag "${sample_id}"
    label 'process_medium'
    
    container { "${params.samtools.container}" }
    
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
        
        return "${params.outdir}/${base_sample_id}/${data_type}/03_samtools/sorted/${sample_type}"
    },
    mode: 'copy',
    saveAs: { filename -> 
        if (filename.endsWith('.bam')) "sorted/${filename}" 
        else null
    }

    input:
    tuple val(sample_id), path(fixmate_bam)

    output:
    tuple val(sample_id), path("${sample_id}_align_sorted.bam")
    
    script:
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
    
    def output_dir = "${params.outdir}/${base_sample_id}/${data_type}/03_samtools/sorted/${sample_type}"
    def sorted_bam = "${output_dir}/sorted/${sample_id}_align_sorted.bam"
    def tmp_dir = "tmp_${sample_id}"
    """
    if [[ -f "${sorted_bam}" ]]; then
        echo "Sorted BAM file already exists for ${sample_id}, creating symlink..."
        ln -s "${sorted_bam}" "${sample_id}_align_sorted.bam"
    else
        echo "Processing ${sample_id} with samtools sort..."
        mkdir -p ${tmp_dir}
        samtools sort -O bam \\
            -o ${sample_id}_align_sorted.bam \\
            -T ${tmp_dir}/${sample_id} \\
            ${fixmate_bam}
        rm -rf ${tmp_dir}
    fi
    """
}

// Samtools Index process
process SAMTOOLS_INDEX {
    tag "${sample_id}"
    label 'process_low'
    
    container { "${params.samtools.container}" }
    
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
        
        return "${params.outdir}/${base_sample_id}/${data_type}/05_picard/read_groups/bam"
    },
    mode: 'copy'

    input:
    tuple val(sample_id), path(final_bam)

    output:
    tuple val(sample_id), 
          path("${final_bam}"),
          path("${final_bam}.bai"),
          emit: indexed_bam

    script:
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
    
    def output_dir = "${params.outdir}/${base_sample_id}/${data_type}/05_picard/read_groups/bam"
    def bai_file = "${output_dir}/${final_bam}.bai"

    """
    if [[ -f "${bai_file}" ]]; then
        echo "Index file already exists for ${final_bam}, creating symlink..."
        ln -s "${bai_file}" "${final_bam}.bai"
    else
        echo "Creating index for ${final_bam}..."
        samtools index ${final_bam}
    fi
    """
}
