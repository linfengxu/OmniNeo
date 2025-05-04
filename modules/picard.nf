// Picard AddOrReplaceReadGroups process
process PICARD_ADD_OR_REPLACE_READ_GROUPS {
    tag "${sample_id}"
    label 'process_medium'
    
    container { "${params.picard.container}" }
    
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
        
        return "${params.outdir}/${base_sample_id}/${data_type}/05_picard/read_groups"
    },
    mode: 'copy',
    saveAs: { filename ->
        if (filename.endsWith('.bam')) "bam/${filename}"
        else null
    }

    input:
    tuple val(sample_id), path(recal_bam)

    output:
    tuple val(sample_id), 
          path("${sample_id}_final.bam"),
          emit: final_bam

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
    
    def output_dir = "${params.outdir}/${base_sample_id}/${data_type}/05_picard/read_groups"
    def output_bam = "${output_dir}/bam/${sample_id}_final.bam"

    """
    if [[ -f "${output_bam}" ]]; then
        echo "Picard AddOrReplaceReadGroups output file already exists for ${sample_id}, creating symlink..."
        ln -s "${output_bam}" "${sample_id}_final.bam"
    else
        echo "Processing ${sample_id} with Picard AddOrReplaceReadGroups..."
        
        picard AddOrReplaceReadGroups \\
            I=${recal_bam} \\
            O=${sample_id}_final.bam \\
            RGID=${sample_id} \\
            RGLB=${params.picard.read_groups.library} \\
            RGPL=${params.picard.read_groups.platform} \\
            RGPU=${params.picard.read_groups.platform_unit} \\
            SORT_ORDER=coordinate \\
            RGSM=${sample_id}
    fi
    """
}
