// Kallisto module，用于RNA-seq表达量定量

// Kallisto index create
process KALLISTO_INDEX {
    tag "kallisto_index"
    label 'process_medium'
    
    container "${params.kallisto.container}"
    
    publishDir(
        path: "${params.outdir}/references/kallisto",
        mode: 'copy'
    )
    
    input:
    path transcriptome_fasta
    
    output:
    path "${params.kallisto.index_name}", emit: kallisto_index
    
    script:
    """
    # 检查索引是否已存在
    if [[ -f "${params.kallisto.index_file}" ]]; then
        echo "Kallisto index already exists at ${params.kallisto.index_file}, creating symlink..."
        ln -s "${params.kallisto.index_file}" "${params.kallisto.index_name}"
    else
        echo "Building Kallisto index..."
        kallisto index -i "${params.kallisto.index_name}" "${transcriptome_fasta}"
    fi
    """
}

// Kallisto定量
process KALLISTO_QUANT {
    tag "${sample_id}"
    label 'process_medium'
    
    container "${params.kallisto.container}"
    
    publishDir(
        path: {
            def base_sample_id = sample_id.replaceAll('_(dna|rna)_(normal|tumor)$', '')
            return "${params.outdir}/${base_sample_id}/rna/09_kallisto"
        },
        mode: 'copy',
        saveAs: { filename ->
            if (filename == "abundance.tsv") return "${sample_id}/abundance.tsv"
            else if (filename == "abundance.h5") return "${sample_id}/abundance.h5"
            else if (filename == "run_info.json") return "${sample_id}/run_info.json"
            else return "${sample_id}/${filename}"
        }
    )
    
    input:
    tuple val(sample_id), path(read1), path(read2)
    path kallisto_index
    
    output:
    tuple val(sample_id), 
          path("abundance.tsv"),
          path("abundance.h5"),
          path("run_info.json"),
          emit: quant_results
    
    script:
    def base_sample_id = sample_id.replaceAll('_(dna|rna)_(normal|tumor)$', '')
    def output_dir = "${params.outdir}/${base_sample_id}/rna/07_kallisto/${sample_id}"
    
    """
    # 检查输出是否已经存在
    if [[ -f "${output_dir}/abundance.tsv" && -f "${output_dir}/abundance.h5" && -f "${output_dir}/run_info.json" ]]; then
        echo "Kallisto quant output already exists for ${sample_id}, creating symlinks..."
        ln -s "${output_dir}/abundance.tsv" "abundance.tsv"
        ln -s "${output_dir}/abundance.h5" "abundance.h5"
        ln -s "${output_dir}/run_info.json" "run_info.json"
    else
        echo "Running Kallisto quant on ${sample_id}..."
        kallisto quant \\
            -i "${kallisto_index}" \\
            -o ./ \\
            -b ${params.kallisto.bootstrap_samples} \\
            -t ${task.cpus} \\
            "${read1}" "${read2}"
    fi
    """
}
