// add create index
process BWA_INDEX {
    tag "Creating BWA index"
    label 'process_high'
    
    container { "${params.bwa.container}" }
    
    input:
    path reference
    
    output:
    path "${reference}*", emit: index
    
    script:
    """
    # Get the actual path of the reference genome
    ref_path=\$(readlink -f ${reference})
    ref_dir=\$(dirname \$ref_path)
    ref_base=\$(basename \$ref_path)
    
    cd \$ref_dir
    
    # Check BWA index files
    missing_bwa_index=0
    for suffix in bwt pac ann amb sa; do
        if [[ ! -f "\${ref_base}.\${suffix}" ]]; then
            missing_bwa_index=1
            break
        fi
    done
    
    # Create BWA index if needed
    if [[ \$missing_bwa_index -eq 1 ]]; then
        echo "Creating BWA index in \$ref_dir..."
        bwa index \$ref_base
    else
        echo "BWA index files already exist in \$ref_dir, skipping..."
    fi
    
    # Check GATK required index files
    if [[ ! -f "\${ref_base}.fai" ]]; then
        echo "Creating FASTA index file..."
        samtools faidx \$ref_base
    fi
    
    if [[ ! -f "\${ref_base}.dict" ]]; then
        echo "Creating sequence dictionary..."
    fi
    
    # Create symbolic links in working directory
    cd -
    for f in \${ref_dir}/\${ref_base}*; do
        ln -sf \$f ./\$(basename \$f)
    done
    """
}

// Modify BWA process to use pre-created index
process BWA {
    tag "${sample_id}"
    label 'process_medium'
    
    container { "${params.bwa.container}" }
    
    publishDir path: {
        def base_sample_id = sample_id.replaceAll('_(dna|rna)_(normal|tumor)$', '')
        
        // Determine data type (DNA or RNA)
        def data_type = "dna"
        if (sample_id.contains("_rna_")) {
            data_type = "rna"
        }
        
        // Determine sample type (normal or tumor)
        def sample_type = "normal"
        if (sample_id.contains("_tumor")) {
            sample_type = "tumor"
        }
        
        return "${params.outdir}/${base_sample_id}/${data_type}/02_bwa/${sample_type}"
    }, mode: 'copy', saveAs: { filename ->
        if (filename.endsWith('.sam')) "alignments/${filename}"
        else null
    }
    
    input:
    tuple val(sample_id), path(clean1), path(clean2)
    path reference  // Explicitly specify reference parameter
    path '*'       // Receive other index files
    
    output:
    tuple val(sample_id), path("${sample_id}_cut.sam")
    
    script:
    """
    # Output debug information
    echo "===== BWA Debug Information ====="
    echo "Current directory: \$(pwd)"
    echo "Sample ID: ${sample_id}"
    echo "Clean reads 1: ${clean1}"
    echo "Clean reads 2: ${clean2}"
    echo "Reference file: ${reference}"
    echo "Directory contents:"
    ls -la
    
    # Run BWA
    bwa mem -t ${params.bwa.threads} -R "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tLB:library1\\tPL:illumina" \\
        ${reference} \\
        ${clean1} ${clean2} \\
        > ${sample_id}_cut.sam
    """
}