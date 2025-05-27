process ANNOVAR {
    tag "$sample_id"
    label 'process_high'
    
    publishDir(
        path: {
            def base_sample_id = sample_id.replaceAll('_(dna|rna)_(normal|tumor)$', '')
            
            //  data type（DNA or RNA） 
            def data_type = "dna"
            if (sample_id.contains("_rna_")) {
                data_type = "rna"
            }
            
            return "${params.outdir}/${base_sample_id}/${data_type}/07_anovar"
        },
        mode: 'copy',
        saveAs: { filename ->
            if (filename.endsWith('.txt')) return "results/${filename}"
            else if (filename.endsWith('.vcf')) return "vcf/${filename}"
            else null
        }
    )
    
    input:
    tuple val(sample_id), path(input_vcf)
    
    output:
    tuple val(sample_id),
          path("mut_result/${sample_id}_anno_out.hg38_multianno.txt"),
          emit: anovar_txt
    tuple val(sample_id),
          path("mut_result/${sample_id}_anno_out.hg38_multianno.vcf"),
          optional: true,
          emit: anovar_vcf
    
    script:
    """
    # Create output directory
    mkdir -p mut_result
    
    # Print sample ID for debugging
    echo "Processing sample ID: ${sample_id}, Input VCF: ${input_vcf}"
    
    # Copy the input file to the output directory with the name expected by the ANNOVAR command
    # Using cp instead of symlink to ensure file exists in the expected location
    cp ${input_vcf} mut_result/${sample_id}_filter_somatic.vcf
    
    # Run ANNOVAR using the installed location
    perl ${params.annovar_path}/table_annovar.pl \\
        mut_result/${sample_id}_filter_somatic.vcf \\
        ${params.annovar_db_path} \\
        -buildver ${params.genome_build} \\
        -out mut_result/${sample_id}_anno_out \\
        -remove \\
        -protocol refGene \\
        -operation g \\
        -nastring . \\
        -vcfinput
    """
}