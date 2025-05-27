#!/usr/bin/env nextflow

process FORMAT_HLA_I_TYPES {
    tag "$sample_id"
    label 'process_low'
    
    publishDir(
        path: { 
            def base_sample_id = sample_id.replaceAll('_(dna|rna)_(normal|tumor)$', '')
            return "${params.outdir}/${base_sample_id}/hla_typing/class_i" 
        },
        mode: 'copy'
    )
    
    input:
    tuple val(sample_id), path(hla_i_result), path(coverage_plot)
    
    output:
    tuple val(sample_id), path("${sample_id}_hla_i_formatted.txt"), emit: formatted_hla_i
    
    script:
    """
    # Read OPTITYPE HLA-I results and format
    awk -F'\t' '
    BEGIN {
        print "# HLA-I alleles from OPTITYPE"
        print "# Format: HLA-A02:01,HLA-A24:02,..."
    }
    NR==2 {
        # A alleles: \$2, \$3
        # B alleles: \$4, \$5
        # C alleles: \$6, \$7
        printf("HLA-A%s,HLA-A%s,HLA-B%s,HLA-B%s,HLA-C%s,HLA-C%s\\n", 
            \$2, \$3, \$4, \$5, \$6, \$7)
    }
    ' ${hla_i_result} > ${sample_id}_hla_i_formatted.txt

    # Fix format: remove extra A
    sed -i \
        -e 's/HLA-AA/HLA-A/g' \
        -e 's/HLA-BB/HLA-B/g' \
        -e 's/HLA-CC/HLA-C/g' \
        "${sample_id}_hla_i_formatted.txt"
    """
}

process FORMAT_HLA_II_TYPES {
    tag "$sample_id"
    label 'process_low'
    
    publishDir(
        path: { 
            def base_sample_id = sample_id.replaceAll('_(dna|rna)_(normal|tumor)$', '')
            return "${params.outdir}/${base_sample_id}/hla_typing/class_ii" 
        },
        mode: 'copy'
    )
    
    input:
    tuple val(sample_id), path(hla_predictions), path(hla_ii_result), path(hla_log)
    
    output:
    tuple val(sample_id), path("${sample_id}_hla_ii_formatted.txt"), emit: formatted_hla_ii
    
    script:
    """
    awk -F',' '
    BEGIN {
        print "# HLA-II alleles for netMHCIIpan"
        print "# Format: DRB1_0101,DRB3_0210,DRB4_0101"
    }
    /DRB1\\*[0-9]+:[0-9]+/ {
        sub(/P\$/,"",\$1)
        if (match(\$1, /DRB1\\*([0-9]+):([0-9]+)/, arr)) {
            drb1_alleles = drb1_alleles (drb1_alleles=="" ? "" : ",") sprintf("DRB1_%02d%02d", arr[1], arr[2])
        }
    }
    /DRB3\\*[0-9]+:[0-9]+/ {
        sub(/P\$/,"",\$1)
        if (match(\$1, /DRB3\\*([0-9]+):([0-9]+)/, arr)) {
            drb3_alleles = drb3_alleles (drb3_alleles=="" ? "" : ",") sprintf("DRB3_%02d%02d", arr[1], arr[2])
        }
    }
    /DRB4\\*[0-9]+:[0-9]+/ {
        sub(/P\$/,"",\$1)
        if (match(\$1, /DRB4\\*([0-9]+):([0-9]+)/, arr)) {
            drb4_alleles = drb4_alleles (drb4_alleles=="" ? "" : ",") sprintf("DRB4_%02d%02d", arr[1], arr[2])
        }
    }
    END {
        print drb1_alleles (drb1_alleles!="" && drb3_alleles!="" ? "," : "") drb3_alleles (((drb1_alleles!="" || drb3_alleles!="") && drb4_alleles!="") ? "," : "") drb4_alleles
    }
    ' ${hla_ii_result} > ${sample_id}_hla_ii_formatted.txt
    """
}

workflow HLA_TYPING {
    take:
    hla_i_typing_results    // [sample_id, result.tsv, coverage_plot.pdf]
    hla_ii_typing_results   // [sample_id, predictions.csv, HLAminer_HPRA.csv, HLAminer.log]
    
    main:
    formatted_hla_i = FORMAT_HLA_I_TYPES(hla_i_typing_results)
    formatted_hla_ii = FORMAT_HLA_II_TYPES(hla_ii_typing_results)
    
    emit:
    formatted_hla_i = formatted_hla_i.formatted_hla_i
    formatted_hla_ii = formatted_hla_ii.formatted_hla_ii
}

