#!/usr/bin/env nextflow

workflow NETMHCPAN1_WORKFLOW {
    take:
    peptides_8  // Channel: tuple val(sample_id), path(peptide_file)
    peptides_9  // Channel: tuple val(sample_id), path(peptide_file)
    peptides_10 // Channel: tuple val(sample_id), path(peptide_file)
    peptides_11 // Channel: tuple val(sample_id), path(peptide_file)
    hla_data    // Channel: tuple val(sample_id), path(hla_file)

    main:
    // Normalize sample_id
    def normalized_8mer = peptides_8.map { sample_id, file ->
        def norm_id = sample_id.replaceAll('_(dna|rna)_(normal|tumor)$', '')
        tuple(norm_id, file)
    }
    def normalized_9mer = peptides_9.map { sample_id, file ->
        def norm_id = sample_id.replaceAll('_(dna|rna)_(normal|tumor)$', '')
        tuple(norm_id, file)
    }
    def normalized_10mer = peptides_10.map { sample_id, file ->
        def norm_id = sample_id.replaceAll('_(dna|rna)_(normal|tumor)$', '')
        tuple(norm_id, file)
    }
    def normalized_11mer = peptides_11.map { sample_id, file ->
        def norm_id = sample_id.replaceAll('_(dna|rna)_(normal|tumor)$', '')
        tuple(norm_id, file)
    }
    def normalized_hla = hla_data.map { sample_id, file ->
        def norm_id = sample_id.replaceAll('_(dna|rna)_(normal|tumor)$', '')
        tuple(norm_id, file)
    }

    // Run predictions in parallel for different lengths
    predict_8mer = NETMHCPAN1_PREDICT_8MER(
        normalized_8mer.combine(normalized_hla, by: 0)
    )
    predict_9mer = NETMHCPAN1_PREDICT_9MER(
        normalized_9mer.combine(normalized_hla, by: 0)
    )
    predict_10mer = NETMHCPAN1_PREDICT_10MER(
        normalized_10mer.combine(normalized_hla, by: 0)
    )
    predict_11mer = NETMHCPAN1_PREDICT_11MER(
        normalized_11mer.combine(normalized_hla, by: 0)
    )

    emit:
    predictions_8mer = predict_8mer.predictions
    predictions_9mer = predict_9mer.predictions
    predictions_10mer = predict_10mer.predictions
    predictions_11mer = predict_11mer.predictions
}

process NETMHCPAN1_PREDICT_8MER {
    tag "${sample_id} - 8mer"
    label 'process_medium'
    
    publishDir(
        path: { 
            def base_sample_id = sample_id.split('_')[0]
            return "${params.outdir}/${base_sample_id}/epitope_prediction/netmhcpan1/8mer" 
        },
        mode: 'copy',
        overwrite: true
    )
    
    input:
    tuple val(sample_id), path(peptide_input_f), path(hla_input_f)
    
    output:
    tuple val(sample_id), path("${sample_id}_netmhcpan1_8mer.xls"), emit: predictions
    
    script:
    """
    #!/bin/bash
    cp "${peptide_input_f}" "peptides.fasta"
    
    # 读取HLA等位基因并移除星号
    hla_alleles=\$(cat "${hla_input_f}" | grep -v "^#" | head -n 1 | sed 's/\\*//g')
    
    # 运行netMHCpan预测
    if [ -s "peptides.fasta" ]; then
        netMHCpan \\
            -f "peptides.fasta" \\
            -a "\$hla_alleles" \\
            -l 8 \\
            -BA \\
            -xls \\
            -xlsfile "${sample_id}_netmhcpan1_8mer.xls"
        
        # 添加预测信息
        echo "Peptide length: 8" >> "${sample_id}_netmhcpan1_8mer.xls"
        echo "Prediction completed at: \$(date)" >> "${sample_id}_netmhcpan1_8mer.xls"
        echo "netMHCpan version: \$(netMHCpan -v 2>&1)" >> "${sample_id}_netmhcpan1_8mer.xls"
    else
        echo "# No 8mer peptides found in ${peptide_input_f} for sample ${sample_id}" > "${sample_id}_netmhcpan1_8mer.xls"
        echo "HLA alleles used: \$hla_alleles" >> "${sample_id}_netmhcpan1_8mer.xls"
    fi
    """
}

process NETMHCPAN1_PREDICT_9MER {
    tag "${sample_id} - 9mer"
    label 'process_medium'
    
    publishDir(
        path: { 
            def base_sample_id = sample_id.split('_')[0]
            return "${params.outdir}/${base_sample_id}/epitope_prediction/netmhcpan1/9mer" 
        },
        mode: 'copy',
        overwrite: true
    )
    
    input:
    tuple val(sample_id), path(peptide_input_f), path(hla_input_f)
    
    output:
    tuple val(sample_id), path("${sample_id}_netmhcpan1_9mer.xls"), emit: predictions
    
    script:
    """
    #!/bin/bash
    cp "${peptide_input_f}" "peptides.fasta"
    
    # 读取HLA等位基因并移除星号
    hla_alleles=\$(cat "${hla_input_f}" | grep -v "^#" | head -n 1 | sed 's/\\*//g')
    
    # 运行netMHCpan预测
    if [ -s "peptides.fasta" ]; then
        netMHCpan \\
            -f "peptides.fasta" \\
            -a "\$hla_alleles" \\
            -l 9 \\
            -BA \\
            -xls \\
            -xlsfile "${sample_id}_netmhcpan1_9mer.xls"
        
        # 添加预测信息
        echo "Peptide length: 9" >> "${sample_id}_netmhcpan1_9mer.xls"
        echo "Prediction completed at: \$(date)" >> "${sample_id}_netmhcpan1_9mer.xls"
        echo "netMHCpan version: \$(netMHCpan -v 2>&1)" >> "${sample_id}_netmhcpan1_9mer.xls"
    else
        echo "# No 9mer peptides found in ${peptide_input_f} for sample ${sample_id}" > "${sample_id}_netmhcpan1_9mer.xls"
        echo "HLA alleles used: \$hla_alleles" >> "${sample_id}_netmhcpan1_9mer.xls"
    fi
    """
}

process NETMHCPAN1_PREDICT_10MER {
    tag "${sample_id} - 10mer"
    label 'process_medium'
    
    publishDir(
        path: { 
            def base_sample_id = sample_id.split('_')[0]
            return "${params.outdir}/${base_sample_id}/epitope_prediction/netmhcpan1/10mer" 
        },
        mode: 'copy',
        overwrite: true
    )
    
    input:
    tuple val(sample_id), path(peptide_input_f), path(hla_input_f)
    
    output:
    tuple val(sample_id), path("${sample_id}_netmhcpan1_10mer.xls"), emit: predictions
    
    script:
    """
    #!/bin/bash
    cp "${peptide_input_f}" "peptides.fasta"
    
    # 读取HLA等位基因并移除星号
    hla_alleles=\$(cat "${hla_input_f}" | grep -v "^#" | head -n 1 | sed 's/\\*//g')
    
    # 运行netMHCpan预测
    if [ -s "peptides.fasta" ]; then
        netMHCpan \\
            -f "peptides.fasta" \\
            -a "\$hla_alleles" \\
            -l 10 \\
            -BA \\
            -xls \\
            -xlsfile "${sample_id}_netmhcpan1_10mer.xls"
        
        # 添加预测信息
        echo "Peptide length: 10" >> "${sample_id}_netmhcpan1_10mer.xls"
        echo "Prediction completed at: \$(date)" >> "${sample_id}_netmhcpan1_10mer.xls"
        echo "netMHCpan version: \$(netMHCpan -v 2>&1)" >> "${sample_id}_netmhcpan1_10mer.xls"
    else
        echo "# No 10mer peptides found in ${peptide_input_f} for sample ${sample_id}" > "${sample_id}_netmhcpan1_10mer.xls"
        echo "HLA alleles used: \$hla_alleles" >> "${sample_id}_netmhcpan1_10mer.xls"
    fi
    """
}

process NETMHCPAN1_PREDICT_11MER {
    tag "${sample_id} - 11mer"
    label 'process_medium'
    
    publishDir(
        path: { 
            def base_sample_id = sample_id.split('_')[0]
            return "${params.outdir}/${base_sample_id}/epitope_prediction/netmhcpan1/11mer" 
        },
        mode: 'copy',
        overwrite: true
    )
    
    input:
    tuple val(sample_id), path(peptide_input_f), path(hla_input_f)
    
    output:
    tuple val(sample_id), path("${sample_id}_netmhcpan1_11mer.xls"), emit: predictions
    
    script:
    """
    #!/bin/bash
    cp "${peptide_input_f}" "peptides.fasta"
    
    # 读取HLA等位基因并移除星号
    hla_alleles=\$(cat "${hla_input_f}" | grep -v "^#" | head -n 1 | sed 's/\\*//g')
    
    # 运行netMHCpan预测
    if [ -s "peptides.fasta" ]; then
        netMHCpan \\
            -f "peptides.fasta" \\
            -a "\$hla_alleles" \\
            -l 11 \\
            -BA \\
            -xls \\
            -xlsfile "${sample_id}_netmhcpan1_11mer.xls"
        
        # 添加预测信息
        echo "Peptide length: 11" >> "${sample_id}_netmhcpan1_11mer.xls"
        echo "Prediction completed at: \$(date)" >> "${sample_id}_netmhcpan1_11mer.xls"
        echo "netMHCpan version: \$(netMHCpan -v 2>&1)" >> "${sample_id}_netmhcpan1_11mer.xls"
    else
        echo "# No 11mer peptides found in ${peptide_input_f} for sample ${sample_id}" > "${sample_id}_netmhcpan1_11mer.xls"
        echo "HLA alleles used: \$hla_alleles" >> "${sample_id}_netmhcpan1_11mer.xls"
    fi
    """
}
