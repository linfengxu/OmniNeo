#!/usr/bin/env nextflow

workflow NETMHCPAN2_WORKFLOW {
    take:
    peptides_12  // Channel: tuple val(sample_id), path(peptide_file)
    peptides_13  // Channel: tuple val(sample_id), path(peptide_file)
    peptides_14  // Channel: tuple val(sample_id), path(peptide_file)
    peptides_15  // Channel: tuple val(sample_id), path(peptide_file)
    peptides_16  // Channel: tuple val(sample_id), path(peptide_file)
    peptides_17  // Channel: tuple val(sample_id), path(peptide_file)
    peptides_18  // Channel: tuple val(sample_id), path(peptide_file)
    peptides_19  // Channel: tuple val(sample_id), path(peptide_file)
    peptides_20  // Channel: tuple val(sample_id), path(peptide_file)
    peptides_21  // Channel: tuple val(sample_id), path(peptide_file)
    peptides_22  // Channel: tuple val(sample_id), path(peptide_file)
    peptides_23  // Channel: tuple val(sample_id), path(peptide_file)
    peptides_24  // Channel: tuple val(sample_id), path(peptide_file)
    peptides_25  // Channel: tuple val(sample_id), path(peptide_file)
    hla_data    // Channel: tuple val(sample_id), path(hla_file)

    main:
    log.info "[NETMHCPAN2] Starting workflow"

    // 规范化sample_id
    def normalized_12mer = peptides_12.map { sample_id, file ->
        def norm_id = sample_id.replaceAll('_(dna|rna)_(normal|tumor)$', '')
        tuple(norm_id, file)
    }
    def normalized_13mer = peptides_13.map { sample_id, file ->
        def norm_id = sample_id.replaceAll('_(dna|rna)_(normal|tumor)$', '')
        tuple(norm_id, file)
    }
    def normalized_14mer = peptides_14.map { sample_id, file ->
        def norm_id = sample_id.replaceAll('_(dna|rna)_(normal|tumor)$', '')
        tuple(norm_id, file)
    }
    def normalized_15mer = peptides_15.map { sample_id, file ->
        def norm_id = sample_id.replaceAll('_(dna|rna)_(normal|tumor)$', '')
        tuple(norm_id, file)
    }
    def normalized_16mer = peptides_16.map { sample_id, file ->
        def norm_id = sample_id.replaceAll('_(dna|rna)_(normal|tumor)$', '')
        tuple(norm_id, file)
    }
    def normalized_17mer = peptides_17.map { sample_id, file ->
        def norm_id = sample_id.replaceAll('_(dna|rna)_(normal|tumor)$', '')
        tuple(norm_id, file)
    }
    def normalized_18mer = peptides_18.map { sample_id, file ->
        def norm_id = sample_id.replaceAll('_(dna|rna)_(normal|tumor)$', '')
        tuple(norm_id, file)
    }
    def normalized_19mer = peptides_19.map { sample_id, file ->
        def norm_id = sample_id.replaceAll('_(dna|rna)_(normal|tumor)$', '')
        tuple(norm_id, file)
    }
    def normalized_20mer = peptides_20.map { sample_id, file ->
        def norm_id = sample_id.replaceAll('_(dna|rna)_(normal|tumor)$', '')
        tuple(norm_id, file)
    }
    def normalized_21mer = peptides_21.map { sample_id, file ->
        def norm_id = sample_id.replaceAll('_(dna|rna)_(normal|tumor)$', '')
        tuple(norm_id, file)
    }
    def normalized_22mer = peptides_22.map { sample_id, file ->
        def norm_id = sample_id.replaceAll('_(dna|rna)_(normal|tumor)$', '')
        tuple(norm_id, file)
    }
    def normalized_23mer = peptides_23.map { sample_id, file ->
        def norm_id = sample_id.replaceAll('_(dna|rna)_(normal|tumor)$', '')
        tuple(norm_id, file)
    }
    def normalized_24mer = peptides_24.map { sample_id, file ->
        def norm_id = sample_id.replaceAll('_(dna|rna)_(normal|tumor)$', '')
        tuple(norm_id, file)
    }
    def normalized_25mer = peptides_25.map { sample_id, file ->
        def norm_id = sample_id.replaceAll('_(dna|rna)_(normal|tumor)$', '')
        tuple(norm_id, file)
    }
    def normalized_hla = hla_data.map { sample_id, file ->
        def norm_id = sample_id.replaceAll('_(dna|rna)_(normal|tumor)$', '')
        tuple(norm_id, file)
    }

    // 并行运行不同长度的预测
    predict_12mer = NETMHCPAN2_PREDICT_12MER(
        normalized_12mer.combine(normalized_hla, by: 0)
    )
    predict_13mer = NETMHCPAN2_PREDICT_13MER(
        normalized_13mer.combine(normalized_hla, by: 0)
    )
    predict_14mer = NETMHCPAN2_PREDICT_14MER(
        normalized_14mer.combine(normalized_hla, by: 0)
    )
    predict_15mer = NETMHCPAN2_PREDICT_15MER(
        normalized_15mer.combine(normalized_hla, by: 0)
    )
    predict_16mer = NETMHCPAN2_PREDICT_16MER(
        normalized_16mer.combine(normalized_hla, by: 0)
    )
    predict_17mer = NETMHCPAN2_PREDICT_17MER(
        normalized_17mer.combine(normalized_hla, by: 0)
    )
    predict_18mer = NETMHCPAN2_PREDICT_18MER(
        normalized_18mer.combine(normalized_hla, by: 0)
    )
    predict_19mer = NETMHCPAN2_PREDICT_19MER(
        normalized_19mer.combine(normalized_hla, by: 0)
    )
    predict_20mer = NETMHCPAN2_PREDICT_20MER(
        normalized_20mer.combine(normalized_hla, by: 0)
    )
    predict_21mer = NETMHCPAN2_PREDICT_21MER(
        normalized_21mer.combine(normalized_hla, by: 0)
    )
    predict_22mer = NETMHCPAN2_PREDICT_22MER(
        normalized_22mer.combine(normalized_hla, by: 0)
    )
    predict_23mer = NETMHCPAN2_PREDICT_23MER(
        normalized_23mer.combine(normalized_hla, by: 0)
    )
    predict_24mer = NETMHCPAN2_PREDICT_24MER(
        normalized_24mer.combine(normalized_hla, by: 0)
    )
    predict_25mer = NETMHCPAN2_PREDICT_25MER(
        normalized_25mer.combine(normalized_hla, by: 0)
    )

    emit:
    predictions_12mer = predict_12mer.predictions
    predictions_13mer = predict_13mer.predictions
    predictions_14mer = predict_14mer.predictions
    predictions_15mer = predict_15mer.predictions
    predictions_16mer = predict_16mer.predictions
    predictions_17mer = predict_17mer.predictions
    predictions_18mer = predict_18mer.predictions
    predictions_19mer = predict_19mer.predictions
    predictions_20mer = predict_20mer.predictions
    predictions_21mer = predict_21mer.predictions
    predictions_22mer = predict_22mer.predictions
    predictions_23mer = predict_23mer.predictions
    predictions_24mer = predict_24mer.predictions
    predictions_25mer = predict_25mer.predictions
}

process NETMHCPAN2_PREDICT_12MER {
    tag "${sample_id} - 12mer"
    label 'process_medium'
    debug true
    
    publishDir(
        path: { 
            def base_sample_id = sample_id.split('_')[0]
            return "${params.outdir}/${base_sample_id}/epitope_prediction/netmhcpan2/12mer" 
        },
        mode: 'copy',
        overwrite: true
    )
    
    input:
    tuple val(sample_id), path(peptide_input_f), path(hla_input_f)
    
    output:
    tuple val(sample_id), path("${sample_id}_netmhcpan2_12mer.xls"), emit: predictions
    
    script:
    """
    #!/bin/bash
    
    
    echo "[DEBUG] Processing sample: ${sample_id}"
    echo "[DEBUG] Peptide file: ${peptide_input_f}"
    echo "[DEBUG] HLA file: ${hla_input_f}"
    
    # 检查输入文件
    if [ ! -s "${peptide_input_f}" ]; then
        echo "[ERROR] Peptide file is empty or does not exist: ${peptide_input_f}"
        exit 1
    fi
    if [ ! -s "${hla_input_f}" ]; then
        echo "[ERROR] HLA file is empty or does not exist: ${hla_input_f}"
        exit 1
    fi
    
    # 显示当前目录和文件信息
    echo "[DEBUG] Current directory: \$(pwd)"
    echo "[DEBUG] Directory contents before copy:"
    ls -la
    
    # 复制输入文件到工作目录
    echo "[DEBUG] Copying peptide file..."
    cp -v "${peptide_input_f}" "peptides.fasta"
    if [ ! -s "peptides.fasta" ]; then
        echo "[ERROR] Failed to copy peptide file"
        exit 1
    fi
    
    echo "[DEBUG] Directory contents after copy:"
    ls -la
    
    # 读取HLA等位基因并移除星号
    hla_alleles=\$(cat "${hla_input_f}" \\
        | grep -v "^#" \\
        | head -n 1 \\
        | sed 's/\\*//g')
    echo "[DEBUG] HLA alleles: \${hla_alleles}"
    
    # netMHCIIpan
    if [ -s "peptides.fasta" ]; then
        echo "[DEBUG] Running netMHCIIpan prediction..."
        # 确保文件存在且可读
        ls -l "peptides.fasta"
        cat "peptides.fasta" | head -n 5
        netMHCIIpan \\
            -f "peptides.fasta" \\
            -a "\$hla_alleles" \\
            -l 12 \\
            -BA \\
            -xls \\
            -xlsfile "${sample_id}_netmhcpan2_12mer.xls" || {
            echo "[ERROR] netMHCpan prediction failed"
            exit 1
        }
        
        # 添加预测信息
        echo "Peptide length: 12" >> "${sample_id}_netmhcpan2_12mer.xls"
        echo "Prediction completed at: \$(date)" >> "${sample_id}_netmhcpan2_12mer.xls"
        echo "netMHCpan version: \$(netMHCIIpan -v 2>&1)" >> "${sample_id}_netmhcpan2_12mer.xls"
        
        echo "[DEBUG] Prediction completed"
    else
        echo "[DEBUG] No peptides found in input file"
        echo "# No 12mer peptides found in ${peptide_input_f} for sample ${sample_id}" > "${sample_id}_netmhcpan2_12mer.xls"
        echo "HLA alleles used: \$hla_alleles" >> "${sample_id}_netmhcpan2_12mer.xls"
    fi
    """
}

// 为13-25mer创建类似的进程，只需要修改mer的数字
process NETMHCPAN2_PREDICT_13MER {
    tag "${sample_id} - 13mer"
    label 'process_medium'
    
    publishDir(
        path: { 
            def base_sample_id = sample_id.split('_')[0]
            return "${params.outdir}/${base_sample_id}/epitope_prediction/netmhcpan2/13mer" 
        },
        mode: 'copy',
        overwrite: true
    )
    
    input:
    tuple val(sample_id), path(peptide_input_f), path(hla_input_f)
    
    output:
    tuple val(sample_id), path("${sample_id}_netmhcpan2_13mer.xls"), emit: predictions
    
    script:
    """
    #!/bin/bash
    
    
    # 显示当前目录和文件信息
    echo "[DEBUG] Current directory: \$(pwd)"
    echo "[DEBUG] Directory contents before copy:"
    ls -la
    
    # 复制输入文件到工作目录
    echo "[DEBUG] Copying peptide file..."
    cp -v "${peptide_input_f}" "peptides.fasta"
    if [ ! -s "peptides.fasta" ]; then
        echo "[ERROR] Failed to copy peptide file"
        exit 1
    fi
    
    echo "[DEBUG] Directory contents after copy:"
    ls -la
    
    # 读取HLA等位基因并移除星号
    hla_alleles=\$(cat "${hla_input_f}" \\
        | grep -v "^#" \\
        | head -n 1 \\
        | sed 's/\\*//g')
    
    # netMHCIIpan
    if [ -s "peptides.fasta" ]; then
        # 确保文件存在且可读
        ls -l "peptides.fasta"
        cat "peptides.fasta" | head -n 5
        netMHCIIpan \\
            -f "peptides.fasta" \\
            -a "\$hla_alleles" \\
            -l 13 \\
            -BA \\
            -xls \\
            -xlsfile "${sample_id}_netmhcpan2_13mer.xls" || {
            echo "[ERROR] netMHCIIpan prediction failed"
            exit 1
        }
        
        # 添加预测信息
        echo "Peptide length: 13" >> "${sample_id}_netmhcpan2_13mer.xls"
        echo "Prediction completed at: \$(date)" >> "${sample_id}_netmhcpan2_13mer.xls"
        echo "netMHCIIpan version: \$(netMHCIIpan -v 2>&1)" >> "${sample_id}_netmhcpan2_13mer.xls"
    else
        echo "# No 13mer peptides found in ${peptide_input_f} for sample ${sample_id}" > "${sample_id}_netmhcpan2_13mer.xls"
        echo "HLA alleles used: \$hla_alleles" >> "${sample_id}_netmhcpan2_13mer.xls"
    fi
    """
}

// 14mer process
process NETMHCPAN2_PREDICT_14MER {
    tag "${sample_id} - 14mer"
    label 'process_medium'
    
    publishDir(
        path: { 
            def base_sample_id = sample_id.split('_')[0]
            return "${params.outdir}/${base_sample_id}/epitope_prediction/netmhcpan2/14mer" 
        },
        mode: 'copy',
        overwrite: true
    )
    
    input:
    tuple val(sample_id), path(peptide_input_f), path(hla_input_f)
    
    output:
    tuple val(sample_id), path("${sample_id}_netmhcpan2_14mer.xls"), emit: predictions
    
    script:
    """
    #!/bin/bash
    
    
    # 显示当前目录和文件信息
    echo "[DEBUG] Current directory: \$(pwd)"
    echo "[DEBUG] Directory contents before copy:"
    ls -la
    
    # 复制输入文件到工作目录
    echo "[DEBUG] Copying peptide file..."
    cp -v "${peptide_input_f}" "peptides.fasta"
    if [ ! -s "peptides.fasta" ]; then
        echo "[ERROR] Failed to copy peptide file"
        exit 1
    fi
    
    echo "[DEBUG] Directory contents after copy:"
    ls -la
    
    # 读取HLA等位基因并移除星号
    hla_alleles=\$(cat "${hla_input_f}" \\
        | grep -v "^#" \\
        | head -n 1 \\
        | sed 's/\\*//g')
    
    # netMHCIIpan
    if [ -s "peptides.fasta" ]; then
        # 确保文件存在且可读
        ls -l "peptides.fasta"
        cat "peptides.fasta" | head -n 5
        netMHCIIpan \\
            -f "peptides.fasta" \\
            -a "\$hla_alleles" \\
            -l 14 \\
            -BA \\
            -xls \\
            -xlsfile "${sample_id}_netmhcpan2_14mer.xls" || {
            echo "[ERROR] netMHCpan prediction failed"
            exit 1
        }
        
        # 添加预测信息
        echo "Peptide length: 14" >> "${sample_id}_netmhcpan2_14mer.xls"
        echo "Prediction completed at: \$(date)" >> "${sample_id}_netmhcpan2_14mer.xls"
        echo "netMHCIIpan version: \$(netMHCIIpan -v 2>&1)" >> "${sample_id}_netmhcpan2_14mer.xls"
    else
        echo "# No 14mer peptides found in ${peptide_input_f} for sample ${sample_id}" > "${sample_id}_netmhcpan2_14mer.xls"
        echo "HLA alleles used: \$hla_alleles" >> "${sample_id}_netmhcpan2_14mer.xls"
    fi
    """
}

// 15mer process
process NETMHCPAN2_PREDICT_15MER {
    tag "${sample_id} - 15mer"
    label 'process_medium'
    
    publishDir(
        path: { 
            def base_sample_id = sample_id.split('_')[0]
            return "${params.outdir}/${base_sample_id}/epitope_prediction/netmhcpan2/15mer" 
        },
        mode: 'copy',
        overwrite: true
    )
    
    input:
    tuple val(sample_id), path(peptide_input_f), path(hla_input_f)
    
    output:
    tuple val(sample_id), path("${sample_id}_netmhcpan2_15mer.xls"), emit: predictions
    
    script:
    """
    #!/bin/bash
    
    
    # 显示当前目录和文件信息
    echo "[DEBUG] Current directory: \$(pwd)"
    echo "[DEBUG] Directory contents before copy:"
    ls -la
    
    # 复制输入文件到工作目录
    echo "[DEBUG] Copying peptide file..."
    cp -v "${peptide_input_f}" "peptides.fasta"
    if [ ! -s "peptides.fasta" ]; then
        echo "[ERROR] Failed to copy peptide file"
        exit 1
    fi
    
    echo "[DEBUG] Directory contents after copy:"
    ls -la
    
    # 读取HLA等位基因并移除星号
    hla_alleles=\$(cat "${hla_input_f}" \\
        | grep -v "^#" \\
        | head -n 1 \\
        | sed 's/\\*//g')
    
    # netMHCIIpan
    if [ -s "peptides.fasta" ]; then
        # 确保文件存在且可读
        ls -l "peptides.fasta"
        cat "peptides.fasta" | head -n 5
        netMHCIIpan \\
            -f "peptides.fasta" \\
            -a "\$hla_alleles" \\
            -l 15 \\
            -BA \\
            -xls \\
            -xlsfile "${sample_id}_netmhcpan2_15mer.xls" || {
            echo "[ERROR] netMHCpan prediction failed"
            exit 1
        }
        
        # 添加预测信息
        echo "Peptide length: 15" >> "${sample_id}_netmhcpan2_15mer.xls"
        echo "Prediction completed at: \$(date)" >> "${sample_id}_netmhcpan2_15mer.xls"
        echo "netMHCIIpan version: \$(netMHCIIpan -v 2>&1)" >> "${sample_id}_netmhcpan2_15mer.xls"
    else
        echo "# No 15mer peptides found in ${peptide_input_f} for sample ${sample_id}" > "${sample_id}_netmhcpan2_15mer.xls"
        echo "HLA alleles used: \$hla_alleles" >> "${sample_id}_netmhcpan2_15mer.xls"
    fi
    """
}

// 16mer process
process NETMHCPAN2_PREDICT_16MER {
    tag "${sample_id} - 16mer"
    label 'process_medium'
    
    publishDir(
        path: { 
            def base_sample_id = sample_id.split('_')[0]
            return "${params.outdir}/${base_sample_id}/epitope_prediction/netmhcpan2/16mer" 
        },
        mode: 'copy',
        overwrite: true
    )
    
    input:
    tuple val(sample_id), path(peptide_input_f), path(hla_input_f)
    
    output:
    tuple val(sample_id), path("${sample_id}_netmhcpan2_16mer.xls"), emit: predictions
    
    script:
    """
    #!/bin/bash
    
    
    # 显示当前目录和文件信息
    echo "[DEBUG] Current directory: \$(pwd)"
    echo "[DEBUG] Directory contents before copy:"
    ls -la
    
    # 复制输入文件到工作目录
    echo "[DEBUG] Copying peptide file..."
    cp -v "${peptide_input_f}" "peptides.fasta"
    if [ ! -s "peptides.fasta" ]; then
        echo "[ERROR] Failed to copy peptide file"
        exit 1
    fi
    
    echo "[DEBUG] Directory contents after copy:"
    ls -la
    
    # 读取HLA等位基因并移除星号
    hla_alleles=\$(cat "${hla_input_f}" \\
        | grep -v "^#" \\
        | head -n 1 \\
        | sed 's/\\*//g')
    
    # netMHCIIpan
    if [ -s "peptides.fasta" ]; then
        # 确保文件存在且可读
        ls -l "peptides.fasta"
        cat "peptides.fasta" | head -n 5
        netMHCIIpan \\
            -f "peptides.fasta" \\
            -a "\$hla_alleles" \\
            -l 16 \\
            -BA \\
            -xls \\
            -xlsfile "${sample_id}_netmhcpan2_16mer.xls" || {
            echo "[ERROR] netMHCIIpan prediction failed"
            exit 1
        }
        
        # 添加预测信息
        echo "Peptide length: 16" >> "${sample_id}_netmhcpan2_16mer.xls"
        echo "Prediction completed at: \$(date)" >> "${sample_id}_netmhcpan2_16mer.xls"
        echo "netMHCIIpan version: \$(netMHCIIpan -v 2>&1)" >> "${sample_id}_netmhcpan2_16mer.xls"
    else
        echo "# No 16mer peptides found in ${peptide_input_f} for sample ${sample_id}" > "${sample_id}_netmhcpan2_16mer.xls"
        echo "HLA alleles used: \$hla_alleles" >> "${sample_id}_netmhcpan2_16mer.xls"
    fi
    """
}

// 17mer process
process NETMHCPAN2_PREDICT_17MER {
    tag "${sample_id} - 17mer"
    label 'process_medium'
    
    publishDir(
        path: { 
            def base_sample_id = sample_id.split('_')[0]
            return "${params.outdir}/${base_sample_id}/epitope_prediction/netmhcpan2/17mer" 
        },
        mode: 'copy',
        overwrite: true
    )
    
    input:
    tuple val(sample_id), path(peptide_input_f), path(hla_input_f)
    
    output:
    tuple val(sample_id), path("${sample_id}_netmhcpan2_17mer.xls"), emit: predictions
    
    script:
    """
    #!/bin/bash
    
    
    # 显示当前目录和文件信息
    echo "[DEBUG] Current directory: \$(pwd)"
    echo "[DEBUG] Directory contents before copy:"
    ls -la
    
    # 复制输入文件到工作目录
    echo "[DEBUG] Copying peptide file..."
    cp -v "${peptide_input_f}" "peptides.fasta"
    if [ ! -s "peptides.fasta" ]; then
        echo "[ERROR] Failed to copy peptide file"
        exit 1
    fi
    
    echo "[DEBUG] Directory contents after copy:"
    ls -la
    
    # 读取HLA等位基因并移除星号
    hla_alleles=\$(cat "${hla_input_f}" \\
        | grep -v "^#" \\
        | head -n 1 \\
        | sed 's/\\*//g')
    
    # netMHCIIpan
    if [ -s "peptides.fasta" ]; then
        # 确保文件存在且可读
        ls -l "peptides.fasta"
        cat "peptides.fasta" | head -n 5
        netMHCIIpan \\
            -f "peptides.fasta" \\
            -a "\$hla_alleles" \\
            -l 17 \\
            -BA \\
            -xls \\
            -xlsfile "${sample_id}_netmhcpan2_17mer.xls" || {
            echo "[ERROR] netMHCpan prediction failed"
            exit 1
        }
        
        # 添加预测信息
        echo "Peptide length: 17" >> "${sample_id}_netmhcpan2_17mer.xls"
        echo "Prediction completed at: \$(date)" >> "${sample_id}_netmhcpan2_17mer.xls"
        echo "netMHCIIpan version: \$(netMHCIIpan -v 2>&1)" >> "${sample_id}_netmhcpan2_17mer.xls"
    else
        echo "# No 17mer peptides found in ${peptide_input_f} for sample ${sample_id}" > "${sample_id}_netmhcpan2_17mer.xls"
        echo "HLA alleles used: \$hla_alleles" >> "${sample_id}_netmhcpan2_17mer.xls"
    fi
    """
}

// 18mer process
process NETMHCPAN2_PREDICT_18MER {
    tag "${sample_id} - 18mer"
    label 'process_medium'
    
    publishDir(
        path: { 
            def base_sample_id = sample_id.split('_')[0]
            return "${params.outdir}/${base_sample_id}/epitope_prediction/netmhcpan2/18mer" 
        },
        mode: 'copy',
        overwrite: true
    )
    
    input:
    tuple val(sample_id), path(peptide_input_f), path(hla_input_f)
    
    output:
    tuple val(sample_id), path("${sample_id}_netmhcpan2_18mer.xls"), emit: predictions
    
    script:
    """
    #!/bin/bash
    
    
    # 显示当前目录和文件信息
    echo "[DEBUG] Current directory: \$(pwd)"
    echo "[DEBUG] Directory contents before copy:"
    ls -la
    
    # 复制输入文件到工作目录
    echo "[DEBUG] Copying peptide file..."
    cp -v "${peptide_input_f}" "peptides.fasta"
    if [ ! -s "peptides.fasta" ]; then
        echo "[ERROR] Failed to copy peptide file"
        exit 1
    fi
    
    echo "[DEBUG] Directory contents after copy:"
    ls -la
    
    # 读取HLA等位基因并移除星号
    hla_alleles=\$(cat "${hla_input_f}" \\
        | grep -v "^#" \\
        | head -n 1 \\
        | sed 's/\\*//g')
    
    # 运行netMHCpan预测
    if [ -s "peptides.fasta" ]; then
        # 确保文件存在且可读
        ls -l "peptides.fasta"
        cat "peptides.fasta" | head -n 5
        netMHCIIpan \\
            -f "peptides.fasta" \\
            -a "\$hla_alleles" \\
            -l 18 \\
            -BA \\
            -xls \\
            -xlsfile "${sample_id}_netmhcpan2_18mer.xls" || {
            echo "[ERROR] netMHCpan prediction failed"
            exit 1
        }
        
        # 添加预测信息
        echo "Peptide length: 18" >> "${sample_id}_netmhcpan2_18mer.xls"
        echo "Prediction completed at: \$(date)" >> "${sample_id}_netmhcpan2_18mer.xls"
        echo "netMHnetMHCIIpanpan version: \$(netMHCIIpan -v 2>&1)" >> "${sample_id}_netmhcpan2_18mer.xls"
    else
        echo "# No 18mer peptides found in ${peptide_input_f} for sample ${sample_id}" > "${sample_id}_netmhcpan2_18mer.xls"
        echo "HLA alleles used: \$hla_alleles" >> "${sample_id}_netmhcpan2_18mer.xls"
    fi
    """
}

// 19mer process
process NETMHCPAN2_PREDICT_19MER {
    tag "${sample_id} - 19mer"
    label 'process_medium'
    
    publishDir(
        path: { 
            def base_sample_id = sample_id.split('_')[0]
            return "${params.outdir}/${base_sample_id}/epitope_prediction/netmhcpan2/19mer" 
        },
        mode: 'copy',
        overwrite: true
    )
    
    input:
    tuple val(sample_id), path(peptide_input_f), path(hla_input_f)
    
    output:
    tuple val(sample_id), path("${sample_id}_netmhcpan2_19mer.xls"), emit: predictions
    
    script:
    """
    #!/bin/bash
    
    
    # 显示当前目录和文件信息
    echo "[DEBUG] Current directory: \$(pwd)"
    echo "[DEBUG] Directory contents before copy:"
    ls -la
    
    # 复制输入文件到工作目录
    echo "[DEBUG] Copying peptide file..."
    cp -v "${peptide_input_f}" "peptides.fasta"
    if [ ! -s "peptides.fasta" ]; then
        echo "[ERROR] Failed to copy peptide file"
        exit 1
    fi
    
    echo "[DEBUG] Directory contents after copy:"
    ls -la
    
    # 读取HLA等位基因并移除星号
    hla_alleles=\$(cat "${hla_input_f}" \\
        | grep -v "^#" \\
        | head -n 1 \\
        | sed 's/\\*//g')
    
    # netMHCIIpan
    if [ -s "peptides.fasta" ]; then
        # 确保文件存在且可读
        ls -l "peptides.fasta"
        cat "peptides.fasta" | head -n 5
        netMHCIIpan \\
            -f "peptides.fasta" \\
            -a "\$hla_alleles" \\
            -l 19 \\
            -BA \\
            -xls \\
            -xlsfile "${sample_id}_netmhcpan2_19mer.xls" || {
            echo "[ERROR] netMHCIIpan prediction failed"
            exit 1
        }
        
        # 添加预测信息
        echo "Peptide length: 19" >> "${sample_id}_netmhcpan2_19mer.xls"
        echo "Prediction completed at: \$(date)" >> "${sample_id}_netmhcpan2_19mer.xls"
        echo "netMHCIIpan version: \$(netMHCIIpan -v 2>&1)" >> "${sample_id}_netmhcpan2_19mer.xls"
    else
        echo "# No 19mer peptides found in ${peptide_input_f} for sample ${sample_id}" > "${sample_id}_netmhcpan2_19mer.xls"
        echo "HLA alleles used: \$hla_alleles" >> "${sample_id}_netmhcpan2_19mer.xls"
    fi
    """
}

// 20mer process
process NETMHCPAN2_PREDICT_20MER {
    tag "${sample_id} - 20mer"
    label 'process_medium'
    
    publishDir(
        path: { 
            def base_sample_id = sample_id.split('_')[0]
            return "${params.outdir}/${base_sample_id}/epitope_prediction/netmhcpan2/20mer" 
        },
        mode: 'copy',
        overwrite: true
    )
    
    input:
    tuple val(sample_id), path(peptide_input_f), path(hla_input_f)
    
    output:
    tuple val(sample_id), path("${sample_id}_netmhcpan2_20mer.xls"), emit: predictions
    
    script:
    """
    #!/bin/bash
    
    
    # 显示当前目录和文件信息
    echo "[DEBUG] Current directory: \$(pwd)"
    echo "[DEBUG] Directory contents before copy:"
    ls -la
    
    # 复制输入文件到工作目录
    echo "[DEBUG] Copying peptide file..."
    cp -v "${peptide_input_f}" "peptides.fasta"
    if [ ! -s "peptides.fasta" ]; then
        echo "[ERROR] Failed to copy peptide file"
        exit 1
    fi
    
    echo "[DEBUG] Directory contents after copy:"
    ls -la
    
    # 读取HLA等位基因并移除星号
    hla_alleles=\$(cat "${hla_input_f}" \\
        | grep -v "^#" \\
        | head -n 1 \\
        | sed 's/\\*//g')
    
    # netMHCIIpan
    if [ -s "peptides.fasta" ]; then
        # 确保文件存在且可读
        ls -l "peptides.fasta"
        cat "peptides.fasta" | head -n 5
        netMHCIIpan \\
            -f "peptides.fasta" \\
            -a "\$hla_alleles" \\
            -l 20 \\
            -BA \\
            -xls \\
            -xlsfile "${sample_id}_netmhcpan2_20mer.xls" || {
            echo "[ERROR] netMHCpan prediction failed"
            exit 1
        }
        
        # 添加预测信息
        echo "Peptide length: 20" >> "${sample_id}_netmhcpan2_20mer.xls"
        echo "Prediction completed at: \$(date)" >> "${sample_id}_netmhcpan2_20mer.xls"
        echo "netMHCIIpan version: \$(netMHCIIpan -v 2>&1)" >> "${sample_id}_netmhcpan2_20mer.xls"
    else
        echo "# No 20mer peptides found in ${peptide_input_f} for sample ${sample_id}" > "${sample_id}_netmhcpan2_20mer.xls"
        echo "HLA alleles used: \$hla_alleles" >> "${sample_id}_netmhcpan2_20mer.xls"
    fi
    """
}

// 21mer process
process NETMHCPAN2_PREDICT_21MER {
    tag "${sample_id} - 21mer"
    label 'process_medium'
    
    publishDir(
        path: { 
            def base_sample_id = sample_id.split('_')[0]
            return "${params.outdir}/${base_sample_id}/epitope_prediction/netmhcpan2/21mer" 
        },
        mode: 'copy',
        overwrite: true
    )
    
    input:
    tuple val(sample_id), path(peptide_input_f), path(hla_input_f)
    
    output:
    tuple val(sample_id), path("${sample_id}_netmhcpan2_21mer.xls"), emit: predictions
    
    script:
    """
    #!/bin/bash
    
    
    # 显示当前目录和文件信息
    echo "[DEBUG] Current directory: \$(pwd)"
    echo "[DEBUG] Directory contents before copy:"
    ls -la
    
    # 复制输入文件到工作目录
    echo "[DEBUG] Copying peptide file..."
    cp -v "${peptide_input_f}" "peptides.fasta"
    if [ ! -s "peptides.fasta" ]; then
        echo "[ERROR] Failed to copy peptide file"
        exit 1
    fi
    
    echo "[DEBUG] Directory contents after copy:"
    ls -la
    
    # 读取HLA等位基因并移除星号
    hla_alleles=\$(cat "${hla_input_f}" \\
        | grep -v "^#" \\
        | head -n 1 \\
        | sed 's/\\*//g')
    
    # netMHCIIpan
    if [ -s "peptides.fasta" ]; then
        # 确保文件存在且可读
        ls -l "peptides.fasta"
        cat "peptides.fasta" | head -n 5
        netMHCIIpan \\
            -f "peptides.fasta" \\
            -a "\$hla_alleles" \\
            -l 21 \\
            -BA \\
            -xls \\
            -xlsfile "${sample_id}_netmhcpan2_21mer.xls" || {
            echo "[ERROR] netMHCpan prediction failed"
            exit 1
        }
        
        # 添加预测信息
        echo "Peptide length: 21" >> "${sample_id}_netmhcpan2_21mer.xls"
        echo "Prediction completed at: \$(date)" >> "${sample_id}_netmhcpan2_21mer.xls"
        echo "netMHCIIpan version: \$(netMHCIIpan -v 2>&1)" >> "${sample_id}_netmhcpan2_21mer.xls"
    else
        echo "# No 21mer peptides found in ${peptide_input_f} for sample ${sample_id}" > "${sample_id}_netmhcpan2_21mer.xls"
        echo "HLA alleles used: \$hla_alleles" >> "${sample_id}_netmhcpan2_21mer.xls"
    fi
    """
}

// 22mer process
process NETMHCPAN2_PREDICT_22MER {
    tag "${sample_id} - 22mer"
    label 'process_medium'
    
    publishDir(
        path: { 
            def base_sample_id = sample_id.split('_')[0]
            return "${params.outdir}/${base_sample_id}/epitope_prediction/netmhcpan2/22mer" 
        },
        mode: 'copy',
        overwrite: true
    )
    
    input:
    tuple val(sample_id), path(peptide_input_f), path(hla_input_f)
    
    output:
    tuple val(sample_id), path("${sample_id}_netmhcpan2_22mer.xls"), emit: predictions
    
    script:
    """
    #!/bin/bash
    
    
    # 显示当前目录和文件信息
    echo "[DEBUG] Current directory: \$(pwd)"
    echo "[DEBUG] Directory contents before copy:"
    ls -la
    
    # 复制输入文件到工作目录
    echo "[DEBUG] Copying peptide file..."
    cp -v "${peptide_input_f}" "peptides.fasta"
    if [ ! -s "peptides.fasta" ]; then
        echo "[ERROR] Failed to copy peptide file"
        exit 1
    fi
    
    echo "[DEBUG] Directory contents after copy:"
    ls -la
    
    # 读取HLA等位基因并移除星号
    hla_alleles=\$(cat "${hla_input_f}" \\
        | grep -v "^#" \\
        | head -n 1 \\
        | sed 's/\\*//g')
    
    # 运行netMHCpan预测
    if [ -s "peptides.fasta" ]; then
        # 确保文件存在且可读
        ls -l "peptides.fasta"
        cat "peptides.fasta" | head -n 5
        netMHCIIpan \\
            -f "peptides.fasta" \\
            -a "\$hla_alleles" \\
            -l 22 \\
            -BA \\
            -xls \\
            -xlsfile "${sample_id}_netmhcpan2_22mer.xls" || {
            echo "[ERROR] netMHCpan prediction failed"
            exit 1
        }
        
        # 添加预测信息
        echo "Peptide length: 22" >> "${sample_id}_netmhcpan2_22mer.xls"
        echo "Prediction completed at: \$(date)" >> "${sample_id}_netmhcpan2_22mer.xls"
        echo "netMHCpan version: \$(netMHCpan -v 2>&1)" >> "${sample_id}_netmhcpan2_22mer.xls"
    else
        echo "# No 22mer peptides found in ${peptide_input_f} for sample ${sample_id}" > "${sample_id}_netmhcpan2_22mer.xls"
        echo "HLA alleles used: \$hla_alleles" >> "${sample_id}_netmhcpan2_22mer.xls"
    fi
    """
}

// 23mer process
process NETMHCPAN2_PREDICT_23MER {
    tag "${sample_id} - 23mer"
    label 'process_medium'
    
    publishDir(
        path: { 
            def base_sample_id = sample_id.split('_')[0]
            return "${params.outdir}/${base_sample_id}/epitope_prediction/netmhcpan2/23mer" 
        },
        mode: 'copy',
        overwrite: true
    )
    
    input:
    tuple val(sample_id), path(peptide_input_f), path(hla_input_f)
    
    output:
    tuple val(sample_id), path("${sample_id}_netmhcpan2_23mer.xls"), emit: predictions
    
    script:
    """
    #!/bin/bash
    
    
    # 显示当前目录和文件信息
    echo "[DEBUG] Current directory: \$(pwd)"
    echo "[DEBUG] Directory contents before copy:"
    ls -la
    
    # 复制输入文件到工作目录
    echo "[DEBUG] Copying peptide file..."
    cp -v "${peptide_input_f}" "peptides.fasta"
    if [ ! -s "peptides.fasta" ]; then
        echo "[ERROR] Failed to copy peptide file"
        exit 1
    fi
    
    echo "[DEBUG] Directory contents after copy:"
    ls -la
    
    # 读取HLA等位基因并移除星号
    hla_alleles=\$(cat "${hla_input_f}" \\
        | grep -v "^#" \\
        | head -n 1 \\
        | sed 's/\\*//g')
    
    # netMHCIIpan
    if [ -s "peptides.fasta" ]; then
        # 确保文件存在且可读
        ls -l "peptides.fasta"
        cat "peptides.fasta" | head -n 5
        netMHCIIpan \\
            -f "peptides.fasta" \\
            -a "\$hla_alleles" \\
            -l 23 \\
            -BA \\
            -xls \\
            -xlsfile "${sample_id}_netmhcpan2_23mer.xls" || {
            echo "[ERROR] netMHCIIpan prediction failed"
            exit 1
        }
        
        # 添加预测信息
        echo "Peptide length: 23" >> "${sample_id}_netmhcpan2_23mer.xls"
        echo "Prediction completed at: \$(date)" >> "${sample_id}_netmhcpan2_23mer.xls"
        echo "netMHCIIpan version: \$(netMHCIIpan -v 2>&1)" >> "${sample_id}_netmhcpan2_23mer.xls"
    else
        echo "# No 23mer peptides found in ${peptide_input_f} for sample ${sample_id}" > "${sample_id}_netmhcpan2_23mer.xls"
        echo "HLA alleles used: \$hla_alleles" >> "${sample_id}_netmhcpan2_23mer.xls"
    fi
    """
}

// 24mer process
process NETMHCPAN2_PREDICT_24MER {
    tag "${sample_id} - 24mer"
    label 'process_medium'
    
    publishDir(
        path: { 
            def base_sample_id = sample_id.split('_')[0]
            return "${params.outdir}/${base_sample_id}/epitope_prediction/netmhcpan2/24mer" 
        },
        mode: 'copy',
        overwrite: true
    )
    
    input:
    tuple val(sample_id), path(peptide_input_f), path(hla_input_f)
    
    output:
    tuple val(sample_id), path("${sample_id}_netmhcpan2_24mer.xls"), emit: predictions
    
    script:
    """
    #!/bin/bash
    
    
    # 显示当前目录和文件信息
    echo "[DEBUG] Current directory: \$(pwd)"
    echo "[DEBUG] Directory contents before copy:"
    ls -la
    
    # 复制输入文件到工作目录
    echo "[DEBUG] Copying peptide file..."
    cp -v "${peptide_input_f}" "peptides.fasta"
    if [ ! -s "peptides.fasta" ]; then
        echo "[ERROR] Failed to copy peptide file"
        exit 1
    fi
    
    echo "[DEBUG] Directory contents after copy:"
    ls -la
    
    # 读取HLA等位基因并移除星号
    hla_alleles=\$(cat "${hla_input_f}" \\
        | grep -v "^#" \\
        | head -n 1 \\
        | sed 's/\\*//g')
    
    # netMHCIIpan
    if [ -s "peptides.fasta" ]; then
        # 确保文件存在且可读
        ls -l "peptides.fasta"
        cat "peptides.fasta" | head -n 5
        netMHCIIpan \\
            -f "peptides.fasta" \\
            -a "\$hla_alleles" \\
            -l 24 \\
            -BA \\
            -xls \\
            -xlsfile "${sample_id}_netmhcpan2_24mer.xls" || {
            echo "[ERROR] netMHCpan prediction failed"
            exit 1
        }
        
        # 添加预测信息
        echo "Peptide length: 24" >> "${sample_id}_netmhcpan2_24mer.xls"
        echo "Prediction completed at: \$(date)" >> "${sample_id}_netmhcpan2_24mer.xls"
        echo "netMHCIIpan version: \$(netMHCIIpan -v 2>&1)" >> "${sample_id}_netmhcpan2_24mer.xls"
    else
        echo "# No 24mer peptides found in ${peptide_input_f} for sample ${sample_id}" > "${sample_id}_netmhcpan2_24mer.xls"
        echo "HLA alleles used: \$hla_alleles" >> "${sample_id}_netmhcpan2_24mer.xls"
    fi
    """
}

// 25mer process
process NETMHCPAN2_PREDICT_25MER {
    tag "${sample_id} - 25mer"
    label 'process_medium'
    
    publishDir(
        path: { 
            def base_sample_id = sample_id.split('_')[0]
            return "${params.outdir}/${base_sample_id}/epitope_prediction/netmhcpan2/25mer" 
        },
        mode: 'copy',
        overwrite: true
    )
    
    input:
    tuple val(sample_id), path(peptide_input_f), path(hla_input_f)
    
    output:
    tuple val(sample_id), path("${sample_id}_netmhcpan2_25mer.xls"), emit: predictions
    
    script:
    """
    #!/bin/bash
    
    
    # 显示当前目录和文件信息
    echo "[DEBUG] Current directory: \$(pwd)"
    echo "[DEBUG] Directory contents before copy:"
    ls -la
    
    # 复制输入文件到工作目录
    echo "[DEBUG] Copying peptide file..."
    cp -v "${peptide_input_f}" "peptides.fasta"
    if [ ! -s "peptides.fasta" ]; then
        echo "[ERROR] Failed to copy peptide file"
        exit 1
    fi
    
    echo "[DEBUG] Directory contents after copy:"
    ls -la
    
    # 读取HLA等位基因并移除星号
    hla_alleles=\$(cat "${hla_input_f}" \\
        | grep -v "^#" \\
        | head -n 1 \\
        | sed 's/\\*//g')
    
    # netMHCIIpan
    if [ -s "peptides.fasta" ]; then
        # 确保文件存在且可读
        ls -l "peptides.fasta"
        cat "peptides.fasta" | head -n 5
        netMHCIIpan \\
            -f "peptides.fasta" \\
            -a "\$hla_alleles" \\
            -l 25 \\
            -BA \\
            -xls \\
            -xlsfile "${sample_id}_netmhcpan2_25mer.xls" || {
            echo "[ERROR] netMHCIIpan prediction failed"
            exit 1
        }
        
        # 添加预测信息
        echo "Peptide length: 25" >> "${sample_id}_netmhcpan2_25mer.xls"
        echo "Prediction completed at: \$(date)" >> "${sample_id}_netmhcpan2_25mer.xls"
        echo "netMHCIIpan version: \$(netMHCIIpan -v 2>&1)" >> "${sample_id}_netmhcpan2_25mer.xls"
    else
        echo "# No 25mer peptides found in ${peptide_input_f} for sample ${sample_id}" > "${sample_id}_netmhcpan2_25mer.xls"
        echo "HLA alleles used: \$hla_alleles" >> "${sample_id}_netmhcpan2_25mer.xls"
    fi
    """
}
