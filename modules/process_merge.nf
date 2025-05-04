// modules/process_merge.nf

process MERGE_FASTQ {
    tag "${sample_id}:${read_type}"
    label 'process_medium'

    input:
    tuple val(sample_id), path(reads), val(read_type)

    output:
    tuple val(sample_id), path("${sample_id}_${read_type}_merged.fastq.gz"), val(read_type), emit: merged_fastq

    script:
    // If reads is a single file, just create a symlink
    if (reads instanceof Path) {
        """
        if [[ "${reads}" == *.gz ]]; then
            ln -s ${reads} ${sample_id}_${read_type}_merged.fastq.gz
        else
            gzip -c ${reads} > ${sample_id}_${read_type}_merged.fastq.gz
        fi
        """
    }
    // If reads is a list of files, merge them
    else {
        """
        # Create a temporary directory for processing
        mkdir -p tmp_processing

        # Create file list with compression status
        for f in ${reads.join(' ')}; do
            if [[ "\$f" == *.gz ]]; then
                echo "\$f\tgz"
            else
                echo "\$f\traw"
            fi
        done > file_list.txt

        # Get number of cores and use one less for parallel processing
        CORES=\$(nproc)
        THREADS=\$((CORES - 1))
        if [ \$THREADS -lt 1 ]; then
            THREADS=1
        fi

        # Process files in parallel and concatenate
        cat file_list.txt | parallel --colsep '\\t' -j \$THREADS '
            if [[ "{2}" == "gz" ]]; then
                zcat "{1}" > tmp_processing/\$(basename "{1}" .gz).part
            else
                cat "{1}" > tmp_processing/\$(basename "{1}").part
            fi'

        # Concatenate all parts and compress
        cat tmp_processing/*.part | gzip > ${sample_id}_${read_type}_merged.fastq.gz

        # Clean up
        rm -rf tmp_processing file_list.txt
        """
    }
}

process MERGE_PAIRED_FASTQ {
    tag "${sample_id}"
    label 'process_high'
    container { "${params.parallel.container}" }
    
    publishDir {
        def data_type = ""
        def sample_type = ""
        
        if (sample_id.contains("_dna_")) {
            data_type = "dna"
        } else if (sample_id.contains("_rna_")) {
            data_type = "rna"
        }
        
        if (sample_id.contains("_normal")) {
            sample_type = "normal"
        } else if (sample_id.contains("_tumor")) {
            sample_type = "tumor"
        }
        
        return "${params.outdir}/${sample_id.split('_')[0]}/${data_type}/00_merged/${sample_type}"
    },
    mode: 'copy',
    saveAs: { filename -> filename }

    input:
    tuple val(sample_id), path(r1_file), path(r2_file)

    output:
    tuple val(sample_id), path("${sample_id}_R1.fastq.gz"), path("${sample_id}_R2.fastq.gz"), emit: paired_fastq

    script:
    """
    # Link or rename the input files to standardized output names
    ln -sf ${r1_file} ${sample_id}_R1.fastq.gz
    ln -sf ${r2_file} ${sample_id}_R2.fastq.gz
    """
}

workflow MERGE_READS {
    take:
    // Channel structure: [sample_id, dna_normal_r1, dna_normal_r2, dna_tumor_r1, dna_tumor_r2, 
    //                    rna_normal_r1, rna_normal_r2, rna_tumor_r1, rna_tumor_r2]
    reads_ch

    main:
    // 准备8个不同类型的读数通道
    read_channels = Channel.empty()
    
    // DNA Normal R1
    read_channels = read_channels.mix(
        reads_ch.map { id, dn_r1, dn_r2, dt_r1, dt_r2, rn_r1, rn_r2, rt_r1, rt_r2 -> 
            def files = dn_r1 ? dn_r1.tokenize(',') : []
            if (files) {
                return tuple("${id}_dna_normal", files, 'R1')
            } else {
                return null
            }
        }.filter { it != null }
    )
    
    // DNA Normal R2
    read_channels = read_channels.mix(
        reads_ch.map { id, dn_r1, dn_r2, dt_r1, dt_r2, rn_r1, rn_r2, rt_r1, rt_r2 -> 
            def files = dn_r2 ? dn_r2.tokenize(',') : []
            if (files) {
                return tuple("${id}_dna_normal", files, 'R2')
            } else {
                return null
            }
        }.filter { it != null }
    )
    
    // DNA Tumor R1
    read_channels = read_channels.mix(
        reads_ch.map { id, dn_r1, dn_r2, dt_r1, dt_r2, rn_r1, rn_r2, rt_r1, rt_r2 -> 
            def files = dt_r1 ? dt_r1.tokenize(',') : []
            if (files) {
                return tuple("${id}_dna_tumor", files, 'R1')
            } else {
                return null
            }
        }.filter { it != null }
    )
    
    // DNA Tumor R2
    read_channels = read_channels.mix(
        reads_ch.map { id, dn_r1, dn_r2, dt_r1, dt_r2, rn_r1, rn_r2, rt_r1, rt_r2 -> 
            def files = dt_r2 ? dt_r2.tokenize(',') : []
            if (files) {
                return tuple("${id}_dna_tumor", files, 'R2')
            } else {
                return null
            }
        }.filter { it != null }
    )
    
    // RNA Normal R1
    read_channels = read_channels.mix(
        reads_ch.map { id, dn_r1, dn_r2, dt_r1, dt_r2, rn_r1, rn_r2, rt_r1, rt_r2 -> 
            def files = rn_r1 ? rn_r1.tokenize(',') : []
            if (files) {
                return tuple("${id}_rna_normal", files, 'R1')
            } else {
                return null
            }
        }.filter { it != null }
    )
    
    // RNA Normal R2
    read_channels = read_channels.mix(
        reads_ch.map { id, dn_r1, dn_r2, dt_r1, dt_r2, rn_r1, rn_r2, rt_r1, rt_r2 -> 
            def files = rn_r2 ? rn_r2.tokenize(',') : []
            if (files) {
                return tuple("${id}_rna_normal", files, 'R2')
            } else {
                return null
            }
        }.filter { it != null }
    )
    
    // RNA Tumor R1
    read_channels = read_channels.mix(
        reads_ch.map { id, dn_r1, dn_r2, dt_r1, dt_r2, rn_r1, rn_r2, rt_r1, rt_r2 -> 
            def files = rt_r1 ? rt_r1.tokenize(',') : []
            if (files) {
                return tuple("${id}_rna_tumor", files, 'R1')
            } else {
                return null
            }
        }.filter { it != null }
    )
    
    // RNA Tumor R2
    read_channels = read_channels.mix(
        reads_ch.map { id, dn_r1, dn_r2, dt_r1, dt_r2, rn_r1, rn_r2, rt_r1, rt_r2 -> 
            def files = rt_r2 ? rt_r2.tokenize(',') : []
            if (files) {
                return tuple("${id}_rna_tumor", files, 'R2')
            } else {
                return null
            }
        }.filter { it != null }
    )
    
    // 合并所有FASTQ文件
    merged_reads = MERGE_FASTQ(read_channels)
    
    // 组合R1和R2成对的FASTQ文件
    merged_pairs = merged_reads.merged_fastq
        .map { sample_id, file, read_type -> 
            // 移除末尾的_dna_normal, _dna_tumor等
            def base_id = sample_id
            tuple(base_id, read_type, file)
        }
        .groupTuple(by: [0])
        .map { sample_id, read_types, files ->
            def r1 = null
            def r2 = null
            
            for (int i = 0; i < read_types.size(); i++) {
                if (read_types[i] == 'R1') {
                    r1 = files[i]
                } else if (read_types[i] == 'R2') {
                    r2 = files[i]
                }
            }
            
            if (r1 && r2) {
                return tuple(sample_id, r1, r2)
            } else {
                return null
            }
        }
        .filter { it != null }
    
    // 处理最终的成对FASTQ文件
    final_paired_fastq = MERGE_PAIRED_FASTQ(merged_pairs)
    
    // 将结果按照样本类型分开
    sample_types = final_paired_fastq.paired_fastq
        .branch {
            dna_normal: it[0].contains("_dna_normal")
            dna_tumor: it[0].contains("_dna_tumor")
            rna_normal: it[0].contains("_rna_normal")
            rna_tumor: it[0].contains("_rna_tumor")
        }
    
    // 准备输出通道
    dna_normal_paired = sample_types.dna_normal.map { sample_id, r1, r2 ->
        def base_id = sample_id.toString().replaceAll('_dna_normal$', '')
        tuple(base_id, r1, r2)
    }
    
    dna_tumor_paired = sample_types.dna_tumor.map { sample_id, r1, r2 ->
        def base_id = sample_id.toString().replaceAll('_dna_tumor$', '')
        tuple(base_id, r1, r2)
    }
    
    rna_normal_paired = sample_types.rna_normal.map { sample_id, r1, r2 ->
        def base_id = sample_id.toString().replaceAll('_rna_normal$', '')
        tuple(base_id, r1, r2)
    }
    
    rna_tumor_paired = sample_types.rna_tumor.map { sample_id, r1, r2 ->
        def base_id = sample_id.toString().replaceAll('_rna_tumor$', '')
        tuple(base_id, r1, r2)
    }

    emit:
    dna_normal_paired = dna_normal_paired
    dna_tumor_paired = dna_tumor_paired
    rna_normal_paired = rna_normal_paired
    rna_tumor_paired = rna_tumor_paired
}
