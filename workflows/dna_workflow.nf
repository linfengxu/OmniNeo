#!/usr/bin/env nextflow

// DNA/WES analysis workflow

// import modules
// MERGE_READS
include { MERGE_READS } from '../modules/process_merge'
// FASTP
include { FASTP } from '../modules/fastp'
// bwa
include { BWA_INDEX } from '../modules/bwa'
include { BWA } from '../modules/bwa'
//samtools
include { SAMTOOLS_FIXMATE } from '../modules/samtools'
include { SAMTOOLS_SORT } from '../modules/samtools'
include { SAMTOOLS_INDEX } from '../modules/samtools'
// GATK
include { GATK_MARKDUPLICATES } from '../modules/gatk'
include { GATK_ADD_OR_REPLACE_READ_GROUPS } from '../modules/gatk'
include { GATK_BASE_RECALIBRATOR } from '../modules/gatk'
include { GATK_APPLY_BQSR } from '../modules/gatk'
include { GATK_MUTECT2_NORMAL } from '../modules/gatk'
include { GATK_MUTECT2_PAIRED } from '../modules/gatk'
include { GATK_FILTER_MUTECT_CALLS } from '../modules/gatk'
include { BCFTOOLS_FILTER_PASS } from '../modules/gatk'
// ANNOVAR
include { ANNOVAR } from '../modules/annovar'

// Picard
include { PICARD_ADD_OR_REPLACE_READ_GROUPS } from '../modules/picard'

workflow DNA_WORKFLOW {
    take:
    raw_samples_ch  // original sample channel [sample_id, dna_normal_r1, dna_normal_r2, dna_tumor_r1, dna_tumor_r2, rna_normal_r1, rna_normal_r2, rna_tumor_r1, rna_tumor_r2]
    reference_file  // reference genome file

    main:
    // prepare input channel only containing DNA data
    dna_samples_ch = raw_samples_ch.map { id, dn_r1, dn_r2, dt_r1, dt_r2, rn_r1, rn_r2, rt_r1, rt_r2 ->
        // replace empty strings with null to allow process_merge to handle correctly
        def nullify = { str -> str == '' ? null : str }
        tuple(id, 
              nullify(dn_r1), nullify(dn_r2), nullify(dt_r1), nullify(dt_r2),  // DNA data
              null, null, null, null) // RNA data set to null
    }
    
    // first use MERGE_READS module to merge all split files, but only process DNA data
    merged_reads = MERGE_READS(dna_samples_ch)
    
    // extract DNA data from merged reads
    dna_normal_reads = merged_reads.dna_normal_paired
    dna_tumor_reads = merged_reads.dna_tumor_paired
    
    // prepare input channel for FASTP to split samples
    dna_split_samples = dna_normal_reads
        .map { id, r1, r2 -> tuple("${id}_dna_normal", r1, r2) }
        .mix(
            dna_tumor_reads.map { id, r1, r2 -> tuple("${id}_dna_tumor", r1, r2) }
        )
    
    // prepare input channel for paired analysis
    dna_paired_samples = dna_normal_reads
        .map { id, r1, r2 -> tuple(id, ['normal': [r1, r2]]) }
        .join(
            dna_tumor_reads.map { id, r1, r2 -> tuple(id, ['tumor': [r1, r2]]) }
        )
        .map { id, normal, tumor -> tuple(id, normal + tumor) }
    
    // call FASTP module to process original reads
    fastp_result = FASTP(dna_split_samples)
    
    // extract paired reads channel after quality control from FASTP module
    cleaned_reads_ch = fastp_result.cleaned_reads
    
    // create BWA index
    bwa_index_result = BWA_INDEX(reference_file)

    // pass the output of FASTP and the reference genome and index file to the BWA module for alignment
    bwa_result = BWA(
        cleaned_reads_ch,        // cleaned reads
        reference_file,          // reference genome file
        bwa_index_result.index   // index file
    )

    // pass the output of BWA and the reference genome and index file to the SAMTOOLS_FIXMATE module
    fixmate_result = SAMTOOLS_FIXMATE(bwa_result)

    // pass the output of SAMTOOLS_FIXMATE and the reference genome and index file to the SAMTOOLS_SORT module
    sort_result = SAMTOOLS_SORT(fixmate_result)
    
    // pass the output of SAMTOOLS_SORT and the reference genome and index file to the GATK_MARKDUPLICATES module
    markdup_result = GATK_MARKDUPLICATES(sort_result)
    
    // pass the output of GATK_MARKDUPLICATES and the reference genome and index file to the GATK_ADD_OR_REPLACE_READ_GROUPS module
    readgroups_result = GATK_ADD_OR_REPLACE_READ_GROUPS(markdup_result.marked_bam)
    
    // pass the output of GATK_ADD_OR_REPLACE_READ_GROUPS and the reference genome and index file to the GATK_BASE_RECALIBRATOR module
    recal_result = GATK_BASE_RECALIBRATOR(
        readgroups_result.bam_with_read_groups,
        reference_file
    )
    
    // pass the output of GATK_BASE_RECALIBRATOR and the reference genome and index file to the GATK_APPLY_BQSR module
    bqsr_result = GATK_APPLY_BQSR(
        readgroups_result.bam_with_read_groups.join(recal_result.recal_table),
        reference_file
    )
    
    // pass the output of GATK_BASE_RECALIBRATOR and the reference genome and index file to the PICARD_ADD_OR_REPLACE_READ_GROUPS module
    picard_result = PICARD_ADD_OR_REPLACE_READ_GROUPS(bqsr_result.recalibrated_bam)
    
    // create index for final BAM file
    index_result = SAMTOOLS_INDEX(picard_result.final_bam)
    
    // combine final BAM file and index file
    final_bams_with_index = picard_result.final_bam
        .map { sample_id, bam -> tuple(sample_id, bam) }
        .join(index_result.indexed_bam
            .map { sample_id, bam, bai -> tuple(sample_id, bai) }
        )
        .map { sample_id, bam, bai -> tuple(sample_id, bam, bai) }
    
    // split final BAM file into normal and tumor channels
    final_bams_normal = final_bams_with_index
        .filter { it[0].toString().endsWith('_dna_normal') }
    
    final_bams_tumor = final_bams_with_index
        .filter { it[0].toString().endsWith('_dna_tumor') }
    
    // run Mutect2 for normal samples
    normal_vcf = GATK_MUTECT2_NORMAL(
        final_bams_normal,
        reference_file
    )
    
    // prepare input data for paired analysis
    paired_input = final_bams_tumor
        .map { tumor_id, tumor_bam, tumor_bai -> 
            def base_id = tumor_id.replaceAll('_dna_tumor$', '')
            tuple(base_id, tumor_id, tumor_bam, tumor_bai)
        }
        .join(
            final_bams_normal.map { normal_id, normal_bam, normal_bai ->
                def base_id = normal_id.replaceAll('_dna_normal$', '')
                tuple(base_id, normal_id, normal_bam, normal_bai)
            }
        )
        .join(
            normal_vcf.normal_vcf.map { normal_id, vcf, idx ->
                def base_id = normal_id.replaceAll('_dna_normal$', '')
                tuple(base_id, vcf, idx)
            }
        )
        .map { base_id, 
               tumor_id, tumor_bam, tumor_bai,
               normal_id, normal_bam, normal_bai,
               vcf, idx ->
            tuple(tumor_id, tumor_bam, tumor_bai,
                  normal_id, normal_bam, normal_bai,
                  vcf, idx)
        }
    
    // run paired sample analysis
    paired_vcf = GATK_MUTECT2_PAIRED(
        paired_input,
        reference_file
    )
    
    // filter mutation results
    filtered_vcf = GATK_FILTER_MUTECT_CALLS(
        paired_vcf.paired_vcf,
        reference_file
    )
    
    // further filter mutation results using bcftools
    pass_filtered_vcf = BCFTOOLS_FILTER_PASS(
        filtered_vcf.filtered_vcf
    )
    
    // annotate mutation results using ANOVAR
    anovar_results = ANNOVAR(
        pass_filtered_vcf
    )
    
    
    // output results
    emit:
    annotation_results = anovar_results.anovar_txt
    annotation_vcf = anovar_results.anovar_vcf
}