#!/usr/bin/env nextflow

// Import all required modules
include { FRAMESHIFT_DELETION_ANALYSIS     } from '../modules/downstream/frameshift_deletion'
include { FRAMESHIFT_INSERTION_ANALYSIS    } from '../modules/downstream/frameshift_insertion'
include { NONFRAMESHIFT_DELETION_ANALYSIS  } from '../modules/downstream/nonframeshift_deletion'
include { NONFRAMESHIFT_INSERTION_ANALYSIS } from '../modules/downstream/nonframeshift_insertion'
include { NONFRAMESHIFT_SUBSTITUTION_ANALYSIS } from '../modules/downstream/nonframeshift_substitution'
include { SNV_ANALYSIS                     } from '../modules/downstream/SNV'
include { STOPLOSS_MUTATION_ANALYSIS       } from '../modules/downstream/stoploss_mutation'

// Define a new process to merge all fasta files
process MERGE_FASTA_FILES {
    tag "$sample_id"
    label 'process_low'
    container "${params.snv_analysis.container}"
    
    publishDir(
        path: { 
            def base_sample_id = sample_id.replaceAll('_(dna|rna)_(normal|tumor)$', '')
            return "${params.outdir}/${base_sample_id}/dna/10_peptides/merged" 
        },
        mode: 'copy'
    )
    
    input:
    tuple val(sample_id), path('input_fastas/*')
    
    output:
    tuple val(sample_id), path("merged_peptides.fasta"), emit: merged_fasta
    
    script:
    """
    # Create header with sample info and date
    echo "# Merged peptide sequences for sample: ${sample_id}" > merged_peptides.fasta
    echo "# Created on \$(date)" >> merged_peptides.fasta
    echo "# Contains mutation peptides from all analysis types" >> merged_peptides.fasta
    echo "" >> merged_peptides.fasta
    
    # Concatenate all fasta files
    cat input_fastas/* >> merged_peptides.fasta
    
    # Report the number of sequences merged
    echo "Merged \$(grep -c "^>" merged_peptides.fasta) peptide sequences into one file"
    """
}

workflow DNA_DOWNSTREAM {
    take:
    ch_dna_anno_results  // Channel: [ val(sample_id), path(anno_result) ]

    main:
    // Perform all types of mutation analyses in parallel
    FRAMESHIFT_DELETION_ANALYSIS(ch_dna_anno_results)
    FRAMESHIFT_INSERTION_ANALYSIS(ch_dna_anno_results)
    NONFRAMESHIFT_DELETION_ANALYSIS(ch_dna_anno_results)
    NONFRAMESHIFT_INSERTION_ANALYSIS(ch_dna_anno_results)
    NONFRAMESHIFT_SUBSTITUTION_ANALYSIS(ch_dna_anno_results)
    SNV_ANALYSIS(ch_dna_anno_results)
    STOPLOSS_MUTATION_ANALYSIS(ch_dna_anno_results)

    // Collect all fasta files from fasta_files directory
    ch_peptides_fasta = FRAMESHIFT_DELETION_ANALYSIS.out.fsdel_fasta_files
        .mix(FRAMESHIFT_INSERTION_ANALYSIS.out.fsins_fasta_files)
        .mix(NONFRAMESHIFT_DELETION_ANALYSIS.out.nfsdel_fasta_files)
        .mix(NONFRAMESHIFT_INSERTION_ANALYSIS.out.nfsins_fasta_files)
        .mix(NONFRAMESHIFT_SUBSTITUTION_ANALYSIS.out.nfssub_fasta_files)
        .mix(SNV_ANALYSIS.out.snv_fasta_files)
        .mix(STOPLOSS_MUTATION_ANALYSIS.out.stoploss_fasta_files)

    // Group all fasta files by sample_id for merging
    ch_fasta_by_sample = ch_peptides_fasta.groupTuple()
    
    // Merge all fasta files into one large file
    MERGE_FASTA_FILES(ch_fasta_by_sample)

    ch_gene_summaries = FRAMESHIFT_DELETION_ANALYSIS.out.fsdel_gene_summary
        .mix(FRAMESHIFT_INSERTION_ANALYSIS.out.fsins_gene_summary)
        .mix(NONFRAMESHIFT_DELETION_ANALYSIS.out.nfsdel_gene_summary)
        .mix(NONFRAMESHIFT_INSERTION_ANALYSIS.out.nfsins_gene_summary)
        .mix(NONFRAMESHIFT_SUBSTITUTION_ANALYSIS.out.nfssub_gene_summary)
        .mix(SNV_ANALYSIS.out.snv_gene_summary)
        .mix(STOPLOSS_MUTATION_ANALYSIS.out.stoploss_gene_summary)

    emit:
    peptides_fasta = ch_peptides_fasta
    merged_fasta = MERGE_FASTA_FILES.out.merged_fasta
    gene_summaries = ch_gene_summaries
    
    // 单独emit每个分析类型的输出通道，而不是整个输出对象
    // Frameshift Deletion outputs
    fsdel_fasta_files = FRAMESHIFT_DELETION_ANALYSIS.out.fsdel_fasta_files
    fsdel_csv_files = FRAMESHIFT_DELETION_ANALYSIS.out.fsdel_csv_files
    fsdel_txt_files = FRAMESHIFT_DELETION_ANALYSIS.out.fsdel_txt_files
    fsdel_gene_summary = FRAMESHIFT_DELETION_ANALYSIS.out.fsdel_gene_summary
    fsdel_var_sequence = FRAMESHIFT_DELETION_ANALYSIS.out.fsdel_var_sequence
    
    // Frameshift Insertion outputs
    fsins_fasta_files = FRAMESHIFT_INSERTION_ANALYSIS.out.fsins_fasta_files
    fsins_csv_files = FRAMESHIFT_INSERTION_ANALYSIS.out.fsins_csv_files
    fsins_txt_files = FRAMESHIFT_INSERTION_ANALYSIS.out.fsins_txt_files
    fsins_gene_summary = FRAMESHIFT_INSERTION_ANALYSIS.out.fsins_gene_summary
    fsins_var_sequence = FRAMESHIFT_INSERTION_ANALYSIS.out.fsins_var_sequence
    
    // Nonframeshift Deletion outputs
    nfsdel_fasta_files = NONFRAMESHIFT_DELETION_ANALYSIS.out.nfsdel_fasta_files
    nfsdel_csv_files = NONFRAMESHIFT_DELETION_ANALYSIS.out.nfsdel_csv_files
    nfsdel_txt_files = NONFRAMESHIFT_DELETION_ANALYSIS.out.nfsdel_txt_files
    nfsdel_gene_summary = NONFRAMESHIFT_DELETION_ANALYSIS.out.nfsdel_gene_summary
    nfsdel_var_sequence = NONFRAMESHIFT_DELETION_ANALYSIS.out.nfsdel_var_sequence
    
    // Nonframeshift Insertion outputs
    nfsins_fasta_files = NONFRAMESHIFT_INSERTION_ANALYSIS.out.nfsins_fasta_files
    nfsins_csv_files = NONFRAMESHIFT_INSERTION_ANALYSIS.out.nfsins_csv_files
    nfsins_txt_files = NONFRAMESHIFT_INSERTION_ANALYSIS.out.nfsins_txt_files
    nfsins_gene_summary = NONFRAMESHIFT_INSERTION_ANALYSIS.out.nfsins_gene_summary
    nfsins_var_sequence = NONFRAMESHIFT_INSERTION_ANALYSIS.out.nfsins_var_sequence
    
    // Nonframeshift Substitution outputs
    nfssub_fasta_files = NONFRAMESHIFT_SUBSTITUTION_ANALYSIS.out.nfssub_fasta_files
    nfssub_csv_files = NONFRAMESHIFT_SUBSTITUTION_ANALYSIS.out.nfssub_csv_files
    nfssub_txt_files = NONFRAMESHIFT_SUBSTITUTION_ANALYSIS.out.nfssub_txt_files
    nfssub_gene_summary = NONFRAMESHIFT_SUBSTITUTION_ANALYSIS.out.nfssub_gene_summary
    nfssub_var_sequence = NONFRAMESHIFT_SUBSTITUTION_ANALYSIS.out.nfssub_var_sequence
    
    // SNV outputs
    snv_fasta_files = SNV_ANALYSIS.out.snv_fasta_files
    snv_csv_files = SNV_ANALYSIS.out.snv_csv_files
    snv_txt_files = SNV_ANALYSIS.out.snv_txt_files
    snv_gene_summary = SNV_ANALYSIS.out.snv_gene_summary
    snv_var_sequence = SNV_ANALYSIS.out.snv_var_sequence
    
    // Stoploss outputs
    stoploss_fasta_files = STOPLOSS_MUTATION_ANALYSIS.out.stoploss_fasta_files
    stoploss_csv_files = STOPLOSS_MUTATION_ANALYSIS.out.stoploss_csv_files
    stoploss_txt_files = STOPLOSS_MUTATION_ANALYSIS.out.stoploss_txt_files
    stoploss_gene_summary = STOPLOSS_MUTATION_ANALYSIS.out.stoploss_gene_summary
    stoploss_var_sequence = STOPLOSS_MUTATION_ANALYSIS.out.stoploss_var_sequence
}

// This allows the workflow to be run directly when it's included
workflow {
    // Define the input channel - for standalone execution
    ch_dna_anno_results = Channel.fromPath(params.anno_results, checkIfExists: true)
        .map { file -> tuple(file.baseName, file) }
    
    DNA_DOWNSTREAM(ch_dna_anno_results)
}