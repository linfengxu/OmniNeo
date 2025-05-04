#!/usr/bin/env nextflow

// Workflow to combine DNA and RNA peptide sequences for mass spectrometry analysis

// Import modules from other workflows
include { DNA_DOWNSTREAM } from './dna_downstream'
include { RNA_DOWNSTREAM } from './rna_downstream'

// Process to combine all FASTA files into one MS database
process COMBINE_MS_DATABASE {
    tag "${sample_id}"
    label 'process_medium'
    publishDir "${params.outdir}/${sample_id}/ms_database", mode: 'copy'
    
    input:
    tuple val(sample_id), path(dna_merged_fasta)
    tuple val(sample_id), path(rna_short_peptides)
    tuple val(sample_id), path(rna_long_peptides)
    path(uniprot_fasta)
    path(crap_fasta)
    
    output:
    tuple val(sample_id), path("${sample_id}_combined_ms_database.fasta"), emit: combined_ms_db
    path "${sample_id}_database_stats.txt", emit: database_stats
    
    script:
    """
    # Create empty files if needed
    touch empty.fa
    
    # Create headers for the combined database
    echo "# Combined MS database for sample: ${sample_id}" > ${sample_id}_combined_ms_database.fasta
    echo "# Created on \$(date)" >> ${sample_id}_combined_ms_database.fasta
    echo "# Contains DNA mutations, RNA variants, and reference databases" >> ${sample_id}_combined_ms_database.fasta
    echo "" >> ${sample_id}_combined_ms_database.fasta
    
    # Add DNA mutation peptides
    if [ -s "${dna_merged_fasta}" ]; then
        echo "# === DNA MUTATION PEPTIDES ===" >> ${sample_id}_combined_ms_database.fasta
        cat ${dna_merged_fasta} >> ${sample_id}_combined_ms_database.fasta
        echo "" >> ${sample_id}_combined_ms_database.fasta
    else
        echo "# === NO DNA MUTATION PEPTIDES AVAILABLE ===" >> ${sample_id}_combined_ms_database.fasta
        echo "" >> ${sample_id}_combined_ms_database.fasta
    fi
    
    # Add RNA short peptides
    if [ -s "${rna_short_peptides}" ]; then
        echo "# === RNA SHORT PEPTIDES ===" >> ${sample_id}_combined_ms_database.fasta
        cat ${rna_short_peptides} >> ${sample_id}_combined_ms_database.fasta
        echo "" >> ${sample_id}_combined_ms_database.fasta
    else
        echo "# === NO RNA SHORT PEPTIDES AVAILABLE ===" >> ${sample_id}_combined_ms_database.fasta
        echo "" >> ${sample_id}_combined_ms_database.fasta
    fi
    
    # Add RNA long peptides
    if [ -s "${rna_long_peptides}" ]; then
        echo "# === RNA LONG PEPTIDES ===" >> ${sample_id}_combined_ms_database.fasta
        cat ${rna_long_peptides} >> ${sample_id}_combined_ms_database.fasta
        echo "" >> ${sample_id}_combined_ms_database.fasta
    else
        echo "# === NO RNA LONG PEPTIDES AVAILABLE ===" >> ${sample_id}_combined_ms_database.fasta
        echo "" >> ${sample_id}_combined_ms_database.fasta
    fi
    
    # Add reference databases
    echo "# === REFERENCE PROTEOME (UniProt) ===" >> ${sample_id}_combined_ms_database.fasta
    cat ${uniprot_fasta} >> ${sample_id}_combined_ms_database.fasta
    echo "" >> ${sample_id}_combined_ms_database.fasta
    
    echo "# === CONTAMINANTS DATABASE (cRAP) ===" >> ${sample_id}_combined_ms_database.fasta
    cat ${crap_fasta} >> ${sample_id}_combined_ms_database.fasta
    
    # Generate statistics
    echo "Database Statistics for ${sample_id}" > ${sample_id}_database_stats.txt
    echo "-------------------------------------" >> ${sample_id}_database_stats.txt
    echo "Total entries: \$(grep -c "^>" ${sample_id}_combined_ms_database.fasta)" >> ${sample_id}_database_stats.txt
    echo "DNA mutation peptides: \$(grep -c "^>" ${dna_merged_fasta} 2>/dev/null || echo 0)" >> ${sample_id}_database_stats.txt
    echo "RNA short peptides: \$(grep -c "^>" ${rna_short_peptides} 2>/dev/null || echo 0)" >> ${sample_id}_database_stats.txt
    echo "RNA long peptides: \$(grep -c "^>" ${rna_long_peptides} 2>/dev/null || echo 0)" >> ${sample_id}_database_stats.txt
    echo "UniProt reference: \$(grep -c "^>" ${uniprot_fasta})" >> ${sample_id}_database_stats.txt 
    echo "Contaminants: \$(grep -c "^>" ${crap_fasta})" >> ${sample_id}_database_stats.txt
    echo "File size: \$(du -h ${sample_id}_combined_ms_database.fasta | cut -f1)" >> ${sample_id}_database_stats.txt
    echo "Created: \$(date)" >> ${sample_id}_database_stats.txt
    
    # Log completion
    echo "Combined MS database created successfully for sample ${sample_id}"
    """
}

// Main workflow
workflow COMBINE_MS {
    take:
    dna_merged_fasta        // Channel: [ val(sample_id), path(merged_fasta) ]
    rna_short_peptides      // Channel: [ val(sample_id), path(short_peptides) ]
    rna_long_peptides       // Channel: [ val(sample_id), path(long_peptides) ]
    uniprot_fasta           // Path: reference proteome from UniProt
    crap_fasta              // Path: contaminants database
    
    main:
    // Prepare reference databases
    ch_uniprot_fasta = file(uniprot_fasta)
    ch_crap_fasta = file(crap_fasta)
    
    // Combine all input into one channel for processing - transform to match expected input types
    // First join the DNA and RNA channels by sample_id
    combined_input = dna_merged_fasta
        .join(rna_short_peptides, failOnMismatch: false, remainder: true)
        .join(rna_long_peptides, failOnMismatch: false, remainder: true)
        .map { items ->
            def sample_id = items[0]
            def dna_fasta = items.size() > 1 ? items[1] : file('empty.fa')
            def rna_short = items.size() > 2 ? items[2] : file('empty.fa')
            def rna_long = items.size() > 3 ? items[3] : file('empty.fa')
            
            return [
                tuple(sample_id, dna_fasta),
                tuple(sample_id, rna_short),
                tuple(sample_id, rna_long)
            ]
        }
        .transpose()  // Get each tuple separately
        .branch {
            dna: it[1].name.contains('merged') || it[1].name == 'empty.fa'
            rna_short: it[1].name.contains('short') || it[1].name == 'empty.fa'
            rna_long: it[1].name.contains('long') || it[1].name == 'empty.fa'
        }
    
    // Build the MS database with correctly structured inputs
    combined_ms_database = COMBINE_MS_DATABASE(
        combined_input.dna,
        combined_input.rna_short,
        combined_input.rna_long,
        ch_uniprot_fasta,
        ch_crap_fasta
    )
    
    emit:
    combined_ms_db = combined_ms_database.combined_ms_db
    database_stats = combined_ms_database.database_stats
}

// This allows the workflow to be run directly when it's included
workflow {
    // Get parameters for the reference databases
    uniprot_db = file(params.combine_ms.uniprot_fasta)
    crap_db = file(params.combine_ms.crap_fasta)
    
    // Get sample data from input channels
    ch_dna_anno_results = Channel.fromPath(params.dna.anno_results, checkIfExists: true)
        .map { file -> tuple(file.baseName, file) }
    
    ch_rna_anno_results = Channel.fromPath(params.rna.anno_results, checkIfExists: true)
        .map { file -> tuple(file.baseName, file) }
    
    ch_kallisto = Channel.fromPath(params.rna.kallisto_results, checkIfExists: true)
        .map { file -> tuple(file.baseName, file) }
    
    ch_star_fusion = Channel.fromPath(params.rna.star_fusion_results, checkIfExists: true)
        .map { file -> tuple(file.baseName, file) }
    
    // Empty channels for optional RNA fusion peptides
    ch_fusion_long = Channel.empty()
    ch_fusion_short = Channel.empty()
    
    // Run DNA downstream analysis
    DNA_DOWNSTREAM(ch_dna_anno_results)
    
    // Run RNA downstream analysis
    RNA_DOWNSTREAM(
        ch_rna_anno_results,
        ch_kallisto,
        ch_star_fusion,
        ch_fusion_long,
        ch_fusion_short
    )
    
    // Combine all results for MS analysis
    COMBINE_MS(
        DNA_DOWNSTREAM.out.merged_fasta,
        RNA_DOWNSTREAM.out.all_short_peptides,
        RNA_DOWNSTREAM.out.all_long_peptides,
        uniprot_db,
        crap_db
    )
}
