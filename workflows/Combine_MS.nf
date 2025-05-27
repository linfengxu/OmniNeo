#!/usr/bin/env nextflow

// Workflow to combine DNA and RNA peptide sequences for mass spectrometry analysis

// Import modules from other workflows
include { DNA_DOWNSTREAM } from './dna_downstream'
include { RNA_DOWNSTREAM } from './rna_downstream'

// Process to combine all FASTA files into one MS database
process COMBINE_MS_DATABASE {
    tag "${sample_id}"
    label 'process_medium'
    publishDir(
        path: { 
            def base_sample_id = sample_id.replaceAll('_(dna|rna)_(normal|tumor)$', '')
            return "${params.outdir}/${base_sample_id}/ms_database" 
        },
        mode: 'copy'
    )
    
    input:
    tuple val(sample_id), path(dna_merged_fasta)
    tuple val(sample_id), path(rna_short_peptides)
    tuple val(sample_id), path(rna_long_peptides)
    path(uniprot_fasta)
    path(crap_fasta)
    
    output:
    tuple val(sample_id), path("${sample_id}_combined_ms_database.fasta"), emit: combined_ms_db
    tuple val(sample_id), path("${sample_id}_clean_ms_database.fasta"), emit: clean_ms_db
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
    
    # Add RNA peptides
    echo "# === RNA PEPTIDES ===" >> ${sample_id}_combined_ms_database.fasta
    
    # Add RNA short peptides (both fusion and noncoding)
    echo "# --- RNA SHORT PEPTIDES ---" >> ${sample_id}_combined_ms_database.fasta
    if [ -n "\$(ls -A ${rna_short_peptides} 2>/dev/null)" ]; then
        for file in ${rna_short_peptides}; do
            if [ -s "\$file" ]; then
                echo "# Adding file: \$(basename "\$file")" >> ${sample_id}_combined_ms_database.fasta
                cat "\$file" >> ${sample_id}_combined_ms_database.fasta
                echo "" >> ${sample_id}_combined_ms_database.fasta
            fi
        done
    else
        echo "# No RNA short peptides found" >> ${sample_id}_combined_ms_database.fasta
    fi
    
    # Add RNA long peptides (both fusion and noncoding)
    echo "# --- RNA LONG PEPTIDES ---" >> ${sample_id}_combined_ms_database.fasta
    if [ -n "\$(ls -A ${rna_long_peptides} 2>/dev/null)" ]; then
        for file in ${rna_long_peptides}; do
            if [ -s "\$file" ]; then
                echo "# Adding file: \$(basename "\$file")" >> ${sample_id}_combined_ms_database.fasta
                cat "\$file" >> ${sample_id}_combined_ms_database.fasta
                echo "" >> ${sample_id}_combined_ms_database.fasta
            fi
        done
    else
        echo "# No RNA long peptides found" >> ${sample_id}_combined_ms_database.fasta
    fi
    
    # Add reference databases
    echo "# === REFERENCE PROTEOME (UniProt) ===" >> ${sample_id}_combined_ms_database.fasta
    cat ${uniprot_fasta} >> ${sample_id}_combined_ms_database.fasta
    echo "" >> ${sample_id}_combined_ms_database.fasta
    
    echo "# === CONTAMINANTS DATABASE (cRAP) ===" >> ${sample_id}_combined_ms_database.fasta
    cat ${crap_fasta} >> ${sample_id}_combined_ms_database.fasta
    
    # Generate clean FASTA without comments
    grep -v "^\\#" ${sample_id}_combined_ms_database.fasta | grep -v "^\$" > ${sample_id}_temp_clean_ms_database.fasta
    
    # Use reformat_fasta.py script to reformat the FASTA headers for MaxQuant compatibility
    python ${projectDir}/bin/reformat_fasta.py ${sample_id}_temp_clean_ms_database.fasta ${sample_id}_clean_ms_database.fasta
    rm ${sample_id}_temp_clean_ms_database.fasta
    
    # Generate statistics
    echo "Database Statistics for ${sample_id}" > ${sample_id}_database_stats.txt
    echo "-------------------------------------" >> ${sample_id}_database_stats.txt
    echo "Total entries: \$(grep -c "^>" ${sample_id}_combined_ms_database.fasta)" >> ${sample_id}_database_stats.txt
    echo "DNA mutation peptides: \$(grep -c "^>" ${dna_merged_fasta} 2>/dev/null || echo 0)" >> ${sample_id}_database_stats.txt
    
    # Count RNA peptides
    echo "RNA peptides:" >> ${sample_id}_database_stats.txt
    echo "--- Short peptides:" >> ${sample_id}_database_stats.txt
    short_count=0
    for file in ${rna_short_peptides}; do
        if [ -s "\$file" ]; then
            count=\$(grep -c "^>" "\$file" 2>/dev/null || echo 0)
            short_count=\$((short_count + count))
            echo "  \$(basename "\$file"): \$count" >> ${sample_id}_database_stats.txt
        fi
    done
    echo "  Total short peptides: \$short_count" >> ${sample_id}_database_stats.txt
    
    echo "--- Long peptides:" >> ${sample_id}_database_stats.txt
    long_count=0
    for file in ${rna_long_peptides}; do
        if [ -s "\$file" ]; then
            count=\$(grep -c "^>" "\$file" 2>/dev/null || echo 0)
            long_count=\$((long_count + count))
            echo "  \$(basename "\$file"): \$count" >> ${sample_id}_database_stats.txt
        fi
    done
    echo "  Total long peptides: \$long_count" >> ${sample_id}_database_stats.txt
    echo "  Total RNA peptides: \$((short_count + long_count))" >> ${sample_id}_database_stats.txt
    
    echo "UniProt reference: \$(grep -c "^>" ${uniprot_fasta})" >> ${sample_id}_database_stats.txt 
    echo "Contaminants: \$(grep -c "^>" ${crap_fasta})" >> ${sample_id}_database_stats.txt
    echo "File size: \$(du -h ${sample_id}_combined_ms_database.fasta | cut -f1)" >> ${sample_id}_database_stats.txt
    echo "Clean file size: \$(du -h ${sample_id}_clean_ms_database.fasta | cut -f1)" >> ${sample_id}_database_stats.txt
    echo "Created: \$(date)" >> ${sample_id}_database_stats.txt
    
    # Log completion
    echo "Combined MS database created successfully for sample ${sample_id}"
    """
}

// Main workflow
workflow COMBINE_MS {
    take:
    dna_merged_fasta 
    rna_short_peptides      // Channel: [ val(sample_id), path(short_peptides) ]
    rna_long_peptides       // Channel: [ val(sample_id), path(long_peptides) ]
    uniprot_fasta           // Path: reference proteome from UniProt
    crap_fasta              // Path: contaminants database
    
    main:
    // Add debug logging
    log.info "Starting COMBINE_MS workflow"
    log.info "DNA merged fasta: ${dna_merged_fasta}"
    log.info "RNA short peptides: ${rna_short_peptides}"
    log.info "RNA long peptides: ${rna_long_peptides}"
    
    // Prepare reference databases
    ch_uniprot_fasta = file(uniprot_fasta)
    ch_crap_fasta = file(crap_fasta)
    
    // Build the MS database directly with the input channels
    combined_ms_database = COMBINE_MS_DATABASE(
        dna_merged_fasta,
        rna_short_peptides,
        rna_long_peptides,
        ch_uniprot_fasta,
        ch_crap_fasta
    )
    
    emit:
    combined_ms_db = combined_ms_database.combined_ms_db
    clean_ms_db = combined_ms_database.clean_ms_db
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
        .map { sample_id, file1, file2, file3 -> 
            def base_id = sample_id.replaceAll('_(dna|rna)_(normal|tumor)$', '')
            return [base_id, file1, file2, file3] 
        }
    
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
