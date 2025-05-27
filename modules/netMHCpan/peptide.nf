#!/usr/bin/env nextflow

// Workflow to combine DNA and RNA peptide sequences for mass spectrometry analysis

// Import modules from other workflows
include { DNA_DOWNSTREAM } from '../../workflows/dna_downstream'
include { RNA_DOWNSTREAM } from '../../workflows/rna_downstream'

// Process to combine and split peptides by length
process SPLIT_PEPTIDES {
    tag "${sample_id}"
    label 'process_medium'
    
    publishDir(
        path: { 
            def base_sample_id = sample_id.replaceAll('_(dna|rna)_(normal|tumor)$', '')
            return "${params.outdir}/${base_sample_id}/peptides_by_length" 
        },
        mode: 'copy'
    )
    
    input:
    tuple val(sample_id), path(dna_merged_fasta)
    tuple val(sample_id), path(rna_short_peptides)
    tuple val(sample_id), path(rna_long_peptides)
    
    output:
    tuple val(sample_id), path("${sample_id}_peptides_8.fasta"), emit: peptides_8
    tuple val(sample_id), path("${sample_id}_peptides_9.fasta"), emit: peptides_9
    tuple val(sample_id), path("${sample_id}_peptides_10.fasta"), emit: peptides_10
    tuple val(sample_id), path("${sample_id}_peptides_11.fasta"), emit: peptides_11
    tuple val(sample_id), path("${sample_id}_peptides_12.fasta"), emit: peptides_12
    tuple val(sample_id), path("${sample_id}_peptides_13.fasta"), emit: peptides_13
    tuple val(sample_id), path("${sample_id}_peptides_14.fasta"), emit: peptides_14
    tuple val(sample_id), path("${sample_id}_peptides_15.fasta"), emit: peptides_15
    tuple val(sample_id), path("${sample_id}_peptides_16.fasta"), emit: peptides_16
    tuple val(sample_id), path("${sample_id}_peptides_17.fasta"), emit: peptides_17
    tuple val(sample_id), path("${sample_id}_peptides_18.fasta"), emit: peptides_18
    tuple val(sample_id), path("${sample_id}_peptides_19.fasta"), emit: peptides_19
    tuple val(sample_id), path("${sample_id}_peptides_20.fasta"), emit: peptides_20
    tuple val(sample_id), path("${sample_id}_peptides_21.fasta"), emit: peptides_21
    tuple val(sample_id), path("${sample_id}_peptides_22.fasta"), emit: peptides_22
    tuple val(sample_id), path("${sample_id}_peptides_23.fasta"), emit: peptides_23
    tuple val(sample_id), path("${sample_id}_peptides_24.fasta"), emit: peptides_24
    tuple val(sample_id), path("${sample_id}_peptides_25.fasta"), emit: peptides_25
    path "${sample_id}_length_stats.txt", emit: stats
    
    script:
    def rna_short_files = rna_short_peptides instanceof List ? rna_short_peptides : [rna_short_peptides]
    def rna_long_files = rna_long_peptides instanceof List ? rna_long_peptides : [rna_long_peptides]
    """
    # Merge all FASTA files
    cat $dna_merged_fasta > combined_peptides.fasta
    
    # Add RNA short peptide files
    for file in ${rna_short_files.join(' ')}; do
        cat \$file >> combined_peptides.fasta
    done
    
    # Add RNA long peptide files
    for file in ${rna_long_files.join(' ')}; do
        cat \$file >> combined_peptides.fasta
    done

    # Create statistics file header
    echo "# Peptide length distribution for sample: ${sample_id}" > ${sample_id}_length_stats.txt
    echo "# Generated on \$(date)" >> ${sample_id}_length_stats.txt
    echo "Length\tCount" >> ${sample_id}_length_stats.txt

    # Classify peptides by length and generate statistics
    awk '
    BEGIN {
        for (i=8; i<=25; i++) count[i] = 0
    }
    /^>/ {
        header = \$0
        next
    }
    {
        len = length(\$0)
        if (len >= 8 && len <= 25) {
            print header > "'${sample_id}'_peptides_" len ".fasta"
            print \$0 > "'${sample_id}'_peptides_" len ".fasta"
            count[len]++
        }
    }
    END {
        for (i=8; i<=25; i++) {
            print i "\t" count[i] > "'${sample_id}'_length_stats.txt"
        }
    }' combined_peptides.fasta

    # Handle empty file cases
    for len in {8..25}; do
        if [ ! -f "${sample_id}_peptides_\${len}.fasta" ]; then
            echo "# No peptides of length \${len} found for sample: ${sample_id}" > "${sample_id}_peptides_\${len}.fasta"
        fi
    done

    # Clean up temporary files
    rm -f combined_peptides.fasta
    """
}

workflow PEPTIDE_ANALYSIS {
    take:
    dna_merged_fasta     // Channel: [val(sample_id), path(dna_merged_fasta)]
    rna_short_peptides   // Channel: [val(sample_id), path(rna_short_peptides)]
    rna_long_peptides    // Channel: [val(sample_id), path(rna_long_peptides)]
    
    main:
    peptide_results = SPLIT_PEPTIDES(
        dna_merged_fasta,
        rna_short_peptides,
        rna_long_peptides
    )
    
    emit:
    peptides_8 = peptide_results.peptides_8
    peptides_9 = peptide_results.peptides_9
    peptides_10 = peptide_results.peptides_10
    peptides_11 = peptide_results.peptides_11
    peptides_12 = peptide_results.peptides_12
    peptides_13 = peptide_results.peptides_13
    peptides_14 = peptide_results.peptides_14
    peptides_15 = peptide_results.peptides_15
    peptides_16 = peptide_results.peptides_16
    peptides_17 = peptide_results.peptides_17
    peptides_18 = peptide_results.peptides_18
    peptides_19 = peptide_results.peptides_19
    peptides_20 = peptide_results.peptides_20
    peptides_21 = peptide_results.peptides_21
    peptides_22 = peptide_results.peptides_22
    peptides_23 = peptide_results.peptides_23
    peptides_24 = peptide_results.peptides_24
    peptides_25 = peptide_results.peptides_25
    stats = peptide_results.stats
}

workflow {
    ch_dna_merged = Channel.fromPath(params.dna.merged_fasta, checkIfExists: true)
        .map { file -> 
            def sample_id = file.baseName.replaceAll('_(dna|rna)_(normal|tumor)$', '')
            tuple(sample_id, file) 
        }
    
    ch_rna_short = Channel.fromPath(params.rna.short_peptides, checkIfExists: true)
        .map { file -> 
            def sample_id = file.baseName.replaceAll('_(dna|rna)_(normal|tumor)$', '')
            tuple(sample_id, file) 
        }
    
    ch_rna_long = Channel.fromPath(params.rna.long_peptides, checkIfExists: true)
        .map { file -> 
            def sample_id = file.baseName.replaceAll('_(dna|rna)_(normal|tumor)$', '')
            tuple(sample_id, file) 
        }
    
    PEPTIDE_ANALYSIS(
        ch_dna_merged,
        ch_rna_short,
        ch_rna_long
    )
}

