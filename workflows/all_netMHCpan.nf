#!/usr/bin/env nextflow

// Import required modules
include { HLA_TYPING } from '../modules/netMHCpan/precess'
include { PEPTIDE_ANALYSIS } from '../modules/netMHCpan/peptide'
include { NETMHCPAN1_WORKFLOW } from '../modules/netMHCpan/netMHCpan1'
include { NETMHCPAN2_WORKFLOW } from '../modules/netMHCpan/netMHCpan2'
include { NETCTL_WORKFLOW } from '../modules/netMHCpan/netCTL'
include { DNA_DOWNSTREAM } from './dna_downstream'
include { RNA_DOWNSTREAM } from './rna_downstream'

// Main workflow for netMHCpan analysis
workflow NETMHCPAN_ANALYSIS {
    take:
    ch_dna_merged_fasta     // Channel: [val(sample_id), path(dna_merged_fasta)]
    ch_rna_short_peptides   // Channel: [val(sample_id), path(rna_short_peptides)]
    ch_rna_long_peptides    // Channel: [val(sample_id), path(rna_long_peptides)]
    ch_hla_i_results        // Channel: [val(sample_id), path(hla_i)]
    ch_hla_ii_results       // Channel: [val(sample_id), path(hla_ii)]
    
    main:
    // 1. format HLA typing results
    hla_data = HLA_TYPING(
        ch_hla_i_results,
        ch_hla_ii_results
    )

    // 2. process and split peptides
    peptide_data = PEPTIDE_ANALYSIS(
        ch_dna_merged_fasta,
        ch_rna_short_peptides,
        ch_rna_long_peptides
    )

    // 3. HLA-I class predictions (8-11mer)
    log.info "[DEBUG] Starting Class I predictions (8-11mer)"
    class1_predictions = NETMHCPAN1_WORKFLOW(
        peptide_data.peptides_8,
        peptide_data.peptides_9,
        peptide_data.peptides_10,
        peptide_data.peptides_11,
        hla_data.formatted_hla_i
    )

    // 4. HLA-II class predictions (12-25mer)
    log.info "[DEBUG] Starting Class II predictions (12-25mer)"
    def class2_predictions = NETMHCPAN2_WORKFLOW(
        peptide_data.peptides_12,
        peptide_data.peptides_13,
        peptide_data.peptides_14,
        peptide_data.peptides_15,
        peptide_data.peptides_16,
        peptide_data.peptides_17,
        peptide_data.peptides_18,
        peptide_data.peptides_19,
        peptide_data.peptides_20,
        peptide_data.peptides_21,
        peptide_data.peptides_22,
        peptide_data.peptides_23,
        peptide_data.peptides_24,
        peptide_data.peptides_25,
        hla_data.formatted_hla_ii
    )

    // 5. netCTL predictions (only 9mer)
    def netctl_input = peptide_data.peptides_9
        .combine(hla_data.formatted_hla_i, by: 0)
        .map { sample_id, pep, hla -> 
            log.info "[DEBUG] Creating netCTL input for sample: ${sample_id}"
            tuple(sample_id, pep, hla)
        }
    
    def netctl_predictions = NETCTL_WORKFLOW(netctl_input)

    emit:
    // HLA typing results
    hla_i_types = hla_data.formatted_hla_i
    hla_ii_types = hla_data.formatted_hla_ii
    
    // peptide data (emit each length separately)
    peptides_8 = peptide_data.peptides_8
    peptides_9 = peptide_data.peptides_9
    peptides_10 = peptide_data.peptides_10
    peptides_11 = peptide_data.peptides_11
    peptides_12 = peptide_data.peptides_12
    peptides_13 = peptide_data.peptides_13
    peptides_14 = peptide_data.peptides_14
    peptides_15 = peptide_data.peptides_15
    peptides_16 = peptide_data.peptides_16
    peptides_17 = peptide_data.peptides_17
    peptides_18 = peptide_data.peptides_18
    peptides_19 = peptide_data.peptides_19
    peptides_20 = peptide_data.peptides_20
    peptides_21 = peptide_data.peptides_21
    peptides_22 = peptide_data.peptides_22
    peptides_23 = peptide_data.peptides_23
    peptides_24 = peptide_data.peptides_24
    peptides_25 = peptide_data.peptides_25
    stats = peptide_data.stats
    
    // predictions - output each length
    class1_8mer_results = class1_predictions.predictions_8mer
    class1_9mer_results = class1_predictions.predictions_9mer
    class1_10mer_results = class1_predictions.predictions_10mer
    class1_11mer_results = class1_predictions.predictions_11mer
    
    // Class II predictions
    class2_12mer_results = class2_predictions.predictions_12mer
    class2_13mer_results = class2_predictions.predictions_13mer
    class2_14mer_results = class2_predictions.predictions_14mer
    class2_15mer_results = class2_predictions.predictions_15mer
    class2_16mer_results = class2_predictions.predictions_16mer
    class2_17mer_results = class2_predictions.predictions_17mer
    class2_18mer_results = class2_predictions.predictions_18mer
    class2_19mer_results = class2_predictions.predictions_19mer
    class2_20mer_results = class2_predictions.predictions_20mer
    class2_21mer_results = class2_predictions.predictions_21mer
    class2_22mer_results = class2_predictions.predictions_22mer
    class2_23mer_results = class2_predictions.predictions_23mer
    class2_24mer_results = class2_predictions.predictions_24mer
    class2_25mer_results = class2_predictions.predictions_25mer
    
    // netCTL predictions
    netctl_results = netctl_predictions.netctl_results
}

// main workflow entry
workflow {
    // prepare input channels
    ch_dna_merged_fasta = Channel.fromPath(params.dna.merged_fasta, checkIfExists: true)
        .map { file -> 
            def sample_id = file.baseName.replaceAll('_(dna|rna)_(normal|tumor)$', '')
            tuple(sample_id, file) 
        }
        .ifEmpty { error "No DNA merged FASTA found: ${params.dna.merged_fasta}" }
    
    ch_rna_short_peptides = Channel.fromPath(params.rna.short_peptides, checkIfExists: true)
        .map { file -> 
            def sample_id = file.baseName.replaceAll('_(dna|rna)_(normal|tumor)$', '')
            tuple(sample_id, file) 
        }
        .ifEmpty { error "No RNA short peptides found: ${params.rna.short_peptides}" }
    
    ch_rna_long_peptides = Channel.fromPath(params.rna.long_peptides, checkIfExists: true)
        .map { file -> 
            def sample_id = file.baseName.replaceAll('_(dna|rna)_(normal|tumor)$', '')
            tuple(sample_id, file) 
        }
        .ifEmpty { error "No RNA long peptides found: ${params.rna.long_peptides}" }
    
    ch_dna_anno = Channel.fromPath(params.dna.anno_results, checkIfExists: true)
        .map { file -> 
            def sample_id = file.baseName.replaceAll('_(dna|rna)_(normal|tumor)$', '')
            tuple(sample_id, file) 
        }
        .ifEmpty { error "No DNA annotation results found: ${params.dna.anno_results}" }
    
    ch_rna_anno = Channel.fromPath(params.rna.anno_results, checkIfExists: true)
        .map { file -> 
            def sample_id = file.baseName.replaceAll('_(dna|rna)_(normal|tumor)$', '')
            tuple(sample_id, file) 
        }
        .ifEmpty { error "No RNA annotation results found: ${params.rna.anno_results}" }
    
    ch_kallisto = Channel.fromPath(params.rna.kallisto_results, checkIfExists: true)
        .map { file -> 
            def sample_id = file.baseName.replaceAll('_(dna|rna)_(normal|tumor)$', '')
            tuple(sample_id, file) 
        }
        .ifEmpty { error "No Kallisto results found: ${params.rna.kallisto_results}" }
    
    ch_star_fusion = Channel.fromPath(params.rna.star_fusion_results, checkIfExists: true)
        .map { file -> 
            def sample_id = file.baseName.replaceAll('_(dna|rna)_(normal|tumor)$', '')
            tuple(sample_id, file) 
        }
        .ifEmpty { log.warn "No STAR-Fusion results found: ${params.rna.star_fusion_results}" }
    
    ch_hla_i = Channel.fromPath(params.rna.hla_i_results, checkIfExists: true)
        .map { file -> 
            def sample_id = file.baseName.replaceAll('_(dna|rna)_(normal|tumor)$', '')
            tuple(sample_id, file) 
        }
        .ifEmpty { error "No HLA-I typing results found: ${params.rna.hla_i_results}" }
    
    ch_hla_ii = Channel.fromPath(params.rna.hla_ii_results, checkIfExists: true)
        .map { file -> 
            def sample_id = file.baseName.replaceAll('_(dna|rna)_(normal|tumor)$', '')
            tuple(sample_id, file) 
        }
        .ifEmpty { error "No HLA-II typing results found: ${params.rna.hla_ii_results}" }
    
    // run netMHCpan analysis workflow
    NETMHCPAN_ANALYSIS(
        ch_dna_merged_fasta,
        ch_rna_short_peptides,
        ch_rna_long_peptides,
        ch_hla_i,
        ch_hla_ii
    )
}



