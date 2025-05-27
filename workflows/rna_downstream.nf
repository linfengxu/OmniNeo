#!/usr/bin/env nextflow

// downstream analysis workflow for RNA analysis, process annotated results and kallisto quantification results, perform subsequent expression filtering and peptide generation

// import modules
include { RNA_TPM_FILTER } from '../modules/downstream/rnatpmfilter'
include { RNA_NONCODING_ANALYSIS } from '../modules/downstream/rnanoncoding'
include { RNA_STAR_FUSION_PEPTIDES } from '../modules/downstream/rnastarfusion'

workflow RNA_DOWNSTREAM {
    take:
    annovar_results         // RNA mutation annotation results from ANNOVAR [sample_id, annovar_txt]
    kallisto_abundance      // Kallisto quantification results [sample_id, abundance_tsv]
    star_fusion_results     // STAR-Fusion results [sample_id, fusion_results]
    
    main:
    // step 1: filter RNA mutations based on gene expression values
    // need transcript-gene mapping file - use the path defined in the configuration file
    gene_ref_file = file(params.kallisto.gene_transcript_map)
    
    // apply TPM filtering
    tpm_filtered_results = RNA_TPM_FILTER(
        kallisto_abundance,
        gene_ref_file,
        annovar_results
    )
    
    // step 2: analyze noncoding mutations
    noncoding_results = RNA_NONCODING_ANALYSIS(
        tpm_filtered_results.filtered_variants
    )
    
    // normalize noncoding results
    noncoding_short_normalized = noncoding_results.noncoding_short_peptides
        .map { sample_id, file -> 
            def base_id = sample_id.replaceAll('_(dna|rna)_(normal|tumor)$', '')
            return [base_id, file] 
        }
     noncoding_long_normalized = noncoding_results.noncoding_long_peptides
        .map { sample_id, file -> 
            def base_id = sample_id.replaceAll('_(dna|rna)_(normal|tumor)$', '')
            return [base_id, file] 
        }
    // step 3: generate fusion peptides
    // check if star_fusion_results is empty
    star_fusion_results_exists = star_fusion_results
        .map { sample_id, file1,file2,file3 -> 
            def base_id = sample_id.replaceAll('_(dna|rna)_(normal|tumor)$', '')
            return [base_id, file3] 
        }
        .ifEmpty { log.warn("No STAR-Fusion results provided, skipping fusion peptide generation"); return Channel.empty() }
    
    // run peptide generation process only when fusion results exist
    fusion_results = RNA_STAR_FUSION_PEPTIDES(
        star_fusion_results_exists
    )
    
    // prepare empty channels to ensure output channels work even when no fusion results exist
    empty_fusion_short = Channel.empty()
    empty_fusion_long = Channel.empty()
    empty_fusion_short_csv = Channel.empty()
    empty_fusion_long_csv = Channel.empty()
    
    // merge actual results with empty channels
    fusion_short_peptides = fusion_results.fusion_short_peptides.ifEmpty { empty_fusion_short }
    fusion_long_peptides = fusion_results.fusion_long_peptides.ifEmpty { empty_fusion_long }
    fusion_short_csv = fusion_results.fusion_short_csv.ifEmpty { empty_fusion_short_csv }
    fusion_long_csv = fusion_results.fusion_long_csv.ifEmpty { empty_fusion_long_csv }
    
    // integrate all results
    
    emit:
    // TPM filtering results
    filtered_variants = tpm_filtered_results.filtered_variants
    gene_expression = tpm_filtered_results.gene_expression
    
    // fusion peptides
    fusion_short_peptides = fusion_short_peptides
    fusion_long_peptides = fusion_long_peptides
    
    // peptide channels for MS analysis
    all_short_peptides = noncoding_short_normalized
        .map { id, files -> [id, [source: 'noncoding', files: files]] }
        .mix(
            fusion_short_peptides.map { id, files -> [id, [source: 'fusion', files: files]] }
        )
        .groupTuple()
        .map { id, data_list ->
            def all_files = []
            data_list.each { item ->
                if (item.files instanceof List) {
                    all_files.addAll(item.files)
                } else {
                    all_files.add(item.files)
                }
            }
            return tuple(id, all_files)
        }
    
    all_long_peptides = noncoding_long_normalized
        .map { id, files -> [id, [source: 'noncoding', files: files]] }
        .mix(
            fusion_long_peptides.map { id, files -> [id, [source: 'fusion', files: files]] }
        )
        .groupTuple()
        .map { id, data_list ->
            def all_files = []
            data_list.each { item ->
                if (item.files instanceof List) {
                    all_files.addAll(item.files)
                } else {
                    all_files.add(item.files)
                }
            }
            return tuple(id, all_files)
        }
        .ifEmpty { log.warn "No RNA long peptides found"; return Channel.empty() }
}
