#!/usr/bin/env nextflow

process NETCTL_PREDICT {
    tag "$sample_id"
    label 'process_medium'
    
    publishDir(
        path: { 
            def base_sample_id = sample_id.replaceAll('_(dna|rna)_(normal|tumor)$', '')
            return "${params.outdir}/${base_sample_id}/epitope_prediction/netCTL" 
        },
        mode: 'copy'
    )
    
    input:
    tuple val(sample_id), path(peptide_file), path(hla_file)  // Combined input with peptide and HLA files
    
    output:
    tuple val(sample_id), path("${sample_id}_netCTL_predictions.txt"), emit: netctl_results
    
    script:
    """
    # Read HLA alleles and convert to supertypes
    awk '
    BEGIN {
        # Define supertype mapping
        split("A*01:01 A1 A*02:01 A2 A*03:01 A3 A*24:02 A24 B*07:02 B7 B*08:01 B8 B*35:01 B35 B*44:02 B44 B*44:03 B44 C*01:02 C1 C*03:04 C3 C*04:01 C4 C*05:01 C5 C*06:02 C6 C*07:01 C7 C*08:01 C8", map)
        for (i=1; i<=length(map); i+=2) {
            supertypes[map[i]] = map[i+1]
        }
    }
    {
        # Read HLA alleles
        if (\$0 ~ /^HLA-/) {
            # Remove HLA- prefix
            allele = substr(\$0, 5)
            # Get supertype
            if (supertypes[allele]) {
                if (!seen[supertypes[allele]]) {
                    print supertypes[allele]
                    seen[supertypes[allele]] = 1
                }
            }
        }
    }' ${hla_file} > supertypes.txt

    # Extract 9mer peptides
    awk '
    BEGIN { seq = ""; header = "" }
    /^>/ {
        if (seq != "" && length(seq) == 9) {
            print header
            print seq
        }
        header = \$0
        seq = ""
    }
    !/^>/ {
        seq = seq \$0
    }
    END {
        if (seq != "" && length(seq) == 9) {
            print header
            print seq
        }
    }' ${peptide_file} > 9mer_sequences.fasta

    # Run prediction for each supertype
    while read supertype; do
        # Run netCTL prediction
        netCTL \\
            -f 9mer_sequences.fasta \\
            -s \$supertype \\
            -wt 0.05 \\
            -we 1.0 \\
            -wc 0.15 \\
            -thr 0.75 \\
            -sort 0 >> ${sample_id}_netCTL_predictions.txt
            
        # Add separator
        echo "----------------------------------------" >> ${sample_id}_netCTL_predictions.txt
    done < supertypes.txt

    # Add prediction information
    echo "Prediction completed at: \$(date)" >> ${sample_id}_netCTL_predictions.txt
    echo "netCTL version: \$(netCTL -v 2>&1)" >> ${sample_id}_netCTL_predictions.txt

    # Clean up temporary files
    rm supertypes.txt 9mer_sequences.fasta
    """
}

workflow NETCTL_WORKFLOW {
    take:
    peptides_and_hla  // Combined channel with peptides and HLA files
    
    main:
    // Run netCTL prediction
    netctl_results = NETCTL_PREDICT(peptides_and_hla)
    
    emit:
    netctl_results = netctl_results.netctl_results
}
