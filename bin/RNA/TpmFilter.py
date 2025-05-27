import pandas as pd
import os
import argparse
import sys

def tpm_filter_with_annovar(kallisto_file, ref_file, annovar_file, output_path, tpm_threshold=0):
    """
    Filter genes based on TPM values from kallisto and compare with Annovar results.
    
    Parameters:
    -----------
    kallisto_file: str
        Path to kallisto abundance.tsv file
    ref_file: str
        Path to gene reference file
    annovar_file: str
        Path to Annovar results file
    output_path: str
        Path to save filtered results
    tpm_threshold: float, default=0
        TPM threshold for filtering genes
    """
    # Read reference file mapping transcript IDs to gene names
    input1 = pd.read_table(ref_file, names=['gene', 'EM'], sep=',')
    input1.dropna(axis=0, how='any', inplace=True)
    em_dic = dict(zip(input1['EM'], input1['gene']))
    
    # Read kallisto abundance file and filter by TPM threshold
    input2 = pd.read_table(kallisto_file, sep='\t')
    input2 = input2[input2['tpm'] > tpm_threshold]
    
    # Get gene list with TPM values
    gene_list = []
    for i in input2.index:
        target_id = input2.at[i, 'target_id']
        if target_id in em_dic:
            gene_list.append([em_dic[target_id], input2.at[i, 'tpm']])
    
    gene_df = pd.DataFrame(gene_list, columns=['gene', 'tpm'])
    gene_df.drop_duplicates(inplace=True)
    
    # Create dictionary of gene to TPM values for quick lookup
    tpm_dict = dict(zip(gene_df['gene'], gene_df['tpm']))
    
    # Read Annovar results and filter based on TPM values
    annovar_df = pd.read_table(annovar_file, sep='\t')
    
    # Print column names for debugging
    print("Columns in Annovar file:", list(annovar_df.columns))
    
    # Try to find the gene column in Annovar file - common names used in Annovar output
    gene_column_candidates = ['Gene', 'gene', 'Gene.refGene', 'Gene_refGene', 'Gene_knownGene', 'Gene_ensGene']
    gene_column = None
    
    for col in gene_column_candidates:
        if col in annovar_df.columns:
            gene_column = col
            print(f"Found gene column: {gene_column}")
            break
    
    if gene_column is None:
        print("Error: Could not find a gene column in the Annovar results file.")
        print("Please check the Annovar file format and specify the gene column manually.")
        sys.exit(1)
    
    # Function to check if any gene in a comma-separated list has TPM > threshold
    def has_gene_with_tpm_above_threshold(gene_str):
        if not pd.notnull(gene_str):
            return False
        
        # Split the gene string by comma, in case there are multiple genes
        genes = gene_str.split(',')
        
        # Check if any of the genes has TPM > threshold
        for gene in genes:
            gene = gene.strip()  # Remove any whitespace
            if gene in tpm_dict and tpm_dict[gene] > tpm_threshold:
                return True
        
        return False
    
    # Filter Annovar results to keep only genes with TPM > threshold
    filtered_annovar = annovar_df[annovar_df[gene_column].apply(has_gene_with_tpm_above_threshold)]
    
    # Save filtered results
    output_dir = os.path.dirname(output_path)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    filtered_annovar.to_csv(output_path, sep='\t', index=False)
    
    # Save gene expression file
    gene_expr_path = os.path.join(os.path.dirname(output_path), 'kallisto_gene_expression.txt')
    gene_df.to_csv(gene_expr_path, sep='\t', index=False)
    
    return filtered_annovar, gene_df

def main():
    """
    Main function to run the TPM filter with command line arguments
    """
    parser = argparse.ArgumentParser(description='Filter genes based on TPM values and compare with Annovar results')
    parser.add_argument('-k', '--kallisto', required=True, help='Path to kallisto abundance.tsv file')
    parser.add_argument('-r', '--reference', required=True, help='Path to gene reference file')
    parser.add_argument('-a', '--annovar', required=True, help='Path to Annovar results file')
    parser.add_argument('-o', '--output', required=True, help='Path to save filtered results')
    parser.add_argument('-t', '--threshold', type=float, default=0, help='TPM threshold for filtering genes (default: 0)')
    parser.add_argument('-g', '--gene-column', help='Name of the gene column in Annovar file (if not automatically detected)')
    
    args = parser.parse_args()
    
    # Make sure output directory exists
    output_dir = os.path.dirname(args.output)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    
    try:
        filtered_annovar, gene_df = tpm_filter_with_annovar(
            kallisto_file=args.kallisto,
            ref_file=args.reference,
            annovar_file=args.annovar,
            output_path=args.output,
            tpm_threshold=args.threshold
        )
        
        print(f"Filtering complete. Processed {len(gene_df)} genes.")
        print(f"Found {len(filtered_annovar)} genes in Annovar results with TPM > {args.threshold}")
        print(f"Results saved to {args.output}")
        print(f"Gene expression saved to {os.path.join(os.path.dirname(args.output), 'kallisto_gene_expression.txt')}")
    
    except Exception as e:
        print(f"Error occurred: {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()

