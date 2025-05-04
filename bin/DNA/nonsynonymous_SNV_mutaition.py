#!/usr/bin/env python
import pandas as pd
import re
import os
import argparse
import numpy as np
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

class SNVPeptideGenerator:
    """
    Generate variant peptides from ANNOVAR annotation results for SNVs
    
    This class processes ANNOVAR mutation annotations and generates variant peptide
    sequences specifically for nonsynonymous single nucleotide variants (SNVs)
    """
    
    def __init__(self, input_file, reference_file, output_dir, generate_peptides=True):
        """
        Initialize peptide generator
        
        Args:
            input_file: Path to ANNOVAR annotation file
            reference_file: Path to reference protein database in FASTA format
            output_dir: Output directory for results
            generate_peptides: Whether to generate peptide files for MHC prediction
        """
        self.input_file = input_file
        self.reference_file = reference_file
        self.out_file = output_dir
        self.generate_peptides = generate_peptides
        
        # Create output directories if they don't exist
        os.makedirs(os.path.join(self.out_file, 'txt_files'), exist_ok=True)
        os.makedirs(os.path.join(self.out_file, 'csv_files'), exist_ok=True)
        if self.generate_peptides:
            os.makedirs(os.path.join(self.out_file, 'fasta_files'), exist_ok=True)
    
    def run(self):
        """Execute the peptide generation workflow for SNVs"""
        try:
            # Read input data
            logger.info(f"Reading input file: {self.input_file}")
            sample = pd.read_table(self.input_file, sep="\t")
            
            # Read reference protein sequences
            logger.info(f"Reading reference file: {self.reference_file}")
            reuniprot = self._read_fasta_reference()
            
            # Process SNV mutations
            self.process_snvs(sample, reuniprot)
            
            # Save final variant sequences
            self.save_var_pro()
            
            logger.info("SNV variant peptide generation completed successfully")
        except Exception as e:
            logger.error(f"Error during peptide generation: {str(e)}")
            raise
    
    def _read_fasta_reference(self):
        """
        Read FASTA format reference file
        
        Returns:
            DataFrame with gene and sequence columns
        """
        logger.info(f"Parsing FASTA reference file: {self.reference_file}")
        
        all_entries = []
        current_genes = []
        current_seq = []
        
        try:
            with open(self.reference_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                        
                    if line.startswith('>'):
                        # Save previous sequence if exists
                        if current_genes and current_seq:
                            sequence = ''.join(current_seq)
                            for gene in current_genes:
                                all_entries.append((gene, sequence))
                        
                        # Start new sequence
                        header = line[1:]  # Remove '>'
                        current_genes = []
                        
                        # Extract all gene identifiers from the header
                        # Format can be: >GeneA|GeneB or >GeneA | GeneB
                        for part in header.split('|'):
                            genes = [g.strip() for g in part.split() if g.strip()]
                            current_genes.extend(genes)
                            
                        logger.debug(f"Found genes in header: {current_genes}")
                        current_seq = []
                    else:
                        # Append sequence line
                        current_seq.append(line)
            
            # Add the last sequence
            if current_genes and current_seq:
                sequence = ''.join(current_seq)
                for gene in current_genes:
                    all_entries.append((gene, sequence))
                
            # Create DataFrame
            ref_df = pd.DataFrame(all_entries, columns=['gene', 'sequence'])
            ref_df.drop_duplicates(inplace=True)  # Remove duplicate gene-sequence pairs
            
            logger.info(f"Loaded {len(ref_df)} reference sequences for {len(set(ref_df['gene']))} genes")
            return ref_df
            
        except Exception as e:
            logger.error(f"Error reading reference file: {str(e)}")
            raise
    
    def process_snvs(self, sample, reuniprot):
        """Process nonsynonymous SNV mutations"""
        logger.info("Processing nonsynonymous SNV mutations")
        tag = 'snv'
        mutation_type = 'nonsynonymous SNV'
        
        # Filter samples with SNV mutations
        sample_filtered = sample[sample['ExonicFunc.refGene'].isin([mutation_type])]
        sample_filtered.reset_index(drop=True, inplace=True)
        
        if sample_filtered.empty:
            logger.info(f"No {mutation_type} found!")
            return
        
        logger.info(f"Found {len(sample_filtered)} {mutation_type} mutations")    
        
        # Extract mutation information
        pattern_c = re.compile(r'c\.([A-Z]\d{1,8}[A-Z]):')
        pattern_p = re.compile(r'p\.([A-Z]\d{1,8}[A-Z])')
        result_df = self.find_gene_c_p(sample_filtered, pattern_c, pattern_p)
        result_df['mut'] = mutation_type
        
        # Save extracted information
        result_df.to_csv(os.path.join(self.out_file, 'txt_files', 
                          f"{mutation_type.replace(' ', '_')}.txt"), sep='\t', index=False)
        
        # Find sequences and generate variants
        self.generate_variant_seqs(result_df, reuniprot, tag)
    
    def find_gene_c_p(self, sample_df, pattern_c, pattern_p):
        """
        Extract gene, cDNA change and protein change from AAChange.refGene column
        
        Args:
            sample_df: DataFrame with ANNOVAR annotations
            pattern_c: Regex pattern for cDNA change
            pattern_p: Regex pattern for protein change
            
        Returns:
            DataFrame with gene, cDNA, protein columns
        """
        result_list = []
        for sample in sample_df['AAChange.refGene']:
            try:
                str_list = sample.split(',')
                for s in str_list:
                    gene = s.split(':')[0]
                    c = re.findall(pattern_c, s)
                    p = re.findall(pattern_p, s)
                    if len(c) and len(p) and len(c) == len(p):
                        for i in range(len(c)):
                            result_list.append([gene, c[i], p[i]])
            except Exception as e:
                logger.warning(f"Error processing annotation: {sample}, error: {str(e)}")
                continue
        
        result_df = pd.DataFrame(result_list, columns=['gene', 'cDNA', 'protein'])
        result_df.drop_duplicates(inplace=True)  # Remove duplicates
        result_df.reset_index(drop=True, inplace=True)  # Reset index
        return result_df
    
    def generate_variant_seqs(self, result_df, reuniprot, tag):
        """
        Generate variant sequences for SNV mutations
        
        Args:
            result_df: DataFrame with mutation information
            reuniprot: Reference protein database
            tag: Mutation type tag (snv)
        """
        try:
            # Join mutation info with reference sequences
            merged_df = pd.merge(result_df, reuniprot, how='inner', left_on='gene', right_on='gene')
            
            if merged_df.empty:
                logger.warning(f"No reference sequences found for {tag} mutations")
                unmatched_genes = set(result_df['gene']) - set(reuniprot['gene'])
                if unmatched_genes:
                    logger.warning(f"Genes without matching reference sequences: {', '.join(list(unmatched_genes)[:10])}" +
                                  (f" and {len(unmatched_genes) - 10} more" if len(unmatched_genes) > 10 else ""))
                return
                
            logger.info(f"Generating variant sequences for {len(merged_df)} {tag} mutations")
            logger.info(f"Found {merged_df['gene'].nunique()} unique genes with mutations")
            
            # Create altered sequences
            alt_df = merged_df.copy()
            alt_df['alt_seq'] = alt_df[['sequence', 'protein']].apply(
                lambda row: self.alter_seq_snv(row['sequence'], row['protein']), axis=1)
            
            # Filter out empty sequences
            empty_count = sum(alt_df['alt_seq'] == '')
            if empty_count > 0:
                logger.warning(f"Filtered out {empty_count} empty altered sequences")
                
            alt_df = alt_df[alt_df['alt_seq'] != '']
            
            if alt_df.empty:
                logger.warning(f"No valid altered sequences generated for {tag} mutations")
                return
                
            # Extract position and amino acid change information
            alt_df['pos'] = alt_df['protein'].apply(lambda p: int(p[1:-1]))
            alt_df['AAc'] = alt_df['protein'].apply(lambda p: p[0] + p[-1])
            
            # Save results
            out_path = os.path.join(self.out_file, 'csv_files', f'alt_{tag}.csv')
            alt_df.to_csv(out_path, index=False)
            logger.info(f"Saved {len(alt_df)} {tag} variant sequences to {out_path}")
            
            # Generate peptides if requested
            if self.generate_peptides:
                self.peps_save(alt_df, tag)
                logger.info(f"{tag} MHC peptide files have been saved")
                
        except Exception as e:
            logger.error(f"Error generating variant sequences for {tag}: {str(e)}")
            raise
    
    def alter_seq_snv(self, sequence, protein):
        """
        Alter sequence for SNVs
        
        Args:
            sequence: Original protein sequence
            protein: Protein change notation (e.g., 'A123B')
            
        Returns:
            Altered sequence
        """
        try:
            original_aa = protein[0]
            new_aa = protein[-1]
            position = int(protein[1:-1])
            
            if position <= len(sequence):
                if sequence[position-1] == original_aa:
                    return sequence[:position-1] + new_aa + sequence[position:]
                else:
                    logger.warning(f"Reference AA mismatch for {protein}: expected {original_aa}, found {sequence[position-1]}")
            else:
                logger.warning(f"Position {position} out of range for sequence length {len(sequence)}")
        except Exception as e:
            logger.warning(f"Error altering sequence for SNV {protein}: {str(e)}")
        
        return ""
    
    def save_var_pro(self):
        """Save all variant sequences to a FASTA file"""
        logger.info("Saving variant sequences to FASTA file")
        alt_df = pd.DataFrame()
        mut_type = 'snv'
        found_mutations = False
        
        # Load SNV mutation data
        try:
            filepath = os.path.join(self.out_file, 'csv_files', f'alt_{mut_type}.csv')
            if os.path.exists(filepath):
                alt_df = pd.read_csv(filepath)
                found_mutations = True
                logger.info(f"Loaded {len(alt_df)} {mut_type} mutations for variant sequence file")
                logger.info(f"Found {alt_df['gene'].nunique()} unique genes with mutations")
            else:
                logger.info(f"No {mut_type} mutations found (file not found)")
        except Exception as e:
            logger.warning(f"Error reading {mut_type} mutations: {str(e)}")
                
        if found_mutations and not alt_df.empty:
            # Create FASTA format
            # Create descriptive titles including both gene IDs and protein change
            alt_df['title'] = alt_df.apply(lambda row: '>' + row['gene'] + '|p.' + row['protein'], axis=1)
            var_pro_df = self.output_peps(alt_df[['title', 'alt_seq']])
            
            # Save to file
            fasta_path = os.path.join(self.out_file, 'SNV_Varsequence.fasta')
            var_pro_df.to_csv(fasta_path, index=False, header=None, sep='\t')
            logger.info(f"Saved {len(var_pro_df)//2} variant sequences to {fasta_path}")
            
            # Also save a summary file with mutation counts per gene
            gene_counts = alt_df['gene'].value_counts().reset_index()
            gene_counts.columns = ['Gene', 'MutationCount']
            summary_path = os.path.join(self.out_file, 'SNV_gene_summary.csv')
            gene_counts.to_csv(summary_path, index=False)
            logger.info(f"Saved gene mutation summary to {summary_path}")
        else:
            logger.warning("No variant sequences found to save")
    
    def peps_save(self, alt_df, tag):
        """
        Generate and save peptides for MHC prediction
        
        Args:
            alt_df: DataFrame with variant sequences
            tag: Mutation type tag
        """
        # Generate short peptides (8-11 amino acids)
        logger.info(f"Generating short peptides (8-11 aa) for {tag}")
        mhc_df1 = self.create_peps(alt_df, 8, 11)
        self.save_fasta(mhc_df1, tag)
        self.save_csv(mhc_df1, tag)
        
        # Generate long peptides (15-30 amino acids)
        logger.info(f"Generating long peptides (15-30 aa) for {tag}")
        mhc_df2 = self.create_peps(alt_df, 15, 30)
        self.save_fasta(mhc_df2, tag, 2)
        self.save_csv(mhc_df2, tag, 2)
    
    def create_peps(self, alt_df, min_len, max_len):
        """
        Create peptides of specified length range
        
        Args:
            alt_df: DataFrame with variant sequences
            min_len: Minimum peptide length
            max_len: Maximum peptide length
            
        Returns:
            DataFrame with peptides
        """
        result_peptides = []
        
        for _, row in alt_df.iterrows():
            try:
                gene = row['gene']
                protein_change = row['protein']
                sequence = row['alt_seq']
                position = row['pos']
                
                # Focus on the region around the mutation
                window_start = max(0, position - max_len)
                window_end = min(len(sequence), position + max_len)
                
                # Generate peptides of different lengths
                for peptide_len in range(min_len, max_len + 1):
                    # Slide through the window
                    for i in range(window_start, window_end - peptide_len + 1):
                        # Check if peptide contains the mutation site
                        if i <= position - 1 < i + peptide_len:
                            peptide = sequence[i:i+peptide_len]
                            title = f">{gene}|p.{protein_change}|{i+1}-{i+peptide_len}"
                            result_peptides.append([title, peptide])
            except Exception as e:
                logger.warning(f"Error generating peptides for {row.get('gene', 'unknown')}: {str(e)}")
                continue
        
        # Create DataFrame
        return pd.DataFrame(result_peptides, columns=['title', 'peptide'])
    
    def save_fasta(self, mhc_df, tag, suffix=1):
        """
        Save peptides to FASTA format
        
        Args:
            mhc_df: DataFrame with peptides
            tag: Mutation type tag
            suffix: Suffix for output filename
        """
        if mhc_df.empty:
            logger.warning(f"No peptides to save for {tag} (suffix {suffix})")
            return
            
        try:
            # Prepare output file path
            fasta_path = os.path.join(
                self.out_file, 
                'fasta_files', 
                f"{tag}_peptides_{suffix}.fasta"
            )
            
            # Format as FASTA
            with open(fasta_path, 'w') as f:
                for _, row in mhc_df.iterrows():
                    f.write(f"{row['title']}\n{row['peptide']}\n")
                    
            logger.info(f"Saved {len(mhc_df)} peptides to {fasta_path}")
        except Exception as e:
            logger.error(f"Error saving FASTA for {tag}: {str(e)}")
    
    def save_csv(self, mhc_df, tag, suffix=1):
        """
        Save peptides to CSV format
        
        Args:
            mhc_df: DataFrame with peptides
            tag: Mutation type tag
            suffix: Suffix for output filename
        """
        if mhc_df.empty:
            logger.warning(f"No peptides to save for {tag} (suffix {suffix})")
            return
            
        try:
            # Prepare output file path
            csv_path = os.path.join(
                self.out_file, 
                'csv_files', 
                f"{tag}_peptides_{suffix}.csv"
            )
            
            # Save as CSV
            mhc_df.to_csv(csv_path, index=False)
            logger.info(f"Saved {len(mhc_df)} peptides to {csv_path}")
        except Exception as e:
            logger.error(f"Error saving CSV for {tag}: {str(e)}")
    
    def output_peps(self, STR_MHC_df):
        """
        Format peptides for output
        
        Args:
            STR_MHC_df: DataFrame with titles and sequences
            
        Returns:
            Formatted DataFrame for FASTA output
        """
        try:
            # Remove duplicates
            STR_MHC_df = STR_MHC_df.drop_duplicates().reset_index(drop=True)
            
            # Convert to list and flatten
            flat_list = []
            for title, seq in STR_MHC_df.values:
                flat_list.append(title)
                flat_list.append(seq)
                
            # Create output DataFrame
            return pd.DataFrame(flat_list, columns=['peptides'])
        except Exception as e:
            logger.error(f"Error formatting peptides for output: {str(e)}")
            return pd.DataFrame(columns=['peptides'])

def main():
    """Parse arguments and run the SNV peptide generator"""
    parser = argparse.ArgumentParser(description='Generate variant peptides for SNVs from ANNOVAR annotations')
    parser.add_argument('-i', '--input', required=True, help='Input ANNOVAR annotation file')
    parser.add_argument('-r', '--reference', required=True, help='Reference protein sequence database in FASTA format')
    parser.add_argument('-o', '--output', required=True, help='Output directory for results')
    parser.add_argument('-p', '--peptides', action='store_true', help='Generate peptide files for MHC prediction')
    parser.add_argument('-v', '--verbose', action='store_true', help='Enable verbose logging')
    
    args = parser.parse_args()
    
    # Set logging level based on verbosity
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
        
    logger.info(f"Starting SNV variant peptide generation")
    logger.info(f"Input file: {args.input}")
    logger.info(f"Reference file: {args.reference}")
    logger.info(f"Output directory: {args.output}")
    logger.info(f"Generate peptides: {args.peptides}")
    
    generator = SNVPeptideGenerator(
        input_file=args.input,
        reference_file=args.reference,
        output_dir=args.output,
        generate_peptides=args.peptides
    )
    
    generator.run()

if __name__ == "__main__":
    main()
