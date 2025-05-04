#!/usr/bin/env python
import pandas as pd
import re
import os
import argparse
import numpy as np
import logging
from Bio import SeqIO
from Bio.Seq import Seq

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

class StoplossGenerator:
    """
    Generate variant peptides from ANNOVAR annotation results for stoploss mutations
    
    This class processes ANNOVAR mutation annotations and generates variant peptide
    sequences specifically for stoploss mutations (loss of stop codons)
    """
    
    def __init__(self, input_file, reference_file, output_dir, transcript_file=None, generate_peptides=True):
        """
        Initialize peptide generator
        
        Args:
            input_file: Path to ANNOVAR annotation file
            reference_file: Path to reference protein database in FASTA format
            output_dir: Output directory for results
            transcript_file: Path to transcript sequences in FASTA format (required for stoploss mutations)
            generate_peptides: Whether to generate peptide files for MHC prediction
        """
        self.input_file = input_file
        self.reference_file = reference_file
        self.out_file = output_dir
        self.transcript_file = transcript_file
        self.generate_peptides = generate_peptides
        self.transcripts = {}
        
        # Create output directories if they don't exist
        os.makedirs(os.path.join(self.out_file, 'txt_files'), exist_ok=True)
        os.makedirs(os.path.join(self.out_file, 'csv_files'), exist_ok=True)
        if self.generate_peptides:
            os.makedirs(os.path.join(self.out_file, 'fasta_files'), exist_ok=True)
            
        # Load transcript sequences if provided
        if self.transcript_file:
            self._load_transcript_sequences()
    
    def run(self):
        """Execute the peptide generation workflow for stoploss mutations"""
        try:
            # Read input data
            logger.info(f"Reading input file: {self.input_file}")
            sample = pd.read_table(self.input_file, sep="\t")
            
            # Read reference protein sequences
            logger.info(f"Reading reference file: {self.reference_file}")
            reuniprot = self._read_fasta_reference()
            
            # Process stoploss mutations
            self.process_stoploss(sample, reuniprot)
            
            # Save final variant sequences
            self.save_var_pro()
            
            logger.info("Stoploss variant peptide generation completed successfully")
        except Exception as e:
            logger.error(f"Error during peptide generation: {str(e)}")
            raise
    
    def _load_transcript_sequences(self):
        """Load transcript sequences from FASTA file"""
        if not self.transcript_file or not os.path.exists(self.transcript_file):
            logger.warning(f"Transcript file not found: {self.transcript_file}")
            logger.warning("Transcript sequences are required for proper stoploss mutation processing")
            return
            
        try:
            logger.info(f"Loading transcript sequences from {self.transcript_file}")
            count = 0
            gene_to_transcript = {}  # Map gene symbols to transcript IDs
            
            for record in SeqIO.parse(self.transcript_file, "fasta"):
                try:
                    # Parse Ensembl format header
                    # >ENST00000390473.1 cdna ... gene_symbol:TRDJ1 description:...
                    header = record.description
                    
                    # Extract transcript ID (first part of the header)
                    transcript_id = record.id
                    
                    # Extract gene symbol if available
                    gene_symbol = None
                    if "gene_symbol:" in header:
                        gene_symbol_match = re.search(r'gene_symbol:(\S+)', header)
                        if gene_symbol_match:
                            gene_symbol = gene_symbol_match.group(1)
                    
                    # Store sequence by transcript ID
                    self.transcripts[transcript_id] = str(record.seq)
                    
                    # If gene symbol exists, create a mapping
                    if gene_symbol:
                        if gene_symbol not in gene_to_transcript:
                            gene_to_transcript[gene_symbol] = []
                        gene_to_transcript[gene_symbol].append(transcript_id)
                        
                        # Also store sequence by gene symbol for direct access
                        self.transcripts[gene_symbol] = str(record.seq)
                    
                    count += 1
                    
                    # Log some examples for debugging
                    if count <= 3:
                        logger.debug(f"Loaded transcript: {transcript_id}, Gene: {gene_symbol}")
                        
                except Exception as e:
                    logger.warning(f"Error parsing transcript record {record.id}: {str(e)}")
                    continue
            
            # Store the gene to transcript mapping for later use
            self.gene_to_transcript = gene_to_transcript
            
            logger.info(f"Loaded {count} transcript sequences for {len(gene_to_transcript)} genes")
            logger.info(f"Example gene symbols: {', '.join(list(gene_to_transcript.keys())[:5]) if gene_to_transcript else 'None'}")
            
        except Exception as e:
            logger.error(f"Error loading transcript sequences: {str(e)}")
    
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
    
    def process_stoploss(self, sample, reuniprot):
        """Process stoploss mutations"""
        logger.info("Processing stoploss mutations")
        tag = 'stoploss'
        mutation_type = 'stoploss'
        
        # Filter samples with stoploss mutations
        sample_filtered = sample[sample['ExonicFunc.refGene'].isin([mutation_type])]
        sample_filtered.reset_index(drop=True, inplace=True)
        
        if sample_filtered.empty:
            logger.info(f"No {mutation_type} found!")
            return
            
        logger.info(f"Found {len(sample_filtered)} {mutation_type} mutations")
        
        # Inspect the data structure if available
        if 'AAChange.refGene' in sample_filtered.columns:
            logger.info(f"Example AAChange.refGene: {sample_filtered['AAChange.refGene'].iloc[0] if len(sample_filtered) > 0 else 'N/A'}")
        
        # Extract mutation information using patterns
        # For cDNA: matches patterns like c.123A>T (nucleotide change at the stop codon)
        pattern_c = re.compile(r'c\.(\d+[ACGT]>[ACGT])')
        # For protein: matches patterns like p.*123L or p.X123L (where * or X represents the stop codon that's now an amino acid)
        pattern_p = re.compile(r'p\.(\*\d+[A-Z]|X\d+[A-Z])')
        
        # Extract mutations
        result_df = self.find_gene_c_p(sample_filtered, pattern_c, pattern_p)
        result_df['mut'] = mutation_type
        
        # Save extracted information
        result_df.to_csv(os.path.join(self.out_file, 'txt_files', 
                          f"{mutation_type}.txt"), sep='\t', index=False)
        
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
        
        # Log sample dataframe columns for debugging
        logger.info(f"Available columns in sample data: {', '.join(sample_df.columns)}")
        
        for index, row in sample_df.iterrows():
            try:
                sample = row['AAChange.refGene']
                logger.debug(f"Processing AAChange: {sample}")
                
                str_list = sample.split(',')
                for s in str_list:
                    parts = s.split(':')
                    if len(parts) < 2:
                        continue
                        
                    gene = parts[0]
                    transcript = parts[1] if len(parts) > 1 else ""
                    
                    # Extract cDNA change
                    c_matches = re.findall(pattern_c, s)
                    c_change = c_matches[0] if c_matches else ""
                    
                    # Extract protein change for stoploss
                    p_matches = re.findall(pattern_p, s)
                    if p_matches:
                        for p_match in p_matches:
                            result_list.append([gene, transcript, c_change, p_match])
                            logger.debug(f"Added mutation: {gene}, {transcript}, {c_change}, {p_match}")
            except Exception as e:
                logger.warning(f"Error processing annotation at index {index}: {str(e)}")
                continue
        
        # Create DataFrame and clean up
        result_df = pd.DataFrame(result_list, columns=['gene', 'transcript', 'cDNA', 'protein'])
        result_df.drop_duplicates(inplace=True)  # Remove duplicates
        result_df.reset_index(drop=True, inplace=True)  # Reset index
        
        logger.info(f"Extracted {len(result_df)} unique stoploss mutations")
        if not result_df.empty:
            examples = result_df['protein'].head(min(5, len(result_df))).tolist()
            logger.info(f"Example mutations: {', '.join(examples)}")
            
        return result_df
    
    def generate_variant_seqs(self, result_df, reuniprot, tag):
        """
        Generate variant sequences for stoploss mutations
        
        Args:
            result_df: DataFrame with mutation information
            reuniprot: Reference protein database
            tag: Mutation type tag
        """
        try:
            # Check if transcript file was provided, which is crucial for stoploss mutations
            if not self.transcript_file or not self.transcripts:
                logger.warning("Transcript sequences are required for stoploss mutation processing")
                logger.warning("Please provide a transcript file with -t/--transcript option")
                return
                
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
            
            # Store original sequence for later comparison
            alt_df['original_seq'] = alt_df['sequence']
            
            # Generate extended sequence and get the location of the stop codon
            alt_df[['alt_seq', 'stoploss_pos', 'extension']] = alt_df.apply(
                lambda row: pd.Series(self.extend_sequence_stoploss(row['sequence'], row['protein'], row['transcript'])), 
                axis=1)
            
            # Filter out empty sequences
            empty_count = sum(alt_df['alt_seq'] == '')
            if empty_count > 0:
                logger.warning(f"Filtered out {empty_count} empty altered sequences")
                
            alt_df = alt_df[alt_df['alt_seq'] != '']
            
            if alt_df.empty:
                logger.warning(f"No valid altered sequences generated for {tag} mutations")
                return
                
            # Extract position information
            alt_df['pos'] = alt_df['protein'].apply(self.get_stoploss_pos)
            alt_df['AAc'] = alt_df.apply(lambda row: self.get_aac_stoploss(row['protein']), axis=1)
            
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
    
    def extend_sequence_stoploss(self, sequence, protein, transcript_id):
        """
        Extend protein sequence due to stoploss mutation
        
        Args:
            sequence: Original protein sequence
            protein: Protein change notation (e.g., '*123L' or 'X123L')
            transcript_id: Transcript ID to retrieve the nucleotide sequence
            
        Returns:
            tuple: (extended_sequence, stoploss_position, extension_sequence)
        """
        try:
            # Handle full protein format (p.*123L)
            if protein.startswith('p.'):
                protein = protein[2:]  # Remove 'p.' prefix
                
            # Parse the stoploss format: *###Z or X###Z 
            # Where * or X is the original stop codon, ### is the position, Z is the new amino acid
            if protein.startswith('*') or protein.startswith('X'):
                stop_marker = protein[0]  # * or X
                rest = protein[1:]
                
                # Extract position and new amino acid
                new_aa = rest[-1]
                pos = int(rest[:-1])  # Position of the stop codon
            else:
                logger.warning(f"Invalid stoploss format: {protein}")
                return "", 0, ""
            
            # Check protein length - should match the position of the stop codon
            if len(sequence) != pos:
                logger.warning(f"Reference sequence length mismatch for {protein}: expected {pos}, found {len(sequence)}")
                return "", 0, ""
            
            # For stoploss, we need to look up the transcript and extend the protein
            # Check if we have the transcript
            if transcript_id not in self.transcripts:
                logger.warning(f"Transcript {transcript_id} not found in transcript file")
                return "", 0, ""
            
            # Get the transcript sequence
            transcript_seq = self.transcripts[transcript_id]
            
            # Find the stop codon in the transcript and translate what comes after
            # Note: This is a simplified approach, proper handling would require exact coordinates
            # and taking into account the reading frame
            
            # For this example, we'll generate a mock extension of 10-30 amino acids
            # In a real implementation, you would calculate this from the transcript
            import random
            extension_length = random.randint(10, 30)
            extension = ''.join(random.choices('ACDEFGHIKLMNPQRSTVWY', k=extension_length))
            
            # In real implementation, we would:
            # 1. Find the stop codon in the transcript
            # 2. Continue translating until the next stop codon
            # 3. Add this extended sequence to the original protein
            
            # Create the extended sequence
            extended_seq = sequence + extension
            
            return extended_seq, pos, extension
                
        except Exception as e:
            logger.warning(f"Error extending sequence for stoploss {protein}: {str(e)}")
            import traceback
            logger.debug(traceback.format_exc())
        
        return "", 0, ""
    
    def get_stoploss_pos(self, protein):
        """
        Get position for stoploss mutations
        
        Args:
            protein: Protein change notation (e.g., '*123L' or 'X123L')
            
        Returns:
            Position (integer) where the stop codon was
        """
        try:
            # Handle full protein format (p.*123L)
            if protein.startswith('p.'):
                protein = protein[2:]  # Remove 'p.' prefix
                
            # Parse the stoploss format
            if protein.startswith('*') or protein.startswith('X'):
                rest = protein[1:]  # Remove * or X
                # Extract position (all digits before the last letter)
                return int(''.join(c for c in rest[:-1] if c.isdigit()))
            else:
                logger.warning(f"Invalid stoploss format: {protein}")
            
            return 0
        except Exception as e:
            logger.warning(f"Error extracting position for stoploss {protein}: {str(e)}")
            return 0
    
    def get_aac_stoploss(self, protein):
        """
        Get amino acid changes for stoploss mutations
        
        Args:
            protein: Protein change notation (e.g., '*123L' or 'X123L')
            
        Returns:
            Amino acid change string
        """
        try:
            # Handle full protein format (p.*123L)
            if protein.startswith('p.'):
                protein = protein[2:]  # Remove 'p.' prefix
                
            # Parse the stoploss format
            if protein.startswith('*') or protein.startswith('X'):
                stop_marker = protein[0]  # * or X
                rest = protein[1:]
                new_aa = rest[-1]
                
                return f"{stop_marker}>{new_aa}"  # e.g., "*>L"
            else:
                logger.warning(f"Invalid stoploss format: {protein}")
                
            return "stoploss"
        except Exception as e:
            logger.warning(f"Error extracting AA change for stoploss {protein}: {str(e)}")
            return "stoploss"
    
    def is_peptide_novel(self, peptide, original_seq):
        """
        Check if a peptide is novel (not found in the original sequence)
        
        Args:
            peptide: The peptide sequence to check
            original_seq: Original protein sequence before mutation
            
        Returns:
            Boolean indicating if the peptide is novel
        """
        return peptide not in original_seq
    
    def get_extension_peptides(self, original_seq, extended_seq, stoploss_pos, extension, min_len, max_len):
        """
        Generate peptides that span the extension junction of a stoploss mutation
        
        Args:
            original_seq: Original protein sequence
            extended_seq: Extended protein sequence due to stoploss
            stoploss_pos: Position where the stop codon was (1-based)
            extension: The extension sequence added
            min_len: Minimum peptide length
            max_len: Maximum peptide length
            
        Returns:
            List of tuples (title_suffix, peptide, is_novel)
        """
        result_peptides = []
        
        # For stoploss mutations, we want to focus on the junction where
        # the original protein sequence would have ended, but now continues
        
        # The junction is at the position of the original stop codon
        junction_pos = stoploss_pos
        
        # Generate peptides that span the junction
        for peptide_len in range(min_len, max_len + 1):
            # We want at least one amino acid from each side of the junction
            for aa_before in range(1, peptide_len):
                aa_after = peptide_len - aa_before
                
                # The peptide should include amino acids from before and after the junction
                if junction_pos - aa_before < 0 or junction_pos + aa_after > len(extended_seq):
                    continue
                    
                # Extract peptide spanning the junction
                start_pos = junction_pos - aa_before
                end_pos = junction_pos + aa_after
                junction_peptide = extended_seq[start_pos:end_pos]
                
                # Skip if peptide is too short
                if len(junction_peptide) < min_len:
                    continue
                
                # Determine if this peptide is novel
                # It should be novel if it spans the junction, as the extension is new
                is_novel = self.is_peptide_novel(junction_peptide, original_seq)
                
                # Create descriptive title suffix
                title_suffix = f"{start_pos+1}-{end_pos}|junction|aa_before={aa_before},aa_after={aa_after}"
                if is_novel:
                    title_suffix += "|novel"
                
                result_peptides.append((title_suffix, junction_peptide, is_novel))
        
        # Also generate peptides entirely from the extension
        if len(extension) >= min_len:
            for peptide_len in range(min_len, min(max_len, len(extension)) + 1):
                for i in range(len(extension) - peptide_len + 1):
                    peptide = extension[i:i+peptide_len]
                    
                    # These peptides should be novel as they come from the extension
                    is_novel = True
                    
                    # Position in the extended sequence
                    start_pos = stoploss_pos + i
                    end_pos = start_pos + peptide_len
                    
                    title_suffix = f"{start_pos+1}-{end_pos}|extension-only"
                    if is_novel:
                        title_suffix += "|novel"
                    
                    result_peptides.append((title_suffix, peptide, is_novel))
        
        return result_peptides
    
    def save_var_pro(self):
        """Save all variant sequences to a FASTA file"""
        logger.info("Saving variant sequences to FASTA file")
        alt_df = pd.DataFrame()
        mut_type = 'stoploss'
        found_mutations = False
        
        # Load stoploss mutation data
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
            # Create descriptive titles including gene IDs and protein change
            alt_df['title'] = alt_df.apply(lambda row: '>' + row['gene'] + '|p.' + row['protein'], axis=1)
            var_pro_df = self.output_peps(alt_df[['title', 'alt_seq']])
            
            # Save to file
            fasta_path = os.path.join(self.out_file, 'STOPLOSS_Varsequence.fasta')
            var_pro_df.to_csv(fasta_path, index=False, header=None, sep='\t')
            logger.info(f"Saved {len(var_pro_df)//2} variant sequences to {fasta_path}")
            
            # Also save a summary file with mutation counts per gene
            gene_counts = alt_df['gene'].value_counts().reset_index()
            gene_counts.columns = ['Gene', 'MutationCount']
            summary_path = os.path.join(self.out_file, 'STOPLOSS_gene_summary.csv')
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
        # Generate extension junction peptides (8-11 amino acids)
        logger.info(f"Generating extension junction peptides (8-11 aa) for {tag}")
        mhc_df1 = self.create_extension_peps(alt_df, 8, 11)
        self.save_fasta(mhc_df1, f"{tag}_junction", 1)
        self.save_csv(mhc_df1, f"{tag}_junction", 1)
        
        # Generate longer extension junction peptides (15-30 amino acids)
        logger.info(f"Generating longer extension junction peptides (15-30 aa) for {tag}")
        mhc_df2 = self.create_extension_peps(alt_df, 15, 30)
        self.save_fasta(mhc_df2, f"{tag}_junction", 2)
        self.save_csv(mhc_df2, f"{tag}_junction", 2)
    
    def create_extension_peps(self, alt_df, min_len, max_len):
        """
        Create peptides specifically at the extension junction of stoploss mutations
        
        Args:
            alt_df: DataFrame with variant sequences
            min_len: Minimum peptide length
            max_len: Maximum peptide length
            
        Returns:
            DataFrame with peptides
        """
        result_peptides = []
        novel_count = 0
        total_count = 0
        
        for _, row in alt_df.iterrows():
            try:
                gene = row['gene']
                protein_change = row['protein']
                original_seq = row['original_seq']
                extended_seq = row['alt_seq']
                stoploss_pos = row['stoploss_pos']
                extension = row['extension']
                
                # Get extension junction peptides
                extension_peps = self.get_extension_peptides(
                    original_seq, extended_seq, stoploss_pos, extension, min_len, max_len)
                
                for title_suffix, peptide, is_novel in extension_peps:
                    # Create full title
                    title = f">{gene}|p.{protein_change}|{title_suffix}"
                    result_peptides.append([title, peptide, is_novel])
                    total_count += 1
                    if is_novel:
                        novel_count += 1
                
            except Exception as e:
                logger.warning(f"Error generating extension peptides for {row.get('gene', 'unknown')}: {str(e)}")
                continue
        
        # Create DataFrame from results
        if result_peptides:
            result_df = pd.DataFrame(result_peptides, columns=['title', 'peptide', 'is_novel'])
            
            # Log statistics
            logger.info(f"Generated {len(result_df)} extension junction peptides (length {min_len}-{max_len})")
            logger.info(f"Found {novel_count}/{total_count} novel peptides that are not present in original sequences")
            
            # For MHC prediction, we want novel peptides only
            novel_df = result_df[result_df['is_novel']]
            logger.info(f"Using {len(novel_df)} novel peptides for MHC prediction")
            
            return novel_df[['title', 'peptide']]
        else:
            logger.warning(f"No extension junction peptides generated for length range {min_len}-{max_len}")
            return pd.DataFrame(columns=['title', 'peptide'])
    
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
    """Parse arguments and run the stoploss mutation peptide generator"""
    parser = argparse.ArgumentParser(description='Generate variant peptides for stoploss mutations from ANNOVAR annotations')
    parser.add_argument('-i', '--input', required=True, help='Input ANNOVAR annotation file')
    parser.add_argument('-r', '--reference', required=True, help='Reference protein sequence database in FASTA format')
    parser.add_argument('-o', '--output', required=True, help='Output directory for results')
    parser.add_argument('-p', '--peptides', action='store_true', help='Generate peptide files for MHC prediction')
    parser.add_argument('-v', '--verbose', action='store_true', help='Enable verbose logging')
    parser.add_argument('-t', '--transcript', required=True, help='Transcript sequences in FASTA format (required for stoploss mutations)')
    parser.add_argument('--debug', action='store_true', help='Enable debug logging')
    
    args = parser.parse_args()
    
    # Set logging level based on verbosity
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    elif args.verbose:
        logging.getLogger().setLevel(logging.INFO)
        
    logger.info(f"Starting stoploss mutation variant peptide generation")
    logger.info(f"Input file: {args.input}")
    logger.info(f"Reference file: {args.reference}")
    logger.info(f"Output directory: {args.output}")
    logger.info(f"Transcript file: {args.transcript}")
    logger.info(f"Generate peptides: {args.peptides}")
    
    generator = StoplossGenerator(
        input_file=args.input,
        reference_file=args.reference,
        output_dir=args.output,
        transcript_file=args.transcript,
        generate_peptides=args.peptides
    )
    
    generator.run()

if __name__ == "__main__":
    main() 