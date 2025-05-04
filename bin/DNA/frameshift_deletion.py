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

class DelPeptideGenerator:
    """
    Generate variant peptides from ANNOVAR annotation results for frameshift deletions
    
    This class processes ANNOVAR mutation annotations and generates variant peptide
    sequences specifically for frameshift deletion mutations
    """
    
    def __init__(self, input_file, reference_file, output_dir, transcript_file=None, generate_peptides=True, 
                use_real_translation=False):
        """
        Initialize peptide generator
        
        Args:
            input_file: Path to ANNOVAR annotation file
            reference_file: Path to reference protein database in FASTA format
            output_dir: Output directory for results
            transcript_file: Path to transcript sequences in FASTA format (for real frameshift translation)
            generate_peptides: Whether to generate peptide files for MHC prediction
            use_real_translation: Whether to use real translation for frameshift mutations instead of X
        """
        self.input_file = input_file
        self.reference_file = reference_file
        self.out_file = output_dir
        self.transcript_file = transcript_file
        self.generate_peptides = generate_peptides
        self.use_real_translation = use_real_translation
        self.transcripts = {}
        
        # Create output directories if they don't exist
        os.makedirs(os.path.join(self.out_file, 'txt_files'), exist_ok=True)
        os.makedirs(os.path.join(self.out_file, 'csv_files'), exist_ok=True)
        if self.generate_peptides:
            os.makedirs(os.path.join(self.out_file, 'fasta_files'), exist_ok=True)
            
        # Load transcript sequences if provided and real translation is requested
        if self.transcript_file and self.use_real_translation:
            self._load_transcript_sequences()
    
    def run(self):
        """Execute the peptide generation workflow for frameshift deletions"""
        try:
            # Read input data
            logger.info(f"Reading input file: {self.input_file}")
            sample = pd.read_table(self.input_file, sep="\t")
            
            # Read reference protein sequences
            logger.info(f"Reading reference file: {self.reference_file}")
            reuniprot = self._read_fasta_reference()
            
            # Process frameshift deletion mutations
            self.process_deletions(sample, reuniprot)
            
            # Save final variant sequences
            self.save_var_pro()
            
            logger.info("Frameshift deletion variant peptide generation completed successfully")
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
    
    def process_deletions(self, sample, reuniprot):
        """Process frameshift deletion mutations"""
        logger.info("Processing frameshift deletion mutations")
        tag = 'fs_del'
        mutation_type = 'frameshift deletion'
        
        # Filter samples with deletion mutations
        sample_filtered = sample[sample['ExonicFunc.refGene'].isin([mutation_type])]
        sample_filtered.reset_index(drop=True, inplace=True)
        
        if sample_filtered.empty:
            logger.info(f"No {mutation_type} found!")
            return
            
        logger.info(f"Found {len(sample_filtered)} {mutation_type} mutations")
        
        # Inspect the data structure if available
        if 'AAChange.refGene' in sample_filtered.columns:
            logger.info(f"Example AAChange.refGene: {sample_filtered['AAChange.refGene'].iloc[0] if len(sample_filtered) > 0 else 'N/A'}")
        
        # Extract mutation information using more flexible pattern
        pattern_c = re.compile(r'c\.(\d+_\d+|\d+)del[A-Z]*:')
        # No need for protein regex as we'll extract directly in find_gene_c_p
        pattern_p = re.compile(r'p\.([A-Z]\d+[A-Z]fs\*\d+)')
        
        # Extract mutations
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
                    
                    # Extract cDNA change
                    c_matches = re.findall(pattern_c, s)
                    c_change = c_matches[0] if c_matches else ""
                    
                    # Extract protein change for frameshift
                    # Look specifically for patterns like p.Q997Vfs*95
                    p_match = None
                    if 'p.' in s and 'fs' in s:
                        # Extract everything after p. until the next : or end of string
                        p_part_match = re.search(r'p\.([^:,]+)', s)
                        if p_part_match:
                            p_part = p_part_match.group(1)
                            # Check if it's a frameshift
                            if 'fs' in p_part:
                                p_match = p_part
                                logger.debug(f"Found frameshift: {p_match}")
                    
                    if p_match:
                        result_list.append([gene, c_change, p_match])
                        logger.debug(f"Added mutation: {gene}, {c_change}, {p_match}")
            except Exception as e:
                logger.warning(f"Error processing annotation at index {index}: {str(e)}")
                continue
        
        # Create DataFrame and clean up
        result_df = pd.DataFrame(result_list, columns=['gene', 'cDNA', 'protein'])
        result_df.drop_duplicates(inplace=True)  # Remove duplicates
        result_df.reset_index(drop=True, inplace=True)  # Reset index
        
        logger.info(f"Extracted {len(result_df)} unique frameshift mutations")
        if not result_df.empty:
            examples = result_df['protein'].head(min(5, len(result_df))).tolist()
            logger.info(f"Example mutations: {', '.join(examples)}")
            
        return result_df
    
    def generate_variant_seqs(self, result_df, reuniprot, tag):
        """
        Generate variant sequences for frameshift deletion mutations
        
        Args:
            result_df: DataFrame with mutation information
            reuniprot: Reference protein database
            tag: Mutation type tag (fs_del)
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
            
            # Store original sequence for later comparison
            alt_df['original_seq'] = alt_df['sequence']
            
            # Apply the alter_seq_del_with_boundaries function to each row, with a transcript ID if available
            def alter_with_transcript(row):
                # Store the gene as the current transcript ID for use in alter_seq_del_with_boundaries
                self.current_transcript_id = row['gene']
                return pd.Series(self.alter_seq_del_with_boundaries(row['sequence'], row['protein']))
                
            alt_df[['alt_seq', 'del_start', 'del_end']] = alt_df.apply(
                lambda row: alter_with_transcript(row), axis=1)
            
            # Filter out empty sequences
            empty_count = sum(alt_df['alt_seq'] == '')
            if empty_count > 0:
                logger.warning(f"Filtered out {empty_count} empty altered sequences")
                
            alt_df = alt_df[alt_df['alt_seq'] != '']
            
            if alt_df.empty:
                logger.warning(f"No valid altered sequences generated for {tag} mutations")
                return
                
            # Extract position and amino acid change information
            alt_df['pos'] = alt_df['protein'].apply(self.get_del_pos)
            alt_df['AAc'] = alt_df.apply(lambda row: self.get_aac_del(row['protein'], row['alt_seq']), axis=1)
            
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
    
    def alter_seq_del_with_boundaries(self, sequence, protein):
        """
        Alter sequence for frameshift deletions and return deletion boundaries
        
        Args:
            sequence: Original protein sequence
            protein: Protein change notation (e.g., 'Q997Vfs*95' or 'p.Q997Vfs*95')
            
        Returns:
            tuple: (altered_sequence, deletion_start, deletion_end)
        """
        # Log the input for debugging
        logger.debug(f"Attempting to alter sequence for frameshift: {protein}")
        
        try:
            # Handle full protein format (p.Q997Vfs*95)
            if protein.startswith('p.'):
                protein = protein[2:]  # Remove 'p.' prefix
            
            # First try standard format with *
            match = re.match(r'([A-Z])(\d+)([A-Z])fs\*(\d+)', protein)
            
            if not match:
                # Try without * (sometimes annotation might use different format)
                match = re.match(r'([A-Z])(\d+)([A-Z])fs(\d*)', protein)
                
            if not match:
                logger.warning(f"Invalid frameshift format: {protein}")
                return "", 0, 0
            
            # Extract components
            groups = match.groups()
            if len(groups) == 4:
                orig_aa, pos_str, new_aa, term_len_str = groups
                term_len = int(term_len_str) if term_len_str else 10
            elif len(groups) == 3:
                orig_aa, pos_str, new_aa = groups
                term_len = 10  # Default
            else:
                logger.warning(f"Unexpected match groups for {protein}: {groups}")
                return "", 0, 0
                
            pos = int(pos_str)
            
            # Check position validity
            if pos > len(sequence):
                logger.warning(f"Position {pos} out of range for frameshift: sequence length {len(sequence)}")
                return "", 0, 0
            
            # Validate reference amino acid
            if sequence[pos-1] != orig_aa:
                logger.warning(f"Reference AA mismatch for {protein}: expected {orig_aa}, found {sequence[pos-1]}")
                return "", 0, 0
            
            # Create modified sequence up to the mutation point
            mod_seq = sequence[:pos-1] + new_aa
            
            if self.use_real_translation and hasattr(self, 'current_transcript_id') and self.current_transcript_id in self.transcripts:
                # If we have the transcript, try to translate the actual frameshift sequence
                frameshift_seq = self.translate_frameshift(self.current_transcript_id, pos, term_len)
                if frameshift_seq:
                    # Start with the new amino acid and add the translated frameshift sequence
                    mod_seq = sequence[:pos-1] + new_aa + frameshift_seq
                    logger.debug(f"Using real translation for {protein}: {frameshift_seq}")
                else:
                    # Fall back to placeholder Xs if translation failed
                    for i in range(min(term_len, 30)):
                        mod_seq += "X"
            else:
                # Add placeholder amino acids for frameshift sequence
                for i in range(min(term_len, 30)):
                    mod_seq += "X"  # Use 'X' to represent unknown/modified amino acids
            
            logger.debug(f"Successfully created altered sequence for {protein}, length: {len(mod_seq)}")
            return mod_seq, pos, pos  # For frameshift, we consider the mutation to be at one position
            
        except Exception as e:
            logger.warning(f"Error altering sequence for frameshift {protein}: {str(e)}")
            import traceback
            logger.debug(traceback.format_exc())
        
        return "", 0, 0
    
    def alter_seq_del(self, sequence, protein):
        """
        Alter sequence for frameshift deletions
        
        Args:
            sequence: Original protein sequence
            protein: Protein change notation (e.g., 'Q997Vfs*95' or 'p.Q997Vfs*95')
            
        Returns:
            Altered sequence with frameshift
        """
        result, _, _ = self.alter_seq_del_with_boundaries(sequence, protein)
        return result
    
    def _load_transcript_sequences(self):
        """Load transcript sequences from FASTA file"""
        if not self.transcript_file or not os.path.exists(self.transcript_file):
            logger.warning(f"Transcript file not found: {self.transcript_file}")
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
    
    def translate_frameshift(self, transcript_id, aa_pos, term_len):
        """
        Translate frameshift sequence based on transcript and amino acid position
        
        Args:
            transcript_id: Transcript ID or gene symbol
            aa_pos: Amino acid position (1-based)
            term_len: Expected termination length
            
        Returns:
            Translated frameshift sequence, or None if translation failed
        """
        try:
            # First try direct match with transcript_id
            if transcript_id in self.transcripts:
                transcript_seq = self.transcripts[transcript_id]
            else:
                # Try to find transcript via gene symbol mapping
                if hasattr(self, 'gene_to_transcript') and transcript_id in self.gene_to_transcript:
                    # Use the first transcript for this gene
                    transcript_list = self.gene_to_transcript[transcript_id]
                    if transcript_list:
                        transcript_seq = self.transcripts[transcript_list[0]]
                        logger.debug(f"Using transcript {transcript_list[0]} for gene {transcript_id}")
                    else:
                        logger.warning(f"No transcripts found for gene {transcript_id}")
                        return None
                else:
                    logger.warning(f"Transcript/gene {transcript_id} not found in transcript file")
                    return None
                
            # Calculate nucleotide position (approximate)
            # This assumes the CDS starts at position 0, which might not be accurate
            # A more precise implementation would require exact CDS coordinates
            nt_pos = (aa_pos - 1) * 3
            
            if nt_pos + 3 > len(transcript_seq):
                logger.warning(f"Position {nt_pos} out of range for transcript {transcript_id} (length {len(transcript_seq)})")
                return None
                
            # Determine reading frame by shifting one position (simulating a frameshift)
            # For a real implementation, you'd need to know the exact deletion that caused the frameshift
            frameshift_seq = transcript_seq[nt_pos + 1:]
            
            # Translate in the new reading frame until stop codon or reaching term_len amino acids
            translated = []
            for i in range(0, len(frameshift_seq) - 2, 3):
                codon = frameshift_seq[i:i+3]
                if len(codon) < 3:
                    break
                    
                try:
                    aa = str(Seq(codon).translate())
                    if aa == "*":  # Stop codon
                        break
                    translated.append(aa)
                    if len(translated) >= term_len:
                        break
                except Exception as e:
                    logger.warning(f"Error translating codon {codon}: {str(e)}")
                    translated.append("X")
                    
            result = ''.join(translated) if translated else None
            logger.debug(f"Translated frameshift for {transcript_id}: {result}")
            return result
            
        except Exception as e:
            logger.warning(f"Error translating frameshift for {transcript_id} at position {aa_pos}: {str(e)}")
            import traceback
            logger.debug(traceback.format_exc())
            return None
    
    def get_del_pos(self, protein):
        """
        Get position for frameshift deletion mutations
        
        Args:
            protein: Protein change notation (e.g., 'Q997Vfs*95')
            
        Returns:
            Position (integer)
        """
        try:
            # Handle full protein format (p.Q997Vfs*95)
            if protein.startswith('p.'):
                protein = protein[2:]  # Remove 'p.' prefix
                
            # Try with termination symbol
            match = re.match(r'([A-Z])(\d+)([A-Z])fs\*(\d+)', protein)
            if not match:
                # Try without termination symbol
                match = re.match(r'([A-Z])(\d+)([A-Z])fs', protein)
                
            if match:
                return int(match.group(2))
            else:
                logger.warning(f"Invalid frameshift format for position extraction: {protein}")
                return 0
        except Exception as e:
            logger.warning(f"Error extracting position for frameshift deletion {protein}: {str(e)}")
            return 0
    
    def get_aac_del(self, protein, alt_seq):
        """
        Get amino acid changes for frameshift deletions
        
        Args:
            protein: Protein change notation (e.g., 'Q997Vfs*95')
            alt_seq: Altered sequence
            
        Returns:
            Amino acid change string
        """
        try:
            # Handle full protein format (p.Q997Vfs*95)
            if protein.startswith('p.'):
                protein = protein[2:]  # Remove 'p.' prefix
                
            # Try with termination symbol
            match = re.match(r'([A-Z])(\d+)([A-Z])fs\*(\d+)', protein)
            if not match:
                # Try without termination symbol
                match = re.match(r'([A-Z])(\d+)([A-Z])fs', protein)
                
            if match:
                orig_aa = match.group(1)
                new_aa = match.group(3)
                return f"{orig_aa}>{new_aa}fs"
            else:
                # Fall back to first character if pattern doesn't match
                return protein[0] + ">fs"
        except Exception as e:
            logger.warning(f"Error extracting AA change for frameshift deletion {protein}: {str(e)}")
            return protein[0] + ">fs" if len(protein) > 0 else "X>fs"
    
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
    
    def get_junction_peptides(self, original_seq, alt_seq, del_start, del_end, min_len, max_len):
        """
        Generate peptides that span the junction of a frameshift
        
        Args:
            original_seq: Original protein sequence
            alt_seq: Altered protein sequence with frameshift
            del_start: Position of frameshift (1-based)
            del_end: Same as del_start for frameshift (1-based)
            min_len: Minimum peptide length
            max_len: Maximum peptide length
            
        Returns:
            List of tuples (title_suffix, peptide, is_novel)
        """
        result_peptides = []
        
        # Calculate the frameshift position in altered sequence (del_start - 1)
        junction_pos = del_start - 1
        
        # For each peptide length
        for peptide_len in range(min_len, max_len + 1):
            # Calculate how many amino acids to include before and after junction
            # At least 1 amino acid on each side of the junction is required
            for amino_before in range(1, peptide_len):
                amino_after = peptide_len - amino_before
                
                # Ensure enough sequence available
                if junction_pos - amino_before + 1 < 0:
                    continue
                if junction_pos + amino_after > len(alt_seq):
                    continue
                
                # Extract peptide spanning the junction
                junction_peptide = alt_seq[junction_pos - amino_before + 1 : junction_pos + amino_after]
                
                # Skip if peptide is too short
                if len(junction_peptide) < min_len:
                    continue
                
                # Determine if this peptide is novel (not found in original sequence)
                is_novel = self.is_peptide_novel(junction_peptide, original_seq)
                
                # Create descriptive title suffix showing:
                # 1. The position range in the altered sequence
                # 2. That it spans the junction
                # 3. Whether it's novel
                # Format: "pos1-pos2|junction|novel" or "pos1-pos2|junction"
                start_pos = junction_pos - amino_before + 2  # +2 for 1-based position and make it inclusive
                end_pos = junction_pos + amino_after
                title_suffix = f"{start_pos}-{end_pos}|junction"
                if is_novel:
                    title_suffix += "|novel"
                
                result_peptides.append((title_suffix, junction_peptide, is_novel))
        
        return result_peptides
    
    def save_var_pro(self):
        """Save all variant sequences to a FASTA file"""
        logger.info("Saving variant sequences to FASTA file")
        alt_df = pd.DataFrame()
        mut_type = 'fs_del'
        found_mutations = False
        
        # Load frameshift deletion mutation data
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
            fasta_path = os.path.join(self.out_file, 'FSDEL_Varsequence.fasta')
            var_pro_df.to_csv(fasta_path, index=False, header=None, sep='\t')
            logger.info(f"Saved {len(var_pro_df)//2} variant sequences to {fasta_path}")
            
            # Also save a summary file with mutation counts per gene
            gene_counts = alt_df['gene'].value_counts().reset_index()
            gene_counts.columns = ['Gene', 'MutationCount']
            summary_path = os.path.join(self.out_file, 'FSDEL_gene_summary.csv')
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
        # Generate junction peptides (8-11 amino acids)
        logger.info(f"Generating junction peptides (8-11 aa) for {tag}")
        mhc_df1 = self.create_junction_peps(alt_df, 8, 11)
        self.save_fasta(mhc_df1, f"{tag}_junction", 1)
        self.save_csv(mhc_df1, f"{tag}_junction", 1)
        
        # Generate longer junction peptides (15-30 amino acids)
        logger.info(f"Generating longer junction peptides (15-30 aa) for {tag}")
        mhc_df2 = self.create_junction_peps(alt_df, 15, 30)
        self.save_fasta(mhc_df2, f"{tag}_junction", 2)
        self.save_csv(mhc_df2, f"{tag}_junction", 2)
    
    def create_junction_peps(self, alt_df, min_len, max_len):
        """
        Create peptides specifically spanning frameshift junctions
        
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
                alt_seq = row['alt_seq']
                del_start = row['del_start']
                del_end = row['del_end']
                
                # Get junction peptides
                junction_peps = self.get_junction_peptides(
                    original_seq, alt_seq, del_start, del_end, min_len, max_len)
                
                for title_suffix, peptide, is_novel in junction_peps:
                    # Create full title
                    title = f">{gene}|p.{protein_change}|{title_suffix}"
                    result_peptides.append([title, peptide, is_novel])
                    total_count += 1
                    if is_novel:
                        novel_count += 1
                
            except Exception as e:
                logger.warning(f"Error generating junction peptides for {row.get('gene', 'unknown')}: {str(e)}")
                continue
        
        # Create DataFrame from results
        if result_peptides:
            result_df = pd.DataFrame(result_peptides, columns=['title', 'peptide', 'is_novel'])
            
            # Log statistics
            logger.info(f"Generated {len(result_df)} junction peptides (length {min_len}-{max_len})")
            logger.info(f"Found {novel_count}/{total_count} novel peptides that are not present in original sequences")
            
            # For MHC prediction, we want novel peptides only
            novel_df = result_df[result_df['is_novel']]
            logger.info(f"Using {len(novel_df)} novel peptides for MHC prediction")
            
            return novel_df[['title', 'peptide']]
        else:
            logger.warning(f"No junction peptides generated for length range {min_len}-{max_len}")
            return pd.DataFrame(columns=['title', 'peptide'])
    
    def create_peps(self, alt_df, min_len, max_len):
        """
        DEPRECATED: Original method for creating peptides.
        This is kept for backward compatibility and will be removed in future versions.
        
        Please use create_junction_peps instead.
        """
        logger.warning("The create_peps method is deprecated. Please use create_junction_peps instead.")
        return self.create_junction_peps(alt_df, min_len, max_len)
    
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
    """Parse arguments and run the frameshift deletion peptide generator"""
    parser = argparse.ArgumentParser(description='Generate variant peptides for frameshift deletions from ANNOVAR annotations')
    parser.add_argument('-i', '--input', required=True, help='Input ANNOVAR annotation file')
    parser.add_argument('-r', '--reference', required=True, help='Reference protein sequence database in FASTA format')
    parser.add_argument('-o', '--output', required=True, help='Output directory for results')
    parser.add_argument('-p', '--peptides', action='store_true', help='Generate peptide files for MHC prediction')
    parser.add_argument('-v', '--verbose', action='store_true', help='Enable verbose logging')
    parser.add_argument('-t', '--transcript', help='Transcript sequences in FASTA format (for real translation)')
    parser.add_argument('--real-translation', action='store_true', 
                       help='Use real translation for frameshift instead of X placeholders')
    parser.add_argument('--debug', action='store_true', help='Enable debug logging')
    
    args = parser.parse_args()
    
    # Set logging level based on verbosity
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    elif args.verbose:
        logging.getLogger().setLevel(logging.INFO)
        
    logger.info(f"Starting frameshift deletion variant peptide generation")
    logger.info(f"Input file: {args.input}")
    logger.info(f"Reference file: {args.reference}")
    logger.info(f"Output directory: {args.output}")
    logger.info(f"Generate peptides: {args.peptides}")
    
    if args.real_translation and not args.transcript:
        logger.warning("Real translation requested but no transcript file provided. Falling back to X placeholders.")
        args.real_translation = False
    
    if args.real_translation:
        logger.info(f"Using real translation with transcript file: {args.transcript}")
    
    generator = DelPeptideGenerator(
        input_file=args.input,
        reference_file=args.reference,
        output_dir=args.output,
        transcript_file=args.transcript,
        generate_peptides=args.peptides,
        use_real_translation=args.real_translation
    )
    
    generator.run()

if __name__ == "__main__":
    main() 