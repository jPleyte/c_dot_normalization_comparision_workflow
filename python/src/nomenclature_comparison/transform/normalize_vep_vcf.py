'''
Created on May 14, 2026

@author: pleyte
'''
import logging
import argparse
from nomenclature_comparison.util.pdot import PDot
import pandas as pd
import re
from nomenclature_comparison.util import chromosome_map
import urllib
from nomenclature_comparison.variant_transcript import VariantTranscript
from collections import Counter
from nomenclature_comparison.io import normalized_writer
import pysam

class NormalizeVepVcf(object):
    '''
    Transform VEP results into csv with common headings
    '''
    def __init__(self):
        '''
        Constructor
        '''
        self._logger = logging.getLogger(f"nomenclature_comparison.{__name__}")
        self._p_dot_mapper = PDot()
        self._variant_transcripts = []

    def read_vep_file(self, vep_vcf_file):
        """
        Read the VEP VCF file.   
        """
        vcf_in = pysam.VariantFile(vep_vcf_file)
        csq_header = vcf_in.header.info['CSQ'].description.split("Format: ")[-1].strip('"').split('|')
        
        # Create a lookup map (e.g., {'Feature': 2, 'HGVSc': 3})
        field_to_idx = {field: idx for idx, field in enumerate(csq_header)}
        
        parsed_records = []
        
        for rec in vcf_in:        
            # Keep your pristine VCF/CSV identity columns!
            chrom = rec.contig
            pos = rec.pos
            ref = rec.ref
            alt = rec.alleles[1] # Assumes mono-allelic rows
            
            # 4. Unpack the CSQ annotations
            if 'CSQ' in rec.info:
                for transcript_annotation in rec.info['CSQ']:
                    csq_parts = transcript_annotation.split('|')
                    
                    # Extract fields using our index map
                    # VEP uses 'Feature' to hold the Transcript ID (e.g., NM_xxx)
                    cdna_transcript = csq_parts[field_to_idx.get('Feature')] if 'Feature' in field_to_idx else None
                    
                    hgvs_g = csq_parts[field_to_idx.get('HGVSg')] if 'HGVSg' in field_to_idx else None
                    
                    ccds_accession = csq_parts[field_to_idx.get('CCDS')] if 'CCDS' in field_to_idx else None 
                    if ccds_accession:
                        cdna_transcript = ccds_accession
                        
                    # HGVSc holds the c. variant (e.g., ENST00000xxx:c.1A>G or NM_xxx:c.1A>G)
                    hgvsc_full = csq_parts[field_to_idx.get('HGVSc')] if 'HGVSc' in field_to_idx else None
                    
                    # HGVSp holds the p. variant (e.g., NP_xxx:p.Arg22Ter)
                    hgvsp_full = csq_parts[field_to_idx.get('HGVSp')] if 'HGVSp' in field_to_idx else None
                    
                    # ENSP or NP protein ID is tucked inside the 'ENSP' field or parseable from HGVSp
                    protein_transcript = csq_parts[field_to_idx.get('ENSP')] if 'ENSP' in field_to_idx else None
                    
                    gene_symbol = csq_parts[field_to_idx.get('SYMBOL')] if 'SYMBOL' in field_to_idx else None
                    
                    # If you ran VEP with --hgvs, HGVSc/p includes the transcript prefix.
                    # Let's clean them so you just get the "c." and "p." segments:
                    c_dot = hgvsc_full.split(':')[-1] if ':' in hgvsc_full else hgvsc_full
                    p_dot = hgvsp_full.split(':')[-1] if ':' in hgvsp_full else hgvsp_full
                                            
                    vt = self._get_variant_transcript(chrom, pos, ref, alt, cdna_transcript, c_dot, p_dot, protein_transcript, hgvs_g, gene_symbol)
                    self._variant_transcripts.append(vt)
                    
                

    def _get_variant_transcript(self, chrom, pos, ref, alt, cdna_transcript, c_dot, p_dot3, protein_transcript, hgvs_g, gene_symbol):
        '''
        '''
        vt = VariantTranscript(chrom, pos, ref, alt, cdna_transcript)
        vt.c_dot = c_dot
        vt.g_dot = hgvs_g
        vt.gene = gene_symbol
        vt.p_dot3 = p_dot3
        
        if p_dot3:
            vt.p_dot1 = self._p_dot_mapper.get_p_dot1(protein_transcript, p_dot3)

        vt.protein_transcript = protein_transcript
        return vt
        
    def _get_variant_values(self, row: pd.Series):
        """
        Return the variant and transcript parts. Pulls values from multiple fields and compares them to make sure they all match
        """
        if row['SOURCE'] == 'RefSeq':            
            transcript = row['Feature']            
        elif row['SOURCE'] == 'Ensembl' and row['CCDS'] != '-':
            transcript = row['CCDS']
        elif row['SOURCE'] == 'Ensembl' and row['Feature'].startswith('ENST'):
            transcript = row['Feature']
        else:
            raise ValueError(f"Unknown feature type: Source={row['SOURCE']}, Feature={row['Feature']}, CCDS={row['CCDS']}")
        
        chromosome, position, reference, alt = self._get_variant_values_from_gdot(row['#Uploaded_variation'])
        return chromosome, position, reference, alt, transcript
    
    def _get_variant_values_from_gdot(self, uploaded_variation: str):
        """
        Parse the variant values out of a Uploaded_variation field (eg 1_100908484_T/C) 
        """
        pattern = r"^(.+)_(\d+)_([ACGTN*.-]+)\/([ACGTN*.-]+)$"
        match = re.match(pattern, uploaded_variation)
                
        chrom, pos, ref, alt = match.groups()
        assert chrom and pos and ref and alt, "Variant parts could not be identified: "
        return chrom, pos, ref, alt
        
    def _get_c_dot(self, row):
        """
        Parse c. from the HGVSc field
        """
        if row['HGVSc'] == '-':
            return None

        transcript, c_dot = row['HGVSc'].split(':')
        
        if transcript.startswith('NM') and transcript != row['Feature']:
            raise ValueError(f"c. transcript does not match Feature: {transcript} != {row['Feature']}")
        elif not c_dot.startswith('c.'):
            raise ValueError(f"c. does not start with c.: {row['HGVSc']}")
        
        return c_dot
        
    def _get_g_dot(self, row):
        """
        Use HGVSg to return g.
        They use ncbi numbers for chromosome rather than refseq accessions so this function converts chromosome to refseq
        """
        chromosome, g_dot = row['HGVSg'].split(':')
        refseq_chromosome = chromosome_map.get_refseq(chromosome)
        return f"{refseq_chromosome}:{g_dot}"
     
    def _get_genomic_region_type(self, row) -> str:
        """
        Map the Consequence field to a one of our four genomic region types. 
        See https://useast.ensembl.org/info/genome/variation/prediction/predicted_data.html
        """
        cons = row['Consequence'].lower()
        biotype = row['BIOTYPE'].lower()
        
        if 'upstream' in cons or 'downstream' in cons or 'intergenic' in cons:
            return 'intergenic'
        
        if "splice_acceptor_variant" in cons or "splice_donor_variant" in cons:
            return "splicing"
    
        if biotype == 'protein_coding':        
            exonic_terms = [ 'missense', 'synonymous', 'stop_gained', 'stop_lost', 
                            'frameshift', 'inframe_insertion', 'inframe_deletion', 
                            'coding_sequence_variant', 'protein_altering_variant', 
                            'start_lost', 'stop_retained_variant' ]
            
            if any(term in cons for term in exonic_terms):
                return 'exon'
        
            # UTR terms
            if 'utr_variant' in cons:
                return 'utr'
                
            # Intronic terms
            if 'intron_variant' in cons:
                # Special check: Essential Splice sites are intronic but critical
                if 'splice_donor' in cons or 'splice_acceptor' in cons:
                    return 'instron/splicing'
                return 'intron'
        
        
        # Non-coding or "Decay" transcripts
        # This handles biotypes like 'nonsense_mediated_decay' or 'antisense'
        if any(x in biotype for x in ['decay', 'antisense', 'rna', 'pseudogene', 'retained']):
            return 'non-coding'
    
        # Non-coding Exons (like in lincRNAs)
        if 'non_coding_transcript_exon_variant' in cons:
            return 'non-coding'
        
        raise ValueError(f"Unable to determine region type for {row['#Uploaded_variation']}: cons={cons}, biotype={biotype}")
    
    def _get_protein_variant_type(self, row):
        """
        Map VEP feilds to a protein variant type        
        """
        cons = row['Consequence'].lower()
        biotype = row['BIOTYPE'].lower()
    
        # 1. Essential Splice Loss (High Priority)
        if 'splice_donor' in cons or 'splice_acceptor' in cons:
            return 'Splice junction loss'
        
        # 2. Start/Stop Changes
        if 'start_lost' in cons: 
            return 'Start loss'
        if 'stop_gained' in cons: 
            return 'Stop gain'
        if 'stop_lost' in cons: 
            return 'Stop loss'
        
        # 3. Coding Changes
        if 'frameshift' in cons: 
            return 'Frameshift'
        
        # ITD Logic: Usually an inframe insertion where the sequence repeats
        # This is a placeholder; real ITD detection often needs the VCF ALT allele string
        if 'inframe_insertion' in cons:
            # If your data has a flag for duplication, use it here
            return 'In-frame'
        
        if 'inframe_deletion' in cons: 
            return 'In-frame'
        
        if 'missense' in cons:
            return 'Missense'
        
        if 'synonymous' in cons: 
            return 'Synonymous'
        
        if cons == 'stop_retained_variant':
            return 'Synonymous'
        
        # 4. MNV (Multi-nucleotide variant)
        # VEP often labels these as 'protein_altering_variant'
        if 'protein_altering' in cons: 
            return 'Multi nucleotide variant'
        
        # This isn't one of the CGD types but what else  to do with downstream_gene_variant?
        if 'downstream_gene_variant' in cons or 'intergenic' in cons:
            return 'Flanking'
    
        # 5. Promoter vs Flanking
        if 'upstream' in cons:
            # Standard convention: Promoter is within 2kb of the start
            # VEP provides 'DISTANCE' if you enable it
            return 'Flanking'
    
        # 6. Intron
        if 'intron' in cons: 
            return 'Intron'
        
        outside_terms = ['utr_variant', 'downstream', 'intergenic', 'non_coding_transcript_exon']
        if any(term in cons for term in outside_terms):
            return 'Flanking'
        
        # I don't know what to do with this one
        if cons == 'coding_sequence_variant' and biotype == 'protein_coding':
            return None     
        
        raise ValueError(f"Unable to determine protein variant type for {row['#Uploaded_variation']}: cons={cons}, biotype={biotype}")
        
        
    def _get_protein_change(self, row):
        """
        Return protein transcript and three letter p.
        """
        if row['HGVSp'] == '-':
            return None, None
        
        protein_transcript, p_dot3_raw = row['HGVSp'].split(':')
        
        # VEP uses "%3D" instead of "="
        p_dot3 = urllib.parse.unquote(p_dot3_raw)
        
        if row['SOURCE'] == 'RefSeq' and not protein_transcript.startswith("NP"):
            raise ValueError(f"Protein transcript does not start with NP: {row['HGVSp']}")
        elif not p_dot3.startswith("p."):
            raise ValueError(f"p. does not start with p.: {row['HGVSp']}")
        
        return protein_transcript, p_dot3
    
    def _get_exon(self, row):
        """
        Return the exon. VEP stores them in a string like '3/9'
        """
        if row['EXON'] == '-' and row['INTRON'] == '-':
            return None
        elif row['INTRON'] == '-':
            value = row['EXON']
        elif row['EXON'] == '-':
            value = row['INTRON']
        elif '/' in row['EXON'] and '/' in row['INTRON']:
            value = row['EXON']
        else: 
            raise ValueError(f"Unknown situation with exon/intron: {row['EXON']} {row['INTRON']}")
        
        return value.split('/')[0]
    
    def get_variant_transcripts(self) -> list[VariantTranscript]:
        """
        Parse the VEP dataframe        
        """
        transcript_counter = Counter()
        variant_transcripts = []
        
        for index, row in self._vep_df.iterrows():
            chromosome, position, reference, alt, transcript = self._get_variant_values(row)

            if transcript.startswith('NR'):
                transcript_counter['skipped_NR_transcript']
                continue
            elif transcript.startswith('ENST'):
                transcript_counter['skipped_ENST_transcript']
                continue
            
            vt = VariantTranscript(chromosome, position, reference, alt, transcript)
            vt.c_dot = self._get_c_dot(row)
            vt.g_dot = self._get_g_dot(row)
            
            protein_transcript, p_dot3 = self._get_protein_change(row)
            vt.protein_transcript = protein_transcript
            vt.p_dot3 = p_dot3            
            if vt.p_dot3:
                vt.p_dot1 = self._p_dot_mapper.get_p_dot1(protein_transcript, p_dot3)
            
            vt.gene = row['SYMBOL'] if row['SYMBOL'] != '-' else None
            vt.strand = row['STRAND']
            
            transcript_counter['valid'] += 1
            variant_transcripts.append(vt)
        
        self._logger.info(f"Vep results: {transcript_counter}")
        return variant_transcripts
            
    def write(self, output_filename: str, keep_enst=True):
        """
        """
        accession_counter = Counter()
        variant_transcripts = [] 
        for x in self._variant_transcripts:
            if x.cdna_transcript.startswith('NM'):                        
                accession_counter['NM'] += 1
                variant_transcripts.append(x)
            elif x.cdna_transcript.startswith('CCDS'):
                accession_counter['CCDS'] += 1
                variant_transcripts.append(x)
            elif x.cdna_transcript.startswith('ENST'):
                accession_counter['ENST'] += 1
                if keep_enst:
                    variant_transcripts.append(x)
            elif x.cdna_transcript.startswith('NR'):
                # Never keep NR
                accession_counter['NR'] += 1                
            else:
                print(f"Other: {x.cdna_transcript}")
                accession_counter['other'] += 1
                
        normalized_writer.write(output_filename, variant_transcripts)
        self._logger.info(f"Transcript types: {accession_counter}")
        self._logger.info(f"Wrote {len(self._variant_transcripts)} variant transcripts to {output_filename}")            


def _parse_args():
    parser = argparse.ArgumentParser(description='Read transcript nomenclature from vep output and write to csv')
    parser.add_argument("--version", action="version", version="0.0.1")    
    parser.add_argument("--in", help="File generated by VEP (tsv)", dest="input", required=True)
    parser.add_argument("--out", help="Normalized VEP output file (csv)", dest="output", required=True)
    args = parser.parse_args()
    return args    

def main():
    from nomenclature_comparison.util.logger import setup_logging
    setup_logging(level=logging.INFO)

    args = _parse_args()

    norm_vep = NormalizeVepVcf()
    
    norm_vep.read_vep_file(args.input)    
    norm_vep.write(args.output, keep_enst=False)
    
if __name__ == '__main__':
    main()