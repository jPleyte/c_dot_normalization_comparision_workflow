'''
Created on May 14, 2026

@author: pleyte
'''
from nomenclature_comparison.variant_transcript import VariantTranscript
import csv
import logging

def write(out_filename: str, variants: list[VariantTranscript]):
    """
    Write a list of VariantTranscript to csv 
    """
    key_headers = ['chromosome', 'position', 'reference', 'alt', 'cdna_transcript' ] 
    nomenclature_headers = ['c_dot', 'g_dot', 'gene', 'p_dot1', 'p_dot3', 'protein_transcript']
    all_headers = key_headers + nomenclature_headers 
    
    rows = 0
    with open(out_filename, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(all_headers)
        
        for v in variants:
            rows += 1 
            row = [v.chromosome,
                   v.position,
                   v.reference,
                   v.alt,
                   v.cdna_transcript,
                   v.c_dot,
                   v.g_dot,
                   v.gene,
                   v.p_dot1,
                   v.p_dot3,
                   v.protein_transcript]        
                
            writer.writerow(row)
        
    logging.info(f"Wrote {rows} to {out_filename}")