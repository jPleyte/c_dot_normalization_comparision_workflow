'''
Convert csv to tab delimited SnpEff input file
Created on May 14, 2026

@author: pleyte
'''
import argparse
import csv
import pysam
from enum import Enum
import functools
import logging
from nomenclature_comparison.util import chromosome_map

VERSION = '0.0.1'

class Field(str, Enum):
    # Update these values if the CSV headers change
    CHROMOSOME = 'chromosome'
    POSITION   = 'position_start'  
    REFERENCE  = 'reference_base'
    ALT  = 'variant_base'
    ID = 'genomic_variant'
    
    def __str__(self):
        return str(self.value)

@functools.total_ordering
class SimpleVariant:
    # Official genomic order for hg19/hg38
    _CHROM_ORDER = {str(i): i for i in range(1, 23)}
    _CHROM_ORDER.update({f"chr{i}": i for i in range(1, 23)})
    _CHROM_ORDER.update({"X": 23, "chrX": 23, "Y": 24, "chrY": 24, "MT": 25, "M": 25, "chrM": 25})

    def __init__(self, chromosome, position, reference, alt):
        self.chromosome = str(chromosome).replace('chr', '')
        self.position = int(position)
        self.reference = str(reference).upper()
        self.alt = str(alt).upper()

    def _sort_key(self):
        """Internal helper to return a (Rank, Position) tuple."""
        # Use 999 as a fallback for unknown scaffolds/contigs
        rank = self._CHROM_ORDER.get(self.chromosome, 999)
        return (rank, self.position, self.reference, self.alt)

    def __hash__(self):
        # Hash the tuple of all fields to ensure uniqueness in a set
        return hash((self.chromosome, self.position, self.reference, self.alt))

    def __eq__(self, other):
        if not isinstance(other, SimpleVariant):
            return NotImplemented
        return (self.chromosome, self.position, self.reference, self.alt) == \
               (other.chromosome, other.position, other.reference, other.alt)

    def __lt__(self, other):
        if not isinstance(other, SimpleVariant):
            return NotImplemented
        return self._sort_key() < other._sort_key()

    def __repr__(self):
        return f"SimpleVariant({self.chromosome}:{self.position} {self.reference}>{self.alt})"
    
class CsvToVcf(object):
    '''
    classdocs
    '''
    def __init__(self):
        '''
        Constructor
        '''
        self._logger = logging.getLogger(f"nomenclature_comparison.{__name__}")
    
    def read(self, in_file_csv):
        """
        Read the csv and return list of variants 
        """         
        variant_set = set()
        
        with open(in_file_csv, mode='r', newline='', encoding='utf-8') as csvfile:
            for row in csv.DictReader(csvfile):
                v = SimpleVariant(row[Field.CHROMOSOME], row[Field.POSITION], row[Field.REFERENCE], row[Field.ALT])
                variant_set.add(v)
                
        return list(variant_set)
    
    def write(self, variants: list, out_file_vcf):
        """
        Write list of variants to VCF file 
        """
        header = pysam.VariantHeader()
        header.add_line("##fileformat=VCFv4.2")        
        header.add_line('##FILTER=<ID=PASS,Description="All filters passed">')
        
        for x in chromosome_map.refseq_to_ncbi.values():            
            header.contigs.add(x)

        with pysam.VariantFile(out_file_vcf, "w", header=header) as vcf_out:
            for x in variants:
                # pysam expects zero based position, and then adds one
                start_pos = x.position - 1
                rec = vcf_out.new_record(
                    contig=x.chromosome,
                    start=start_pos, 
                    stop=start_pos + len(x.reference),
                    alleles=(x.reference, x.alt),
                    id=".",
                    qual=None,
                    filter="PASS")

                vcf_out.write(rec)
        
        print("Done")
        self._logger.info(f"Wrote {len(variants)} variants to {out_file_vcf}")

def _parse_args():
    parser = argparse.ArgumentParser(description='Read variants from a csv formatted file and write out a VCF')

    parser.add_argument('--in',
                dest="input",
                help='csv list of variants',
                required=True)

    parser.add_argument('--out',
                dest="output",
                help='Vcf out file',
                required=True)

    parser.add_argument('--version', action='version', version=VERSION)

    return parser.parse_args()


def main():
    from nomenclature_comparison.util.logger import setup_logging
    setup_logging(level=logging.INFO)
    
    args = _parse_args()

    c2v = CsvToVcf()
    variants = c2v.read(args.input)
    
    # This sort is not effective
    sorted_variants = sorted(variants)
    
    c2v.write(sorted_variants, args.output)

if __name__ == '__main__':
    main()
