'''
Convert csv to tab delimited Annovar avinput file
Created on Mar 4, 2025

@author: pleyte
'''
import argparse
import csv
import logging
from enum import Enum

VERSION = '0.0.1'

logger = logging.getLogger(__name__)

class Field(str, Enum):
    # Update these values if the CSV headers change
    CHROMOSOME = 'chromosome'
    POSITION   = 'position_start'  
    REFERENCE  = 'reference_base'
    ALT  = 'variant_base'
    ID = 'genomic_variant'
    
    def __str__(self):
        return str(self.value)
    
class CsvToAvinput(object):
    '''
    Convert a list of variants in csv format to an Annovar avinput file 
    '''
    def __init__(self):
        '''
        Constructor
        '''
        self._logger = logging.getLogger(__name__)                

    def read(self, in_stream):
        """
        Read the csv and return list of variants 
        """         
        variants = {}
        
        for row in csv.DictReader(in_stream):
            chromosome = row[Field.CHROMOSOME]
            position = row[Field.POSITION]
            reference = row[Field.REFERENCE]
            alt = row[Field.ALT]
            key = f"{chromosome}-{position}-{reference}-{alt}"
            if key not in variants:
                variants[key] = row
        
        return variants.values()

    def write(self, variants: list, out_stream):
        """
        Write list of variants to tab delimited avinput file 
        """
        tsv_writer = csv.writer(out_stream, delimiter='\t')
        for x in variants:
            position_end = int(x[Field.POSITION]) + len(x[Field.REFERENCE]) - 1      
            
            external_id = '' if not x.get(Field.ID) else f"id:{x[Field.ID]}"
            tsv_writer.writerow([self._get_chromosome(x[Field.CHROMOSOME]),
                                 x[Field.POSITION],
                                 position_end,
                                 x[Field.REFERENCE],
                                 x[Field.ALT],
                                 external_id])

    def _get_chromosome(self, chromosome):
        return chromosome.replace('chr', '')


def _parse_args():
    parser = argparse.ArgumentParser(description='Read variants from a csv formatted file and write out an Annovar input file.')

    parser.add_argument('-i', '--in',
                dest="input",
                help='csv list of variants',
                type=argparse.FileType('r'),
                required=True)

    parser.add_argument('-o', '--out',
                dest="output",
                help='avinput file',
                type=argparse.FileType('w'),
                required=True)

    parser.add_argument('--version', action='version', version=VERSION)

    return parser.parse_args()


def main():
    from nomenclature_comparison.util.logger import setup_logging
    setup_logging(level=logging.INFO)
    
    args = _parse_args()

    c2a = CsvToAvinput()
    variants = c2a.read(args.input)

    c2a.write(variants, args.output)

    logger.info(f"Wrote {len(variants)} variants to {args.output.name}")


if __name__ == '__main__':
    main()
