'''
Created on May 14, 2026

@author: pleyte
'''
import argparse
import logging
import pandas as pd

class CompareResults(object):
    '''
    classdocs
    '''
    def __init__(self):
        '''
        Constructor
        '''
        self._logger = logging.getLogger(__name__)
        self._input_dfs = []
        self._df = None

    def add_datasource(self, label, path):
        '''
        Read nomenclature and add it to our dataframe
        '''
        df = pd.read_csv(path)
        df['source'] = label
        self._input_dfs.append(df)
        self._logger.info(f"Read {label} dataframe with {len(df)} rows")
    
    def create_comparison(self):
        self._df = pd.concat(self._input_dfs, ignore_index=True)
    
    def _reconstruct_into_blocks(self):
        """Implements the 'Sheet per Variant' logic."""
        identity_cols = ['chromosome', 'position', 'reference', 'alt', 'cdna_transcript']
        all_sources = self._df['source'].unique()
        
        reconstructed_blocks = []
        
        # Group by the variant identity (like looking at one 'sheet' at a time)
        for identity, variant_df in self._df.groupby(identity_cols):
            # Create a 3-row template for this variant
            template = pd.DataFrame(
                [list(identity) + [s] for s in all_sources], 
                columns=identity_cols + ['source']
            )
            # Merge actual data into the template
            block = pd.merge(template, variant_df, on=identity_cols + ['source'], how='left')
            reconstructed_blocks.append(block)
        
        return pd.concat(reconstructed_blocks, ignore_index=True)
    
    def write(self, output_file):
        # 1. First, transform the flat data into the interleaved blocks
        final_df = self._reconstruct_into_blocks()
        
        # 2. Identify the groups for coloring
        identity_cols = ['chromosome', 'position', 'reference', 'alt', 'cdna_transcript']
        final_df['group_id'] = final_df.groupby(identity_cols).ngroup()

        # 3. Create the styling function
        def style_alternating(row):
            # Light grey for even groups, white for odd
            color = 'background-color: #F2F2F2' if row.group_id % 2 == 0 else ''
            return [color] * len(row)

        # 4. Apply style and export (dropping the helper group_id)
        styled_df = final_df.style.apply(style_alternating, axis=1)
        
        # Remove group_id from the final export
        cols_to_show = [c for c in final_df.columns if c != 'group_id']
        styled_df.to_excel(output_file, index=False, columns=cols_to_show)
        
def _parse_args():
    parser = argparse.ArgumentParser(description='Read results from multiple tools and produce comparison spreadsheet')
    parser.add_argument("--version", action="version", version="0.0.1")    
    parser.add_argument("--input", help="Pair of label and file path (e.g., --input tfx file.csv)", dest="inputs", nargs=2, metavar=('LABEL', 'PATH'), action="append", required=True)
    parser.add_argument("--out", help="Comparison file (xlsx)", dest="output", required=True)
    args = parser.parse_args()
    return args

def main():
    from nomenclature_comparison.util.logger import setup_logging
    setup_logging(level=logging.INFO)

    args = _parse_args()
    cr = CompareResults()
    
    for label, filepath in args.inputs:
        cr.add_datasource(label, filepath)
    
    cr.create_comparison()
    cr.write(args.output)
    
if __name__ == '__main__':
    main()        