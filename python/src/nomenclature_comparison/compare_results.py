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
        
        reconstructed_blocks = []
        
        # Group by the variant identity
        for identity, variant_df in self._df.groupby(identity_cols):
            # 1. Extract the 'c_dot' values for our two specific TFX sources
            # We use .loc to find the rows matching the source, and safely grab the first value if it exists
            cdot_branch_one = variant_df.loc[variant_df['source'] == 'tfx_bio831_two_gdots', 'c_dot'].values
            cdot_branch_two = variant_df.loc[variant_df['source'] == 'tfx_hgvs_normalizer', 'c_dot'].values
            # 2. Determine agreement (1 or 0)
            # Default to 0 unless BOTH branches have an entry and their c_dot strings match exactly
            # 2. Determine agreement (1 or 0)
            if len(cdot_branch_one) > 0 and len(cdot_branch_two) > 0:
                agreement = 1 if cdot_branch_one[0] == cdot_branch_two[0] else 0
            else:
                agreement = 0

            # 3. Add the new field to every row currently inside this block
            # (Using a copy prevents Pandas SettingWithCopy warnings)
            sorted_variant_df = variant_df.sort_values(by='source').copy()
            sorted_variant_df['tfx_dev_branch_agree'] = None
            
            # Then, we target only the rows matching our specific TFX sources and assign the 1 or 0
            tfx_rows_mask = sorted_variant_df['source'].isin(['tfx_bio831_two_gdots', 'tfx_hgvs_normalizer'])
            sorted_variant_df.loc[tfx_rows_mask, 'tfx_dev_branch_agree'] = agreement
            
            # Append the actual data rows that exist for this variant
            reconstructed_blocks.append(sorted_variant_df)
            
            # Add a single entirely blank row to act as a visual spacer 
            # between blocks in your final Excel sheet
            spacer_row = pd.DataFrame([{col: None for col in self._df.columns}])
            reconstructed_blocks.append(spacer_row)
        
        # Combine everything back together, omitting the last spacer row if you want it clean
        if reconstructed_blocks:
            return pd.concat(reconstructed_blocks[:-1], ignore_index=True)
        
        return pd.DataFrame(columns=self._df.columns)

    
    def write(self, output_file):
        from openpyxl.styles import Font

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
        # styled_df.to_excel(output_file, index=False, columns=cols_to_show)

        # 5. Use ExcelWriter with the openpyxl engine to combine Styler and custom formatting
        with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
            # Write the styled dataframe to the workbook
            styled_df.to_excel(writer, sheet_name='Comparisons', index=False, columns=cols_to_show)
            
            # Access openpyxl objects directly from the writer
            worksheet = writer.sheets['Comparisons']
            
            # --- FREEZE THE TOP ROW ---
            # 'A2' tells Excel to freeze everything above row 2 (leaving row 1 fixed)
            worksheet.freeze_panes = 'A2'
            
            # --- MAKE THE TOP ROW BOLD ---
            bold_font = Font(bold=True)
            
            # Loop through all header cells in the first row and apply the bold style
            for col_num in range(1, len(cols_to_show) + 1):
                cell = worksheet.cell(row=1, column=col_num)
                cell.font = bold_font

        
        
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