# this file contains utilities that are not specific to any of the
# submodules like "blast" or "fasta"

import pandas as pd

def load_gene_annotations(gff3_path):
    """
    load a subset of annotations from a gff3 file
    returns a dataframe that may be passed to annotate_blast_result()
    """
    
    df = pd.read_table(gff3_path, skiprows=6, header=None)
    df.columns = ["chrom","1","type","start","end","5","6","7","identifiers"]
    df = df[df['type'] == 'gene']
    
    return df
    