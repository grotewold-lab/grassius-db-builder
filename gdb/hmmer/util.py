# this file contains utility functions to support hmmer

from .hmmscan_result import HmmscanResult
import pandas as pd

def build_minified_hmm( hmm_path, domain_subset, output_path ):
    """
    Create a new hmm file containing only the necessary domains,
    in order to run hmmscan efficiently
    
    
    Arguments:
    ----------
    hmm_path -- (str) the path to the existing hmm file containing hidden-markov models
    domain_subset -- (list of str) a subset of domain names present in the existing hmm file
    output_path -- (str) the path for the resulting file which will be created or replaced
    """
    
    return None
    
    

def read_hmmscan_output(path):
    """
    Parse a file containing output from hmmscan
    return an instance of HmmscanResult
    """
    
    # prepare dataframe to contain data
    colnames = [
        "target name","accession","query name","accession",
        "E-value","score","bias","E-value","score","bias",
        "exp","reg","clu","ov","env","dom","rep","inc"
    ]
    ncols = len(colnames)
    df = pd.DataFrame(columns=colnames)
    
    # start reading the output file
    with open(path) as fin:
        
        # skip 3 lines
        [fin.readline() for i in range(3)]
        
        # parse data
        i = 0
        while True:
            line = fin.readline()
            if line.startswith("#"):
                break
                
            df.loc[i,:] = line.split()[:ncols]
            i += 1
    
        # extract footer
        footer = ""
        while True:
            line = fin.readline()
            if not line:
                break
            footer += line
                
    
    return HmmscanResult(df, footer)