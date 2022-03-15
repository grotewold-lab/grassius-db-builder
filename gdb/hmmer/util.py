# this file contains utility functions to support hmmer

from .hmmscan_result import HmmscanResult
import pandas as pd


def read_family_criteria(path):
    """
    Load an excel sheet containing rules for assigning families to protein sequences
    
    The following columns will be considered in an automated pipeline:
        "GRASSIUS" - the family name
        "Required" - accessions that must be present (one or more accessions)
        "Forbidden" - accessions that must NOT be present (zero or more accessions)
    
    Formatting of "Required" and "Forbidden" columns is based on iTAK
    
    return a dataframe
    """
    return pd.read_excel(path).fillna('').iloc[:,:12]

    
def get_relevant_accessions(family_criteria):
    """
    Get the minimum set of accessions that should be considered 
    in order to assign families to protein sequences
    
    Arguments:
    ----------
    family_criteria -- (DataFrame) criteria returned by read_family_criteria
    """
    
    df = family_criteria
    all_acc = set()
    for row in df.index:
        for col in ["Required","Forbidden"]:
            all_acc.update( [s[:-2] for s in df.loc[row,col].split(":")] )
    all_acc.discard('')
    
    return all_acc

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