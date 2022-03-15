# this file contains utility functions to support hmmer

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
    
    return HmmscanResult()