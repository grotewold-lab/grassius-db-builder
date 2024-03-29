# this file contains utility functions to support hmmer

# local imports
from .hmmscan_result import HmmscanResult

import pandas as pd

def get_acc_dict(hmmscan_result):
    """
    Summarize hmmscan results in the form of a dictionary
    
    Results should be filtered before using this function,
    because the resulting object will not contain scores.
    
    Return a dictionary where keys are transcript IDs, 
    and values are lists of matching accession names
    
    Arguments:
    ----------
    hmmscan_result -- an instance of HmmscanResult
    """
    
    df = hmmscan_result.data
    result = {}
    for row in df.index:        
        tid,acc = df.loc[row,["query name","accession"]]
        acc = df.loc[row,"accession"].split(".")[0]
        if tid not in result.keys():
            result[tid] = []
        result[tid].append( acc )
    return result


def concatenate_hmms( hmm_paths, output_path ):
    """
    Create a new hmm file containing the contents of two existing hmm files
    
    
    Arguments:
    ----------
    hmm_paths -- (list of str) paths to existing hmm files
    output_path -- (str) the path for the resulting file which will be created or replaced
    """
    
    with open(output_path,"w") as fout:
        for input_path in hmm_paths:
            with open(input_path) as fin:
                while True:
                    line = fin.readline()
                    if not line:
                        break
                    fout.write(line)
                
    

def get_accessions( hmm_path ):
    """
    Get a list of accession names present in the given hmm file
    """
    result = []
    
    with open(hmm_path) as fin:
        while True:
            line = fin.readline()
            if not line:
                break
            if line.startswith("ACC"):
                acc = line.strip().split()[1].split(".")[0]
                result.append(acc)

    return result


def build_minified_hmm( hmm_path, domain_subset, output_path ):
    """
    Create a new hmm file containing only the necessary domains,
    in order to be used as input for hmmscan
    
    
    Arguments:
    ----------
    hmm_path -- (str) the path to the existing hmm file containing hidden-markov models
    domain_subset -- (list of str) a subset of "ACC" values present in the existing hmm file
    output_path -- (str) the path for the resulting file which will be created or replaced
    """

    missing_acc = set(domain_subset)

    buffer = []
    with open(hmm_path) as fin:
        with open(output_path,"w") as fout:

            input_line_index = 0
            save_section = True
            
            while True:
                line = fin.readline()

                # reached accession label
                # we may decide to give up on saving the current section
                if line.startswith("ACC"):
                    acc = line.strip().split()[1].split(".")[0]
                    if acc not in domain_subset:
                        save_section = False
                    missing_acc.discard(acc)

                # reached end of input section
                if line.startswith("//"):

                    # save the section that just ended, if necessary
                    if save_section:
                        fout.writelines(buffer)
                        fout.write("//\n")
                            
                    # clear the buffer and assume the next section will be save-able
                    buffer = []
                    save_section = True

                # reached end of input file
                if not line:
                    break

                if save_section and (not line.startswith("//")):
                    buffer.append( line )

                input_line_index += 1

    #print( f"missing accessions: {missing_acc}" )

    

def read_hmmscan_output(path):
    """
    Parse a file containing output from hmmscan
    return an instance of HmmscanResult
    """
    
    # prepare dataframe to contain data
    colnames = [
        "target name","accession","tlen","query name","_accession","qlen",
        "E-value","fs_score","fs_bias", # full sequence
        "#","of","c-Evalue","i-Evalue","td_score","td_bias", # this domain
        "hmm_from","hmm_to", # hmm coord
        "ali_from","ali_to", # ali coord
        "env_from","env_to", # env coord
        "acc"
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
                
    # remove suffixes from accession names
    df["accession"] = [s.split(".")[0] for s in df["accession"]]
                
    # convert some columns to numeric
    for col in colnames[4:]:
        df[col] = pd.to_numeric( df[col], errors="coerce" )
    
    return HmmscanResult(df, footer)