# this file contains utility functions to support hmmer

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
    
    result = {}
    for row in hmm_df.index:        
        tid,acc = hmm_df.loc[row,"query name"]
        acc = hmm_df.loc[row,"accession"].split(".")[0]
        if tid not in result.keys():
            result[tid] = []
        result[tid].append( acc )
    return result

    
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
                if line.startswith("ACC"):
                    acc = line.strip().split()[1].split(".")[0]
                    if acc not in domain_subset:
                        save_section = False
                    missing_acc.discard(acc)

                # reached end of input section
                if line.startswith("//"):

                    # debug
                    #print( f"section break at line {input_line_index}" )
                    #if save_section:
                    #    print( f"saving last section (accession {acc}) to minified output" )

                    if save_section:
                        fout.writelines(buffer)
                    buffer = []
                    save_section = True

                # reached end of input file
                if not line:
                    break

                if save_section:
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
        "target name","accession","query name","_accession",
        "fs_E-value","fs_score","fs_bias",
        "best_E-value","best_score","best_bias",
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
                
    # remove suffixes from accession names
    df["accession"] = [s.split(".")[0] for s in df["accession"]]
                
    # convert some columns to numeric
    for col in colnames[4:]:
        df[col] = pd.to_numeric( df[col], errors="coerce" )
    
    return HmmscanResult(df, footer)