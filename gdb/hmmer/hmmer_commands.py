# this file contains wrappers for hmmer commands like "hmmscan"

# local imports
from .hmmer_util import read_hmmscan_output


import subprocess
import tempfile
import shutil

def run_hmmpress( hmm_path ):
    """
    prepare the given hmm file for use with hmmscan
    """
    
    command = [ "hmmpress", hmm_path ]
    p = subprocess.Popen(command)
    p.wait()
    

def run_hmmscan( hmm_path, fasta_path ):
    """
    identify domains in protein sequences
    by executing the command "hmmscan" and parsing the results
    
    Arguments:
    ----------
    hmm_path -- (str) the path to an hmm file containing hidden-markov models
    fasta_path -- (str) the path to a fasta file containing protein sequences
    
    return an instance of HmmscanResult
    """
    
    # make a temporary folder
    folder = tempfile.mkdtemp()
    out_path = folder+"/hmmscan_output.txt"
    
     
    command = [
        "hmmscan",
        "--domtblout", out_path,
        hmm_path,
        fasta_path
    ]
    
    p = subprocess.Popen(command)
    p.wait()
    
    result = read_hmmscan_output(out_path)
    shutil.rmtree(folder)
    
    return result