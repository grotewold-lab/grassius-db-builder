# this file contains wrappers for hmmer commands like "hmmscan"

from .util import read_hmmscan_output
import subprocess
import tempfile

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
        "--tblout", out_path,
        hmm_path,
        fasta_path
    ]
    
    p = subprocess.Popen(command)
    p.wait()
    
    return read_hmmscan_output(out_path)