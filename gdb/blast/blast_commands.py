# this file contains wrappers for the blast commands "makeblastdb" and "tblastn"

from .util import read_blast_output
import subprocess
from os.path import dirname,basename
import tempfile
import shutil


def prepare_blast_db( fasta_path ):
    """
    run makeblastdb so that the given fasta file can be used
    with the tblastn argument "-db"
    """
    
    fasta_dir = dirname(fasta_path)
    fasta_name = basename(fasta_path)
    
    command = [
        "makeblastdb",
        "-in", fasta_name,
        "-parse_seqids",
        "-dbtype", "nucl",
        #"-dbtype", "prot",
    ]
    
    p = subprocess.Popen(command, cwd=fasta_dir )
    p.wait()
    
    
    
def run_blastn( dna_sequence, fasta_path ):
    """
    Search for fasta DNA entries matching the given DNA sequence
    
    
    Arguments:
    ----------
    dna_sequence -- (str) one protein sequence to search for
    fasta_filepath -- (str) the path to a dna fasta file to search in
    
    
    return an instance of BlastResult
    """
    
    # make a temprorary folder
    folder = tempfile.mkdtemp()
    seq_path = folder+"/seq.fa"
    out_path = folder+"/blast_output.txt"
    
    # build fasta file
    with open(seq_path, "w") as f:
        f.writelines([">\n", dna_sequence])

    # run blast
    command = [
        "blastn",
        "-query", seq_path,
        "-out", out_path,
        "-db", fasta_path
    ]
    p = subprocess.Popen(command)
    p.wait()
    
    result = read_blast_output(out_path)

    shutil.rmtree(folder)
    
    return result