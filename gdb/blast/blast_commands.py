# this file contains wrappers for the blast commands "makeblastdb" and "tblastn"

from .util import read_blast_output
import subprocess
from os.path import dirname,basename
import os
import tempfile


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
    
    
    
def run_tblastn( protein_sequence, fasta_path ):
    """
    Search for fasta DNA entries matching the given protein sequence
    
    
    Arguments:
    ----------
    protein_sequence -- (str) one protein sequence to search for
    fasta_filepath -- (str) the path to a dna fasta file to search in
    
    
    return an instance of BlastResult
    """
    
    # make a temprorary folder
    folder = tempfile.mkdtemp()
    seq_path = folder+"/seq.fa"
    out_path = folder+"/blast_output.txt"
    
    # build fasta file
    with open(seq_path, "w") as f:
        f.writelines([">\n", protein_sequence])

    # run blast
    os.system(f"tblastn -query {seq_path} -out {out_path} -db {fasta_path}")

    return read_blast_output(out_path)