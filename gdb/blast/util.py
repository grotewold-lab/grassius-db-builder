from .blast_result import BlastResult
import pandas as pd
import tempfile
import os

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
        
        
        

def read_blast_output(path):
    """
    Parse a file containing output from TBLASTN
    return an instance of BlastResult
    """
    
    # load the whole file
    with open(path,'r') as fin:
        lines = fin.readlines()
        
    # identify section breaks
    data_start_index, data_end_index = _find_section_breaks(lines)
    
    # return empty result if there is no data
    if data_start_index is None:
        return BlastResult("\n".join(lines),None,None)
    
    # parse data
    data_lines = lines[data_start_index:data_end_index]
    data = _parse_data(data_lines)
    
    # prepare text
    header = "\n".join(lines[:data_start_index])
    footer = "\n".join(lines[data_end_index:])
    
    return BlastResult(header,data,footer)



def _find_section_breaks(lines):
    data_start_index = None
    data_end_index = None
    for i,line in enumerate(lines):
        if (data_start_index is None) and line.startswith(">"):
            data_start_index = i
        if (data_start_index is not None) and line.startswith("lambda"):
            data_end_index = i
            
    return data_start_index,data_end_index
    
    
def _parse_data(data_lines):
    result = pd.DataFrame(columns=["Chrom","Positions","Identities%"])
    entry_index = -1
    chrom = ""
    
    for line in data_lines:
    
        # start of a new chromosome
        if line.startswith(">"):
            chrom = line.strip()[1:]
    
        # start of a new entry
        if line.startswith(" Score"):
            entry_index += 1
            result.at[entry_index,"Chrom"] = chrom
            result.at[entry_index,"Positions"] = []
            
        # second line of new entry
        if line.startswith(" Identities"):
            result.at[entry_index,"Identities%"] = int(line.split()[3][1:-3])
            
        # start of genomic coordinates 
        if line.startswith("Sbjct"):
            parts = line.split()
            for pi in 1,3:
                pos = int(parts[pi])
                result.at[entry_index,"Positions"].append(pos)
                
    return result
                
            
        
    