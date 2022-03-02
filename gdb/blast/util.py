from .blast_result import BlastResult
import pandas as pd
import os


def annotate_blast_result( blast_result, gff3_data ):
    """
    Add annotations to a blast result 
    
    This will add a new column "Gene_ID" to blast_result.data 
    """

    df = blast_result.data
    
    for row in df.index:
        chrom,start_pos,stop_pos = df.loc[row,["Chrom","Start_Pos","Stop_Pos"]]
        df.loc[row,"Gene ID"] = _get_matching_gene_id( chrom,start_pos,stop_pos, gff3_data )
        
    
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
                
    # now the result may have >2 genomic positions for each row
    # summarize all genomic coordinates into two values: start and stop
    result['Start_Pos'] = [min(result.at[row,"Positions"]) for row in result.index]
    result['Stop_Pos'] = [max(result.at[row,"Positions"]) for row in result.index]
    result = result[["Chrom","Start_Pos","Stop_Pos","Identities%"]]
                
    return result
                
    
        
def _get_matching_gene_id( chrom, start_pos, stop_pos, df ):
    df = df[df['chrom'] == chrom]
    df_sub = df[
        ((start_pos >= df['start']) & (stop_pos <= df['end'])) |  # blast region inside gff region
        ((df['start'] >= start_pos) & (df['end'] <= stop_pos)) |  # gff region inside blast region
        ((df['start'] >= start_pos) & (df['end'] <= start_pos))| # regions partially overlap
        ((df['start'] >= stop_pos) & (df['end'] <= stop_pos))     # regions partially overlap
    ]
    if len(df_sub) > 0:
        identifiers = df_sub["identifiers"].values[0]
        assert identifiers.startswith("ID=")
        return identifiers.split(";")[0][3:]
        
    return ""
            
        
    