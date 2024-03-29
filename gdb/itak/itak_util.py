# this file contains utilities that are related to iTAK
#   - building input: criteria for families (TF_Rule.txt)
#   - parsing and interpreting output from iTAK

import pandas as pd


def get_gene_families( transcript_families, transcript_genes, conflict_report_path=None ):
    """
    convert from raw itak results (based on transcript IDs)
    to gene_id -> family classifications
    
    report conflicts
    
    return a dataframe with columns:
        gene_id       -- the gene_id that was categorized
        family        -- the resulting category
        transcript_id -- one arbitrary transcript that is a representative
    
    Arguments:
    ----------
    transcript_families -- (dict) transcript_id -> family
                    should be based on output from ItakRunner.run_itak()
    transcript_genes -- (dict) transcript_id -> gene_id
                    typically output from gdb.fasta.get_transcript_gene_dict()
    conflict_report_path -- (optional) (str) filepath to save a report of 
                    conflicting transcripts
    """
   
    # prepare to store conflicts and results
    conflict_count = 0
    conflict_report = ""
    df = pd.DataFrame(columns=['gene_id','family','transcript_id'])
    
    # iterate through itak results
    for tid,family in transcript_families.items():
        gid = transcript_genes[tid]
        if gid in df['gene_id']:
            ex_fam,ex_tid = df.loc[ df['gene_id'] == gid, ['family','transcript_id'] ].values[0]
            if ex_fam != family:
                conflict_count += 1
                conflict_report += "\n\t".join(["conflict:",
                        f"transcript {ex_tid} has family {ex_fam}",
                        f"transcript {tid} has family {family}\n"])
        else:
            df.loc[gid,:] = [gid,family,tid]
            
    # write report file if necessary
    if conflict_report_path is not None:
        with open(conflict_report_path, "w") as fout:
            fout.write( f"{conflict_count} pairs of transcripts had conflicting families\n\n" )
            fout.write( conflict_report )
            
    return df
    
    


def read_itak_output( output_folder ):
    """
    Get all transcript -> family matches found by iTAK
    
    This involves parsing the output file "tf_all_matches.txt"
    
    This relies on out modified fork of iTAK which should be 
    associated with input name "itak_git_repo"
    https://github.com/grotewold-lab/iTAK
    
    Traditionally iTAK outputs just one family for 
    each input sequence: "tf_classification.txt"
    
    return a dataframe with columns "transcript_id", "family", "final", where the 
    "final" column contains True for transcript-family entries that were also 
    present in the traditional iTAK output.
    """
    
    # load custom output file
    result = pd.read_table(output_folder + "/tf_all_matches.txt", 
                         header=None, names=['transcript_id','family'])
    result['final'] = False
    
    # load traditional output file
    df = pd.read_table(output_folder + "/tf_classification.txt", header = None)
    for row in df.index:
        tid,family = df.loc[row,[0,1]]
        result.loc[(result['transcript_id']==tid) & (result['family']==family),
                   "final"] = True
        
    return result
    


def build_rules_file( database_folder, rules_df ):
    """
    Set the rules that iTAK will use for clasifying transcripts 
    replace the file "TF_Rule.txt"
        
    Arguments:
    ----------
    database_folder -- (str) the path of the folder containing TF_Rule.txt
    rules_df -- (DataFrame) the criteria for classifying transcripts
                       output from gdb.hmmer.get_family_criteria()
                       must contain columns "Required", "Forbidden", and "GRASSIUS"
    """
    
    
    with open( database_folder + "/TF_Rule.txt", "w") as fout:
        for row in rules_df.index:
            req,forb,name = rules_df.loc[row,["Required","Forbidden","GRASSIUS"]]
            
            id = str(row)
            while(len(id) < 4):
                id = "0" + id
                
            if len(forb.strip()) == 0:
                forb = "NA"
                
            fout.write( "\n".join([
                f"ID:T{id}",
                f"Name:{name}",
                f"Family:{name}",
                f"Required:{req}",
                "Auxiiary:NA",
                f"Forbidden:{forb}",
                "Type:TF",
                "Desc:NA",
                "//\n\n"
            ]))