
import pandas as pd

def read_itak_output( output_folder ):
    """
    read itak classification results
    
    return a dictionary where keys are transcript IDs,
    values are classifications
    """
    df = pd.read_table(output_folder + "/tf_classification.txt", header = None)
    return { df.loc[row,0]:df.loc[row,1] for row in df.index }


def build_rules_file( database_folder, rules_df ):
    """
    Set the rules that iTAK will use for clasifying transcripts 
    replace the file "TF_Rule.txt"
        
    Arguments:
    ----------
    database_folder -- (str) the path of the folder containing TF_Rule.txt
    rules_df -- (DataFrame) the criteria for classifying transcripts
                       output from gdb.hmmer.read_family_criteria()
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