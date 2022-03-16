# this file contains functions related to the criteria for assigning families to transcripts

import pandas as pd

def read_family_criteria(path):
    """
    Load an excel sheet containing rules for assigning families to protein sequences
    
    The following columns will be considered in an automated pipeline:
        "GRASSIUS" - the family name
        "Required" - accessions that must be present (one or more accessions)
        "Forbidden" - accessions that must NOT be present (zero or more accessions)
    
    Formatting of "Required" and "Forbidden" columns is based on iTAK
    
    return a dataframe
    """
    return pd.read_excel(path).fillna('').iloc[:,:12]


    
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



def categorize_all_transcripts( acc_dict, family_criteria ):
    """
    categorize transcripts into families
    
    return a dictionary where 
        keys are family names
        values are lists of matching transcript IDs
    
    Arguments:
    ----------
    acc_dict -- (dict) summary of hmmscan results returned by get_acc_dict
                keys are transcript IDs
                values are lists of matching accession names
    family_criteria -- (Dataframe) criteria returned by read_family_criteria
    
    """
    tids_by_family = {}
    for row in crit_df.index:
        required_accs = crit_df.loc[row,"Required"].split(":")
        forbidden_accs = crit_df.loc[row,"Forbidden"].split(":")
        family_name = crit_df.loc[row,"GRASSIUS"]    
        family_tids = find_matching_transcripts( acc_dict, required_accs, forbidden_accs )
        tids_by_family[family_name] = family_tids
        
    return tids_by_family



def find_matching_transcripts( acc_dict, required_accs, forbidden_accs ):
    """
    find transcripts that fit the given criteria 
    
    the criteria must be given in iTAK format (see below)
    
    return a list of transcript IDs, a subset of the keys from acc_dict 
    
    Arguments:
    ----------
    acc_dict -- (dict) a dictionary containing a summary of hmmscan results (see util.get_acc_dict)
    required_accs -- (list of str) the list of required accessions suffixed with "#1" 
                                    or "#2" to indicate that two copies must be present
    forbidden_accs -- (list of str) the list of forbidden accessions suffixed with "#1"
    """
    result = []
    for t_id,t_accs in acc_dict.items():
        if is_match( t_accs, required_accs, forbidden_accs ):
            result.append( t_id )
    return result


def is_match( t_accs, required_accs, forbidden_accs ):
    """
    check if the given list of accessions fits the given criteria
    
    the criteria must be given in iTAK format (see below)
    
    Arguments:
    ----------
    t_accs -- (list of str) the list of accession names to check, which may contain repetition
    required_accs -- (list of str) the list of required accessions suffixed with "#1" 
                                    or "#2" to indicate that two copies must be present
    forbidden_accs -- (list of str) the list of forbidden accessions suffixed with "#1"
    """
    for acc in required_accs:
        if acc.endswith("#1"):
            if acc[:-2] not in t_accs:
                return False
        elif acc.endswith("#2"):
            if t_accs.count(acc[:-2]) < 2:
                return False
        else:
            raise Exception( "unrecognized required accession criteria \"{0}\"".format( acc ) )
    for acc in forbidden_accs:
        if acc == '':
            continue
        elif acc.endswith("#1"):
            if acc[:-2] in t_accs:
                return False
        else:
            raise Exception( "unrecognized forbidden accession criteria \"{0}\"".format( acc ) )
    return True