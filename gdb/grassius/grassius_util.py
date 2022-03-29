# this file contains utility functions that are specific to grassius

# local imports
from ..input_manager import InputManager


import pandas as pd


def assign_protein_names( gene_families, old_grassius_names, mgdb_assoc ):
    """
    Assign protein names to gene IDs in the traditional grassius form
    
    Attempt to match the names on the old grassius website
    
    e.g. "ZmbZIP1" is a protein in the bZIP family with gene ID "GRMZM2G015534"
   
    This function should be called once, so the first argument should contain
    the full set of genes to be included in the database
    
    For assigning brand new protein names, the order of the given genes
    will be used for the numerical ordering of protein names
    
    return a dictionary where keys are gene ids, values are protein names
    
    Arguments:
    ----------
    gene_families -- (dataframe) including columns "gene_id" and "family"
                        concatenated outputs from gdb.itak.get_gene_families
                        row-order may effect newly-assigned protein names
    old_grassius_names -- (dataframe) names from old grassius website
                          output from get_old_grassius_names()
    mgdb_assoc -- (dataframe) relates gene_ids between genome versions
                        output from get_maizegdb_associations()
    """
    result = {}
    
    # parse numeric suffixes from protein names (add columns)
    old_grassius_names = _parse_protein_names(old_grassius_names)
    
    # for each family, find or create a protein-name-prefix
    family_prefixes = _get_family_prefixes( gene_families, old_grassius_names )
    
    # iterate through genes to find agreement between old and new families
    n = len(gene_families.index)
    for i,row in enumerate(gene_families.index):
        
        # report progress 
        if i%100 == 0:
            print( f'checking for agreement with old grassius... ({i}/{n})' )
        
        # where there is agreement between old and new families,
        # keep protein names unchanged
        gene_id,new_family = gene_families.loc[row,["gene_id","family"]]
        if gene_id in result.keys():
            continue
        rel_ids = _get_related_gene_ids(gene_id,mgdb_assoc)
        old_match = old_grassius_names[old_grassius_names["v3_id"].isin(rel_ids)]
        if len(old_match.index)>0:    
            old_family = old_match["family"].values[0]
            if old_family == new_family:
                protein_name = old_match["name"].values[0]
                for gene_id in set(rel_ids):
                    if gene_id in gene_families["gene_id"]:
                        result[gene_id] = protein_name
          
        
    # identify protein names that were in old grassius, but are now unused
    unused_old_grassius_names = old_grassius_names[~old_grassius_names["name"].isin(result.values())]
    unused_old_grassius_names = unused_old_grassius_names.sort_values("suffix", ascending=True)
    
    # iterate through remaining genes to make new protein names
    df = gene_families[~gene_families['gene_id'].isin(result.keys())]
    n = len(df.index)
    for i,row in enumerate(df.index):
        
        # report progress 
        if i%100 == 0:
            print( f'assigning names to remaining genes... ({i}/{n})' )
        
        # check if there is an old unused protein name that can be reassigned
        gene_id,new_family = df.loc[row,["gene_id","family"]]
        if gene_id in result.keys():
            continue
        unused_old_names = unused_old_grassius_names[unused_old_grassius_names["family"] == new_family]
        if len(unused_old_names.index)>0:
            new_name = unused_old_names["name"].values[0]
            old_gene_id = unused_old_names["name"].values[0]
            all_new_gene_ids = _get_related_gene_ids(gene_id,mgdb_assoc)
            print("\n\t".join([
                f're-assigning old protein name "{new_name}"',
                f'old gene id: {old_gene_id}',
                f'new gene id(s): {all_new_gene_ids}'
            ]))
            for gid in all_new_gene_ids:
                result[gid] = new_name
            unused_old_grassius_names = unused_old_grassius_names[unused_old_grassius_names['name'] != new_name]
            continue
        
        # make up a new protein name
        
            
    
    return result


def _get_family_prefixes( gene_families, old_grassius_names ):
    """
    Find or create protein-name-prefixes for all families
    
    return a dictionary where keys are family names, 
    values are protein-name-prefixes
    
    Arguments:
    ----------
    gene_families -- (dataframe) including columns "gene_id" and "family"
                        concatenated outputs from gdb.itak.get_gene_families
    old_grassius_names -- (dataframe) parsed names from old grassius website
                          output from _parse_protein_names()
    """
    result = {}
    
    # assert that protein names have been parsed
    if 'prefix' not in old_grassius_names.columns:
        raise Exception('old_grassius_names should have been parsed using _parse_protein_names()')
              
    all_families = set(gene_families['family'])
        
    # get prefixes from old grassius
    for family_name in all_families:
        old_match = old_grassius_names[old_grassius_names['family'] == family_name]
        if len(old_match.index) > 0:
            result[family_name] = old_match['prefix'].values[0]
        
    # make new prefixes
    for family_name in set(all_families)-set(result.keys()):
        
        if family_name.startswith('AP2/'):
            new_prefix = "Zm" + family_name[4:].replace('-','')
        elif family_name == "MADS-MIKC":
            new_prefix = "ZmMIKC"
        else:
            new_prefix = "Zm" + family_name[:3].upper()
        print(f'using new prefix "{new_prefix}" for new family "{family_name}"')
        
        if new_prefix in result.values():
            raise Exception(f'new prefix "{new_prefix}" already taken')
        result[family_name] = new_prefix
            
    return result
                       
    
def _parse_protein_names( old_grassius_names ):
    """
    Extract the numerical suffixes from protein names
    
    return a new dataframe including all data from the given dataframe, 
    and two added columns: "prefix" (str), "suffix" (int)
    
    Arguments:
    ----------
    old_grassius_names -- (dataframe) names from old grassius website
                          output from get_old_grassius_names()
    """
    
    all_names = list(old_grassius_names['name'])
    all_prefixes = []
    all_suffixes = []
    
    for name in all_names:
        i = len(name)-1
        if name[i] not in '0123456789':
            raise Exception(f'protein name "{name}" does not end with a number')
        while name[i] in '0123456789':
            i -= 1
        all_prefixes.append(name[:(i+1)])
        all_suffixes.append(int(name[(i+1):]))
        
    # return as dataframe
    df = old_grassius_names.copy()
    df['prefix'] = all_prefixes
    df['suffix'] = all_suffixes
    return df
                            
                    
    
def _get_related_gene_ids( gene_id, mgdb_assoc ):
    """
    Get the set of all gene ids that are related to the
    given gene id (according to maizegdb)
    """
    
    df = mgdb_assoc
    if gene_id not in df["gene_id"].values:
        return [gene_id]
        #raise Exception(f'gene id "{gene_id}" not in maizegdb associations')
        
    name = df.loc[df["gene_id"] == gene_id,"name"].values[0]
    return df.loc[df["name"] == name, "gene_id"].values



def get_old_grassius_names():
    """
    get families, protein names, v3 gene ids
    corresponding with the old grassius website
    
    return a dataframe with columns:
    "class","family","name","synonym","v3_id"
    """
    
    return pd.read_excel( InputManager()['old_grassius_names'] )

    
def get_maizegdb_associations():
    """
    parse the uniquely-formatted public input:
    "maizegdb_gene_id_associations"
    
    return a dataframe with columns: "gene_id","name"
    where matching "name" values indicate an association
    """
    
    # read maizegdb gene_id associations    
    all_gene_ids = []
    all_names = []
    line_count = 0
    filepath = InputManager()["maizegdb_gene_id_associations"]
    with open(filepath, "r") as fin:

        # skip first line
        fin.readline()

        while True:

            line = fin.readline()
            if not line:
                break

            parts = line.split("\t")
            for gid in parts[1:]:
                gid = gid.strip()
                if len(gid) == 0:
                    continue
                elif gid.startswith("B73v1_"):
                    continue
                elif gid.startswith("B73v2_"):
                    continue
                elif gid.startswith("B73v3_"):
                    gid = gid[6:]
                elif gid.startswith("GRMZM"):
                    pass
                elif gid.startswith("Zm00001d"):
                    pass
                elif gid.startswith("Zm00001eb"):
                    pass
                else:
                    raise Exception("unrecognized prefix for gene id " + gid)

                all_gene_ids.append( gid )
                all_names.append( parts[0] )


    # return as dataframe
    return pd.DataFrame(data={
        "gene_id": all_gene_ids,
        "name": all_names
    })
    