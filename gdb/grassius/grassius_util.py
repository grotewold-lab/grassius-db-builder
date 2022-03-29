# this file contains utility functions that are specific to grassius

# local imports
from ..input_manager import InputManager


import pandas as pd


def assign_protein_names( gene_families, old_grassius_names, mgdb_assoc ):
    """
    Assign protein names in the traditional grassius form
    
    e.g. "ZmbZIP1" is the name of a protein in the bZIP family
    
    Attempt to match the names on the old grassius website
    
    return a dictionary where keys are gene ids, values are protein names
    
    Arguments:
    ----------
    gene_families -- (dataframe) including columns "gene_id" and "family"
                        output from gdb.itak.get_gene_families
    old_grassius_names -- (dataframe) names from old grassius website
                          output from get_old_grassius_names()
    mgdb_assoc -- (dataframe) relates gene_ids between genome versions
                        output from get_maizegdb_associations()
    """
    
    
    result = {}
    n = len(gene_families.index)
    
    # iterate through genes
    for i,row in enumerate(gene_families.index):
        
        # report progress 
        if i%100 == 0:
            print( f'assigning protein names to genes... ({i}/{n})' )
        
        # where there is agreement between old and new families,
        # keep protein names unchanged
        gene_id,new_family = gene_families.loc[row,["gene_id","family"]]
        rel_ids = _get_related_gene_ids(gene_id,mgdb_assoc)
        old_match = old_grassius_names[old_grassius_names["v3_id"].isin(rel_ids)]
        if len(old_match.index)>0:    
            old_family = old_match["family"].values[0]
            if old_family == new_family:
                result[gene_id] = old_match["name"].values[0]
                continue
                
        # TODO make up new protein names
    
    return result
                     
    
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
    