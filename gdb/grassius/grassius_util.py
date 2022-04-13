# this file contains utility functions that are specific to grassius

# local imports
from ..input_manager import InputManager


import pandas as pd


def assign_protein_names( gene_families, old_grassius_names, mgdb_assoc, report_folder=None ):
    """
    Assign protein names to gene IDs in the traditional grassius form
    
    Attempt to match the names on the old grassius website
    
    e.g. "ZmbZIP1" is a protein in the bZIP family with gene ID "GRMZM2G015534"
   
    This function should be called once, so the first argument should contain
    the full set of genes to be included in the database
    
    The order of the given genes (order of rows in gene_families) may 
    affect the results:
        In cases where gene-family categories conflict with 
        maizegdb gene-gene associations, the gene-family that
        was encountered first takes priority 
        
        numeric suffixes for protein names (only when a brand 
        new protein name is created)
    
    return a dataframe with columns "gene_id", "name", "family"
    
    Arguments:
    ----------
    gene_families -- (dataframe) including columns "gene_id" and "family"
                        concatenated outputs from gdb.itak.get_gene_families
                        row-order may effect newly-assigned protein names
    old_grassius_names -- (dataframe) names from old grassius website
                          output from get_old_grassius_names()
    mgdb_assoc -- (dataframe) relates gene_ids between genome versions
                        output from get_maizegdb_associations()
    report_folder -- (optional) (str) path to a folder to save text reports
    """
    result = pd.DataFrame(columns=["gene_id","name","family"])
    
    # prepare to write report if necessary
    if report_folder is None:
        fout_reas = None
        fout_conf = None
    else:
        fout_reas = open(report_folder + "/reassigned_names.txt", "w")
        fout_conf = open(report_folder + "/assoc_conflicts.txt", "w")
    
    # for each family, find or create a protein-name-prefix
    family_prefixes = _get_family_prefixes( gene_families, old_grassius_names, report_folder )
    
    # for each family, pick the next numerical value for assigning new protein names
    family_next_suffixes = {}
    for family in set(gene_families['family']):
        old_match = old_grassius_names[old_grassius_names["family"] == family]
        if len(old_match.index) > 0:
            next_suffix = old_match["suffix"].max()+1
        else:
            next_suffix = 1
        family_next_suffixes[family] = next_suffix
    
    # iterate through genes to find agreement between old and new families
    n = len(gene_families.index)
    for i,row in enumerate(gene_families.index):
        
        # report progress 
        if i%100 == 0:
            print( f'checking for agreement with old grassius... ({i}/{n})' )
        
        # where there is agreement between old and new families,
        # keep protein names unchanged
        gene_id,new_family = gene_families.loc[row,["gene_id","family"]]
        if gene_id in result["gene_id"]:
            continue
        rel_ids = _get_related_gene_ids(gene_id,mgdb_assoc)
        old_match = old_grassius_names[old_grassius_names["v3_id"].isin(rel_ids)]
        if len(old_match.index)>0:    
            old_family = old_match["family"].values[0]
            if old_family == new_family:
                protein_name = old_match["name"].values[0]
                for gid in set(rel_ids):
                    if gid in result["gene_id"]:
                        conf_name = result.loc[result["gene_id"]==gid,"name"].values[0]
                        _report(fout_conf,"\n\t".join([
                            f'\ngene id "{gene_id}" has traditional name "{protein_name}"',
                            f'but assocated gene_id "{gid}" aleady has name "{conf_name}"'
                        ]))
                    else:
                        result.loc[gid,:] = [gid,protein_name,new_family]
          
    
    # iterate through remaining genes to make new protein names
    df = gene_families[~gene_families['gene_id'].isin(result.keys())]
    n = len(df.index)
    for i,row in enumerate(df.index):
        
        # report progress 
        if i%100 == 0:
            print( f'assigning names to remaining genes... ({i}/{n})' )
        
        # get related gene IDs
        gene_id,new_family = df.loc[row,["gene_id","family"]]
        if gene_id in result["gene_id"]:
            continue
        rel_ids = _get_related_gene_ids(gene_id,mgdb_assoc)
        
        # make up a new protein name
        prefix = family_prefixes[new_family]
        suffix = family_next_suffixes[new_family]
        family_next_suffixes[new_family] = suffix + 1
        new_name = prefix + str(suffix)
        for gid in set(rel_ids):
            if gid in result["gene_id"]:
                conf_name = result.loc[result["gene_id"]==gid,"name"].values[0]
                _report(fout_conf,"\n\t".join([
                    f'\ngene id "{gene_id}" has new name "{new_name}"',
                    f'but assocated gene_id "{gid}" aleady has name "{conf_name}"'
                ]))
            else:
                result.loc[gid,:] = [gid,new_name,new_family]
    
    return result


def _get_family_prefixes( gene_families, old_grassius_names, report_folder=None ):
    """
    Find or create protein-name-prefixes for all families
    
    return a dictionary where keys are family names, 
    values are protein-name-prefixes
    
    Arguments:
    ----------
    gene_families -- (dataframe) including columns "gene_id" and "family"
                        concatenated outputs from gdb.itak.get_gene_families
    old_grassius_names -- (dataframe) parsed names from old grassius website
                          output from parse_protein_names()
    report_folder -- (optional) (str) path to a folder to save text reports
    """
    result = {}
    
    # prepare to write report if necessary
    if report_folder is None:
        fout = None
    else:
        fout = open(report_folder + "/new_prefixes.txt", "w")
    
    # assert that protein names have been parsed
    if 'prefix' not in old_grassius_names.columns:
        raise Exception('old_grassius_names should have been parsed using parse_protein_names()')
              
    all_families = set(gene_families['family'])
        
    # get prefixes from old grassius
    for family_name in all_families:
        old_match = old_grassius_names[old_grassius_names['family'] == family_name]
        if len(old_match.index) > 0:
            result[family_name] = old_match['prefix'].values[0]
        
    # make new prefixes
    for family_name in set(all_families)-set(result.keys()):
        
        if family_name.startswith('AP2/'):
            new_prefix = "Zm" + family_name[4:].replace('-','') + "_"
        elif family_name == "MADS-MIKC":
            new_prefix = "ZmMIKC"
        else:
            new_prefix = "Zm" + family_name[:3].upper()
            
        _report(fout, f'using new prefix "{new_prefix}" for new family "{family_name}"')
        
        if new_prefix in result.values():
            raise Exception(f'new prefix "{new_prefix}" already taken')
        result[family_name] = new_prefix
            
    return result
    
def parse_protein_names( df ):
    """
    Extract the numerical suffixes from protein names
    
    return a new dataframe including all data from the given dataframe, 
    and two added columns: "prefix" (str), "suffix" (int)
    
    Arguments:
    ----------
    df -- (dataframe) with column 'name', containing
                            protein names
    """
    
    all_names = list(df['name'])
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
    df = df.copy()
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
    "class","family","name","accepted","prefix","suffix","synonym","v3_id"
    """
    
    old_grassius_names = pd.read_excel( InputManager()['old_grassius_names'] )
    old_grassius_names = parse_protein_names(old_grassius_names)
    return old_grassius_names


def get_old_grassius_tfomes():
    """
    get tfome sequences and metadata from the old grassius website
    
    return a dataframe with columns including:
    "utname","gene_id","sequence","translation"
    """
    
    return pd.read_table( InputManager()['old_grassius_tfomes'] ).fillna('')

        
    
def get_protein_name_dict( metadata_df ):
    """
    get a dictionary where keys are gene_ids, and
    values are protein names
    
    Arguments:
    ----------
    metadata_df -- (Dataframe) with columns 'gene_id' 
                                and 'name'
    """
    
    df = metadata_df
    return {
        df.loc[row,"gene_id"] : df.loc[row,"name"]
        for row in df.index
    }


def get_maize_v3_uniprot_ids():
    """
    Get uniprot IDs corresponding with maize v3 gene IDs
    
    This involves parsing some parts of a gff3 file,
    which can be done using BCBio but this is faster
    
    return a dictionary where keys are v3 gene IDs,
    and values are uniprot IDs
    """
    
    path = InputManager()['maize_v3_gff3']
    
    result = {}
    gid_search_str = 'gene:'
    uid_search_str = 'UniProtKB/TrEMBL%3BAcc:'
    
    with open(path) as f:
        while True:
            line = f.readline()
            if not line:
                break
                
            if (gid_search_str in line) and (uid_search_str in line):
                i = line.index(gid_search_str) + len(gid_search_str)
                j = line.index(';',i)
                gene_id = line[i:j]
                
                i = line.index(uid_search_str) + len(uid_search_str)
                j = line.index(']', i)
                uniprot_id = line[i:j]
                
                result[gene_id] = uniprot_id
                
    return result
        

    
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
    
                       
    
def _report( fout, text ):
    """
    print some text to the console,
    and also to the given output stream if it is not None
    """
    
    if fout is not None:
        fout.write(text + "\n")
    print(text)