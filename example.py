# standard imports
import pandas as pd

# local imports
import gdb
from gdb.blast import run_blast_and_annotate
from gdb.fasta import get_all_gene_ids,get_gene_id_from_record,get_records_for_gene_ids


# download and/or check integrity of all inputs
im = gdb.InputManager()
#im.prepare_all_inputs()


# get full sets of gene IDs in all maize versions
filepath = im.get_input_filepath("all_maize_gene_ids")
df = pd.read_table(filepath, header=None, names=["version","gene_id"])
id_version = {df.at[row,"gene_id"]:df.at[row,"version"] for row in df.index}
        

# find v3 gene IDs are in grassius but do not have corresponding v4 or v5 IDs
mgid_path = im.get_input_filepath("maizegdb_gene_ids")
df = pd.read_table(mgid_path, index_col=0)
missing_v4,missing_v5 = set(),set()
for row in df.index:
    alleles = df.at[row,"alleles"].split(",")
    v3_ids = [ a for a in alleles if id_version[a] == "v3" ]
    versions_present = set( id_version[a] for a in alleles )
    if "v3" not in versions_present:
        raise Exception( "no v3 id for protein " + df.at[row,"grassius_name"] )
    if "v4" not in versions_present:
        missing_v4.update(v3_ids)
    if "v5" not in versions_present:
        missing_v5.update(v3_ids)
        
missing_v4_and_v5 = missing_v4.intersection(missing_v5)


# blast v3 proteins that are missing v5 IDs
result_df = run_blast_and_annotate( 
    needle_dna_fasta   = im.get_input_filepath("maize_v3_cdna"), 
    haystack_dna_fasta = im.get_input_filepath("maize_v4_full_fasta"), 
    haystack_gff       = im.get_input_filepath("maize_v4_gff3"), 
    gene_id_subset     = missing_v4 )
result_df.to_csv("v3_missing_v4.csv", index=False)


# s
