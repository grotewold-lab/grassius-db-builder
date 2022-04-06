"""
build a database compatible with the grassius website
using new protein naming pipeline discussed Mar31-2022
"""

import tempfile
import shutil
import os
import pickle
import pandas as pd

# local imports
import gdb
from gdb.fasta import *
from gdb.hmmer import *
from gdb.itak import *
from gdb.grassius import *
from gdb.chado import *


im = gdb.InputManager()



        

# load family criteria and descriptions
family_criteria_df = read_family_criteria(im["family_rules"])
family_desc_df = pd.read_csv(im['family_descriptions'])

# compare new family names with old family names
#old_df = pd.read_excel( im['old_grassius_names'] )
#old_df = old_df[old_df['family'] != 'Orphans']
#old_families = set(old_df['family'])
#new_families = set(family_criteria_df["GRASSIUS"])

# run itak against maize proteins
#desired_accessions = get_relevant_accessions(family_criteria_df)
#build_minified_hmm( im["pfam_hmm"], desired_accessions, "pfam_min.hmm" )
#concatenate_hmms( ["pfam_min.hmm",im["selfbuild_hmm"]], "combined.hmm" )
#ir = ItakRunner(reset=True)
#ir.set_database( "combined.hmm", family_criteria_df )
#itak_results = ir.run_itak( im["maize_v3_proteins"] )
#with open( "itak_results.txt", "w" ) as fout:
#    for key,value in itak_results.items():
#        fout.write( f"{key}\t{value}\n" )


# load premade results as if we had just run itak (above)
itak_results = {}
with open( "itak_results.txt" ) as fin:
    while True:
        line = fin.readline()
        if not line:
            break
        parts = line.split("\t")
        itak_results[parts[0]] = parts[1].strip()

        
# convert itak results (based on transcript IDs)
# to gene -> family classifications
transcript_genes = get_transcript_gene_dict( im['maize_v3_proteins'] )
gene_families = get_gene_families( itak_results, transcript_genes, "conflicts.txt" )


# assign protein names
old_grassius_names = get_old_grassius_names()
mgdb_assoc = get_maizegdb_associations()
protein_names = assign_protein_names( gene_families, old_grassius_names, mgdb_assoc, 
                            report_folder = "." )


# build metadata 
df = protein_names
for row in df.index:
    gid,name,family = df.loc[row,["gene_id","name","family"]]
    raw_clazz = family_criteria_df.loc[family_criteria_df["GRASSIUS"]==family,"category"].values[0]
    clazz = "Coreg" if (raw_clazz == "coregulators") else "TF"
    df.loc[row,"class"] = clazz
df.sort_values("name").to_csv("metadata.csv", index=False)
    
    
# build gene_id -> genome_version dictionary
gene_versions = {}
for version in ['v3','v4','v5']:
    all_gene_ids = get_all_gene_ids( im[f'maize_{version}_proteins'] )
    for gid in all_gene_ids:
        gene_versions[gid] = version
    
    
# start building database
cb = ChadoBuilder()
cb.build_grassius_tables( df, gene_versions, family_desc_df )

    
# insert sequences from fasta files
for suffix in ["cdna","proteins"]:
    for version in ["v3","v4","v5"]:
        organism = f"Maize_{version}"
        fasta_filepath = im[f"maize_{version}_{suffix}"]
        cb.insert_sequences( organism, df, fasta_filepath, is_protein=False )

        
# save a snapshot of the database that was built
cb.write_snapshot( "build_db.sql.tar.gz" )