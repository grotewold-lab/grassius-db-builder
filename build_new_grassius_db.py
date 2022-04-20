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
old_grassius_names = get_old_grassius_names()
old_grassius_tfomes = get_old_grassius_tfomes()
mgdb_assoc = get_maizegdb_associations()
transcript_genes = get_transcript_gene_dict( im['maize_v3_proteins'] )
gene_interactions = pd.read_excel(im['gene_interactions'])
        
    

# load family criteria and descriptions
family_criteria_df = get_family_criteria()
family_desc_df = get_family_descriptions()


# compare new family names with old family names
#old_df = pd.read_excel( im['old_grassius_names'] )
#old_df = old_df[old_df['family'] != 'Orphans']
#old_families = set(old_df['family'])
#new_families = set(family_criteria_df["GRASSIUS"])

if False:
    # build subset of maize v3 protein fasta
    # containing only genes that where myb-related in old grassius
    df = old_grassius_names
    mybr_gene_ids = df.loc[ df['family']=='MYB-related', 'v3_id' ].values
    with open('test.fa','w') as fout:
        for rec in read_records_for_gene_ids( im['maize_v3_proteins'], mybr_gene_ids ):
            fout.write( ">{0}\n{1}\n".format( rec.id, str(rec.seq) ) )

    raise Exception('test')

if False:
            
    # investigate issue with myb-related family
    # run itak once and exclude some families that seem to
    # take prioirity over MYB-related
    df = family_criteria_df
    mod_family_criteria_df = df[~df['GRASSIUS'].isin(['MYB','ARR-B'])]
    desired_accessions = get_relevant_accessions(mod_family_criteria_df)
    build_minified_hmm( im["pfam_hmm"], desired_accessions, "pfam_min.hmm" )
    ir = ItakRunner(reset=True)
    ir.set_database( "pfam_min.hmm", mod_family_criteria_df )
    mod_itak_results = ir.run_itak( im['maize_v3_proteins'] )

    # save results
    with open( "mod_itak_results.txt", "w" ) as fout:
        for key,value in mod_itak_results.items():
            fout.write( f"{key}\t{value}\n" )

    raise Exception('test')

if False:
    
    # run itak with all rules
    desired_accessions = get_relevant_accessions(family_criteria_df)
    build_minified_hmm( im["pfam_hmm"], desired_accessions, "pfam_min.hmm" )
    concatenate_hmms( ["pfam_min.hmm",im["selfbuild_hmm"]], "combined.hmm" )
    ir = ItakRunner(reset=True)
    ir.set_database( "combined.hmm", family_criteria_df )
    itak_results = ir.run_itak( im["maize_v3_proteins"] )

    # save results
    with open( "itak_results.txt", "w" ) as fout:
        for key,value in itak_results.items():
            fout.write( f"{key}\t{value}\n" )

    raise Exception('test')


# load premade results as if we had used ItakRunner.run_itak
def load_itak_results( filename ):
    itak_results = {}
    with open( filename ) as fin:
        while True:
            line = fin.readline()
            if not line:
                break
            parts = line.split("\t")
            itak_results[parts[0]] = parts[1].strip()
    return itak_results

mod_itak_results = load_itak_results('mod_itak_results.txt')
itak_results = load_itak_results('itak_results.txt')


# modify itak results
# give myb-related priority over ARR-B and MYB
# except for cases where there is agreement with old grassius

df = old_grassius_names
old_myb_gids = df.loc[ df['family']=='MYB', 'v3_id' ].values
old_arrb_gids = df.loc[ df['family']=='ARR-B', 'v3_id' ].values
new_mybr_tids = [k for k,v in mod_itak_results.items() if v=='MYB-related']
for tid in new_mybr_tids:
    gid = transcript_genes[tid]
    family = itak_results[tid]
    if family not in ['MYB','ARR-B']:
        raise Exception(f'unexpected family "{family}"')
    if (family=='MYB') and (gid in old_myb_gids):
        continue
    if (family=='ARR-B') and (gid in old_arrb_gids):
        continue
    itak_results[tid] = 'MYB-related'
    

        
# convert itak results (based on transcript IDs)
# to gene -> family classifications
gene_families = get_gene_families( itak_results, transcript_genes, "conflicts.txt" )


# assign protein names
protein_names = assign_protein_names( gene_families, old_grassius_names, mgdb_assoc, 
                            report_folder = "." )


# build metadata 
df = protein_names
for row in df.index:
    gid,name,family = df.loc[row,["gene_id","name","family"]]
    if family == 'Orphans':
        clazz = 'Orphans'
    else:
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
    
    
# start a new database or connect to an existing database
cb = ChadoBuilder()

# create non-chado tables
cb.build_grassius_tables( df, gene_versions, family_desc_df, 
                         old_grassius_names, old_grassius_tfomes, 
                         gene_interactions )


# insert sequences from fasta files
for suffix in ["cdna","proteins"]:
    for version in ["v3","v4","v5"]:
        organism = f"Maize_{version}"
        fasta_filepath = im[f"maize_{version}_{suffix}"]
        cb.insert_sequences( organism, df, fasta_filepath, 
                            is_protein=(suffix=='proteins') )

# insert Jan2022 secondary structure
cb.insert_secondary_structures()
        

# add tfome sequences
cb.insert_tfomes( old_grassius_tfomes )


# save a snapshot of the database that was built
cb.write_snapshot( "build_db.sql.tar.gz" )