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
from gdb.patches import *


im = gdb.InputManager()
old_grassius_names = get_old_maize_grassius_names()
old_grassius_tfomes = get_old_grassius_tfomes()
mgdb_assoc = get_maizegdb_associations()
transcript_genes = get_transcript_gene_dict( im['maize_v5_proteins'] )
gene_interactions = pd.read_excel(im['gene_interactions'])
domain_descriptions = get_domain_descriptions()
domain_annotations = get_domain_annotations()
    

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
    
    # run itak with all rules
    desired_accessions = get_relevant_accessions(family_criteria_df)
    build_minified_hmm( im["pfam_hmm"], desired_accessions, "pfam_min.hmm" )
    concatenate_hmms( ["pfam_min.hmm",im["selfbuild_hmm"]], "combined.hmm" )
    ir = ItakRunner(reset=True)
    ir.set_database( "combined.hmm", family_criteria_df )
    itak_results = ir.run_itak( 'test.fa' )
    #itak_results = ir.run_itak( im["maize_v3_proteins"] )

    # save results
    itak_results.to_csv("itak_results.txt", sep="\t", index=False)

    raise Exception('test')

# load MAIZE results as if we had just run iTAK (above)
itak_results = pd.read_csv('applied_rules.csv')

# based on iTAK results, pick one family for each transcript
# give certain families priority
# give MYB-Related priority over ARR-B and MYB
# except for cases where there is agreement with old grassius
all_transcripts = set(itak_results['transcript_id'])
transcript_families = {}
priority_families = ['AP2/ERF-AP2','MYB']
df = old_grassius_names
old_myb_gids = df.loc[ df['family']=='MYB', 'v3_id' ].values
old_arrb_gids = df.loc[ df['family']=='ARR-B', 'v3_id' ].values
for tid in all_transcripts:
    matched_families = itak_results.loc[
            itak_results['transcript_id'] == tid,
            'family'].values
    gid = transcript_genes[tid]
    matching_priority_families = [fam for fam in priority_families if fam in matched_families]
    
    # check for special cases
    if ('MYB' in matched_families) and (gid in old_myb_gids):
        transcript_families[tid] = 'MYB'
    elif ('ARR-B' in matched_families) and (gid in old_arrb_gids):
        transcript_families[tid] = 'ARR-B'
        
    # check for prioritized families
    elif len(matching_priority_families) > 0:
        transcript_families[tid] = matching_priority_families[0]
        
    # default to an arbitrary matching family
    else:
        transcript_families[tid] = list(matched_families)[0]
        
# convert itak results (based on transcript IDs)
# to gene -> family classifications
gene_families = get_gene_families( transcript_families, transcript_genes, "conflicts.txt" )


# set order of RAV family (previously called 'ABI3-VP1')
# based on old-grassius 'ABI3-VP1' protein names
# this will effect the suffixes of the new protein names
df = gene_families
df1 = df[df['family'] == 'RAV'].copy()
df2 = df[df['family'] != 'RAV'].copy()
for row in df1.index:
    gid = df1.loc[row,"gene_id"]
    old_match = old_grassius_names[old_grassius_names["v3_id"] == gid]
    
    if len(old_match.index) == 0:
        order = 9999
    else:
        old_family = old_match['family'].values[0]
        if old_family != 'ABI3-VP1':
            order = 9999
        else:
            order = old_match['suffix'].values[0]
    df1.loc[row,'order'] = order
df1 = df1.sort_values('order').drop(columns=['order'])
gene_families = pd.concat([df1,df2])
gene_families.sort_values("family").to_csv("gene_families.csv", index=False)


# assign maize protein names
protein_names = assign_protein_names( gene_families, old_grassius_names, mgdb_assoc, 
                            report_folder = "." )


# build maize metadata 
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
    
    
# assign non-maize protein names
old_nonmaize_names = get_old_maize_grassius_names()
for input_name in ['2023_brachy_families']:
    
    # load results from Dec2023 grass iTAK pipeline
    # for one species
    new_families = pd.read_table( im[input_name] )
    
    # build metadata for one non-maize species
    brachy_df = assign_protein_names( 
        new_families, old_nonmaize_names, 
        mgdb_assoc=None, report_folder="." )

    # apend metadata for one non-maize species
    df = pd.concat([
        df,
        brachy_df
    ], ignore_index=True).fillna('')
    
    
    

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
                         gene_interactions, domain_descriptions )


    

# insert maize sequences from fasta files
for suffix in ["cdna","proteins"]:
    for version in ["v3","v4","v5"]:
        organism = f"Maize_{version}"
        fasta_filepath = im[f"maize_{version}_{suffix}"]
        cb.insert_sequences( organism, df, fasta_filepath, 
                            is_protein=(suffix=='proteins') )

        
# insert non-maize sequences from fasta files
for organism in ['Rice','Sorghum','Sugarcane','Brachypodium']:
    for suffix in ["cdna","proteins"]:
        fasta_filepath = im[f"{organism}_{suffix}".lower()]
        cb.insert_sequences( organism + "_", df, fasta_filepath, 
                            is_protein=(suffix=='proteins') )
            
# insert Jan2022 secondary structure
cb.insert_secondary_structures()

cb.insert_domain_annotations( domain_annotations )
        

# add tfome sequences
cb.insert_tfomes( old_grassius_tfomes )

# apply patches
apply_all_patches(cb)

# save a snapshot of the database that was built
cb.write_snapshot( "build_db.sql.tar.gz" )