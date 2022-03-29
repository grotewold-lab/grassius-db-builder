
"""
judge the quality of the HMMs and iTAK pipeline

attempt to re-create old grassius family categories
"""

# local imports
import gdb
from gdb.fasta import *
from gdb.hmmer import *
from gdb.itak import *
from gdb.grassius import *


import tempfile
import shutil
import os
import pickle
import pandas as pd



im = gdb.InputManager()



        

# load family criteria
family_criteria_df = read_family_criteria(im["family_rules"])

        
# compare new family names with old family names
old_df = pd.read_excel( im['old_grassius_names'] )
old_df = old_df[old_df['family'] != 'Orphans']
old_families = set(old_df['family'])
new_families = set(family_criteria_df["GRASSIUS"])

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
df = get_gene_families( itak_results, transcript_genes, "conflicts.txt" )

        

# DEBUG test assigning protein names
gene_families = df
old_grassius_names = get_old_grassius_names()
mgdb_assoc = get_maizegdb_associations()
test = assign_protein_names( gene_families, old_grassius_names, mgdb_assoc )
raise Exception('test')

    
    
# compare new results with old grassius families
old_df = pd.read_excel( im['old_grassius_names'] )
old_df = old_df[old_df['family'] != 'Orphans']
correct_count = 0
changed_family_count = 0
changed_family_report = ""
all_old_ids = set(old_df['v3_id'])
all_new_ids = set(df['gene_id'])
missing_ids = (all_old_ids - all_new_ids)
new_ids = (all_new_ids - all_old_ids)
common_ids = all_old_ids.intersection(all_new_ids)
for gid in common_ids:
    old_fam = old_df.loc[ old_df['v3_id'] == gid, 'family'].values[0]
    new_fam = df.loc[ df['gene_id'] == gid, 'family'].values[0]
    if old_fam == new_fam:
        correct_count += 1
    else:
        changed_family_count += 1
        changed_family_report += "\n\t".join([f"changed family for {gid}:",
                f"old family was {old_fam}",
                f"new family is {new_fam}\n\n"])
        
        
# save reports
df.to_csv('all_new_families.csv',index=False)
with open("changed_families.txt", "w") as fout:
    fout.write( f"{changed_family_count} genes changed families\n\n" )
    fout.write( changed_family_report )
with open("missing_genes.txt", "w") as fout:
    fout.write( "{0} genes were not categorized into families\n".format(len(missing_ids)) )
    fout.write( "these are genes that were categorized in old grassius\n\n" )
    for gid in missing_ids:
        old_fam = old_df.loc[ old_df['v3_id'] == gid, 'family'].values[0]
        fout.write( f"{gid} old family was {old_fam}\n" )
with open("new_genes.txt", "w") as fout:
    fout.write( "{0} new genes were categorized\n".format(len(new_ids)) )
    fout.write( "these are genes that were missing or considered orphans in old grassius\n\n")
    for gid in new_ids:
        new_fam = df.loc[ df['gene_id'] == gid, 'family'].values[0]
        fout.write( f"{gid} has family {new_fam}\n" )



    
    