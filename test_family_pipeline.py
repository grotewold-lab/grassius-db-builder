
"""
judge the quality of the HMMs and thresholds

attempt to re-create old grassius family categories
"""


# local imports
import gdb
from gdb.fasta import *
from gdb.hmmer import *


import tempfile
import shutil
import os
import pickle
import pandas as pd





im = gdb.InputManager()


# load family criteria
family_criteria_df = read_family_criteria(im["family_rules"])

# load gene -> transcripts dictionary (in two different forms for convenience)
transcript_gene_dict = get_transcript_gene_dict(im["maize_v3_proteins"])
gene_transcript_dict = get_gene_transcript_dict(im["maize_v3_proteins"])

# load old grassius family -> gene dictionary
df = pd.read_excel( im["old_grassius_names"] )
old_grassius_families = {}
for row in df.index:
    family,gene = df.loc[row,["family","v3_id"]]
    if family not in old_grassius_families.keys():
        old_grassius_families[family] = []
    old_grassius_families[family].append(gene)


# attempt to load pickled hmmscan results and skip some steps
premade_result_path = "hmmscan_result.p"
if os.path.exists(premade_result_path):
    print( "loading premade results and skipping steps..." )
    
    with open( premade_result_path, "rb" ) as fin:
        hmmscan_result = pickle.load( fin )
    
    made_temp_folder = False
    
else:
    
    # create temp folder
    folder = tempfile.mkdtemp()
    made_temp_folder = True
    min_pfam_hmm_path = folder + "/pfam_min.hmm"
    combined_hmm_path = folder + "/combined.hmm"

    # create minified hmm with all necessary accessions
    necessary_accessions = get_relevant_accessions(family_criteria_df)
    build_minified_hmm( im["pfam_hmm"], necessary_accessions, min_pfam_hmm_path )
    concatenate_hmms( [min_pfam_hmm_path,im["selfbuild_hmm"]], combined_hmm_path )
    accessions = get_accessions( combined_hmm_path )
    missing_accessions = set(necessary_accessions) - set(accessions)
    if len(missing_accessions) > 0:
        raise Exception( "combined hmm file is missing necessary accessions!" )


    # run hmmscan against maize v3 proteins
    run_hmmpress( combined_hmm_path )
    hmmscan_result = run_hmmscan( combined_hmm_path, im["maize_v3_proteins"] )


df = hmmscan_result.data

# check if a false-negative gene id is un-fixable using score thresholds
def is_problematic( gid, forbidden ):
    for tid in gene_transcript_dict[gid]:
        accs = df.loc[df["query name"] == tid,"accession"].values
        for acc in accs:
            if acc in forbidden:
                # not problematic because a stricter threshold might fix the issue
                return False 
    return True
    
    
# assign family names using UNFILTERED hmmscan results
# check for false negatives
# which CANNOT be fixed using score tdhresholds
acc_dict = get_acc_dict( hmmscan_result )
families = categorize_all_genes( acc_dict, family_criteria_df, transcript_gene_dict )

for family_name in families.keys():
    if family_name not in old_grassius_families.keys():
        continue
    required = family_criteria_df.loc[family_criteria_df["GRASSIUS"] == family_name,"Required"].values[0]
    forbidden = family_criteria_df.loc[family_criteria_df["GRASSIUS"] == family_name,"Forbidden"].values[0]
    fn_genes = set(old_grassius_families[family_name]) - set(families[family_name])
    fn_genes = [g for g in fn_genes if g in gene_transcript_dict.keys()]
    fn_genes = [g for g in fn_genes if is_problematic(g, forbidden)]
    
    if len(fn_genes) > 0:
        
        print( f"\nrequired for {family_name} family: {required}" )
        print( f"forbidden for {family_name} family: {forbidden}" )
        print( f"problematic false-negatives for {family_name} family:" )
        for gid in fn_genes:
            print( f"\t{gid}" )
            for tid in gene_transcript_dict[gid]:
                accs = df.loc[df["query name"] == tid,"accession"].values
                print( f"\t\t{tid} has accessions: {accs}" )

                

raise Exception("test")


# assign family names using FILTERED hmmscan results

filtered_hmmscan_result = get_filtered_hmmscan_result( hmmscan_result )
filtered_acc_dict = get_acc_dict( filtered_hmmscan_result )
filtered_families = categorize_all_genes( filtered_acc_dict, family_criteria_df, transcript_gene_dict )


# remove temp folder if necessary
if made_temp_folder:
    shutil.rmtree(folder)
    
    