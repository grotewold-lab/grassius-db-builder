
# local imports
import gdb
from gdb.fasta import *
from gdb.hmmer import *


import tempfile
import shutil



"""
judge the quality of the HMMs and thresholds

attempt to re-create old grassius family categories
"""


im = gdb.InputManager()


# create temp folder
folder = tempfile.mkdtemp()
min_pfam_hmm_path = folder + "/pfam_min.hmm"
combined_hmm_path = folder + "/combined.hmm"

# load family criteria
family_criteria = read_family_criteria(im["family_rules"])

# create minified hmm with all necessary accessions
necessary_accessions = get_relevant_accessions(family_criteria)
build_minified_hmm( im["pfam_hmm"], necessary_accessions, min_pfam_hmm_path )
concatenate_hmms( [min_pfam_hmm_path,im["selfbuild_hmm"]], combined_hmm_path )
accessions = get_accessions( combined_hmm_path )
missing_accessions = set(necessary_accessions) - set(accessions)
if len(missing_accessions) > 0:
    raise Exception( "combined hmm file is missing necessary accessions!" )
    

    
# run hmmscan against maize v3 proteins
run_hmmpress( combined_hmm_path )
hmmscan_result = run_hmmscan( combined_hmm_path, im["maize_v3_proteins"] )

# assign family names using UNFILTERED hmmscan results
# check for false negatives
# which CANNOT be fixed using score thresholds
acc_dict = get_acc_dict( hmmscan_result )
families = categorize_all_transcripts( acc_dict, family_criteria )

# assign family names using FILTERED hmmscan results
# check for false positives
# which may be fixable by adjusting score thresholds
filtered_hmmscan_result = get_filtered_hmmscan_result( hmmscan_result )
filtered_acc_dict = get_acc_dict( filtered_hmmscan_result )
filtered_families = categorize_all_transcripts( filtered_acc_dict, family_criteria )


# remove temp folder
shutil.rmtree(folder)
    
    