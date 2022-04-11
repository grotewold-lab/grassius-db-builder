

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
family_desc_df = pd.read_csv(im['family_descriptions'])
old_grassius_names = get_old_grassius_names()
old_grassius_tfomes = get_old_grassius_tfomes()

# load metadata (created in build_new_grassius_db.py)
df = pd.read_csv('metadata.csv')

    
# build gene_id -> genome_version dictionary
gene_versions = {}
for version in ['v3','v4','v5']:
    all_gene_ids = get_all_gene_ids( im[f'maize_{version}_proteins'] )
    for gid in all_gene_ids:
        gene_versions[gid] = version
    
    
# start building database
cb = ChadoBuilder()
cb.build_grassius_tables( df, gene_versions, family_desc_df, 
                         old_grassius_names, old_grassius_tfomes )

    
# insert sequences from fasta files
#for suffix in ["proteins"]:
#    for version in ["v3","v4","v5"]:
#        organism = f"Maize_{version}"
#        fasta_filepath = im[f"maize_{version}_{suffix}"]
#        cb.insert_sequences( organism, df, fasta_filepath, is_protein=False )
        
cb.insert_tfomes( old_grassius_tfomes )