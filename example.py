# standard imports
import pandas as pd

# local imports
import gdb
from gdb.blast import run_tblastn
from gdb.fasta import get_records_for_gene_ids


# download and/or check integrity of all inputs
im = gdb.InputManager()
im.prepare_all_inputs()


# load a fasta file with maize v5 full genome DNA
genome_path = im.get_input_filepath("maize_v5_full_fasta")

# load a list of v3 gene IDs from old grassius
filepath = im.get_input_filepath("old_grassius_names")
v3_ids = list(pd.read_excel(filepath)["v3_id"])

# get all protein sequences related to the list of ids
# blast each protein sequence against maize v5 genome dna
filepath = im.get_input_filepath("maize_v3_proteins")
for r in get_records_for_gene_ids(filepath,v3_ids):
    blast_result = run_tblastn( str(r.seq), genome_path )
    raise Exception("test")
