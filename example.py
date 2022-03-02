# standard imports
import pandas as pd

# local imports
import gdb
from gdb.blast import prepare_blast_db,run_tblastn,annotate_blast_result
from gdb.fasta import get_records_for_gene_ids


# download and/or check integrity of all inputs
im = gdb.InputManager()
im.prepare_all_inputs()


# prepare blast database with maize v5 genome (DNA)
genome_path = im.get_input_filepath("maize_v5_full_fasta")
prepare_blast_db(genome_path)


# load maize v5 gene annotations
gff_path = im.get_input_flepath("maize_v5_gff3")
gff_data = gdb.load_gene_annotations(gff_path)


# load a list of v3 gene IDs from old grassius
filepath = im.get_input_filepath("old_grassius_names")
v3_ids = list(pd.read_excel(filepath)["v3_id"])


# iterate over protein sequences related to the list of ids
filepath = im.get_input_filepath("maize_v3_proteins")
for r in get_records_for_gene_ids(filepath,v3_ids):
    
    # blast each protein sequence against maize v5 genome dna
    blast_result = run_tblastn( str(r.seq), genome_path )
    
    # filter blast results
    df = blast_result.data
    blast_result.data = df[df["Identities%"] >= 99]
    
    # annotate blast results
    annotate_blast_result( blast_result, gff_data )
    
    raise Exception("test")
