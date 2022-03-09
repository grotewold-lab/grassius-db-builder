# standard imports
import pandas as pd

# local imports
import gdb
from gdb.blast import read_blast_output,prepare_blast_db,run_tblastn,annotate_blast_result
from gdb.fasta import get_gene_id_from_record,get_records_for_gene_ids


# download and/or check integrity of all inputs
im = gdb.InputManager()
im.prepare_all_inputs()


# prepare blast database with maize v5 genome (DNA)
genome_path = im.get_input_filepath("maize_v5_full_fasta")
prepare_blast_db(genome_path)


# load maize v5 gene annotations
gff_path = im.get_input_filepath("maize_v5_gff3")
gff_data = gdb.load_gene_annotations(gff_path)


# load a list of v3 gene IDs from old grassius
filepath = im.get_input_filepath("old_grassius_names")
v3_ids = list(pd.read_excel(filepath)["v3_id"])


# prepare final results dataframe
result_df = pd.DataFrame()

# iterate over protein sequences related to the list of ids
filepath = im.get_input_filepath("maize_v3_proteins")
for r in get_records_for_gene_ids(filepath,v3_ids):
    print(r.id)
    
    # blast each protein sequence against maize v5 genome dna
    blast_result = run_tblastn( str(r.seq), genome_path )
    
    # filter blast results
    df = blast_result.data
    blast_result.data = pd.DataFrame(df[df["Identities%"] >= 99])
    
    # annotate blast results
    annotate_blast_result( blast_result, gff_data )
    
    # append to final results
    df = blast_result.data
    for row in df.index:
        new_row_index = len(result_df.index)
        result_df.loc[new_row_index,"v3_gene_id"] = get_gene_id_from_record(r)
        result_df.loc[new_row_index,"v3_transcript_id"] = r.id
        result_df.loc[new_row_index,"matching_v5_chrom"] = df.loc[row,"Chrom"]
        result_df.loc[new_row_index,"matching_v5_start_pos"] = df.loc[row,"Start_Pos"]
        result_df.loc[new_row_index,"matching_v5_stop_pos"] = df.loc[row,"Stop_Pos"]
        result_df.loc[new_row_index,"matching_v5_gene_id"] = df.loc[row,"Gene ID"]
        result_df.loc[new_row_index,"Identities%"] = df.loc[row,"Identities%"]
        
        
# save final results
result_df.to_csv("report.csv", index=False)
