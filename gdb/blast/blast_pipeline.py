from .blast_commands import prepare_blast_db,run_blastn
from .util import annotate_blast_result
from ..util import load_gene_annotations
from ..fasta import get_records_for_gene_ids,get_gene_id_from_record
import pandas as pd

def run_blast_and_annotate( 
    needle_dna_fasta, haystack_dna_fasta, haystack_gff, gene_id_subset ):
    """
    run blast and get a dataframe with annotated results
    
    if this task is interrupted, the partially-complete results are returned
    
    Arguments:
    ----------
    needle_dna_fasta -- (str) path to a fasta file containing DNA sequences to search for
    haystack_dna_fasta -- (str) path to a fasta file containing a full genome of DNA (blast database)
    haystack_gff -- (str) path to a gff3 file containing annotations corresponding with the genome
    gene_id_subset -- (list of str) subset of ids to consider from the needle_fasta
    """
    
    # prepare blast database with maize v5 genome (DNA)
    prepare_blast_db(haystack_dna_fasta)


    # load maize v5 gene annotations
    gff_data = load_gene_annotations(haystack_gff)


    # prepare final results dataframe
    result_df = pd.DataFrame()

    # iterate over protein sequences related to the list of ids
    for r in get_records_for_gene_ids(needle_dna_fasta,gene_id_subset):
        try:
            print(r.id)

            # blast each protein sequence against maize v5 genome dna
            blast_result = run_blastn( str(r.seq), haystack_dna_fasta )

            # filter blast results
            df = blast_result.data
            blast_result.data = pd.DataFrame(df[df["Identities%"] >= 90])

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
                
        except (KeyboardInterrupt,Exception) as e:
            print(e)
            return result_df
            
    return result_df