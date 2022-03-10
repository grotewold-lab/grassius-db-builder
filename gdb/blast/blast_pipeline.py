from .blast_commands import prepare_blast_db,run_tblastn

def run_blast_and_annotate( 
    needle_prot_fasta, haystack_dna_fasta, haystack_gff, gene_id_subset=None ):
    """
    run blast and get a dataframe with annotated results
    
    
    Arguments:
    ----------
    needle_prot_fasta -- (str) path to a fasta file containing sequences to search for
    haystack_dna_fasta -- (str) path to a fasta file containing a full genome (blast database)
    haystack_gff -- (str) path to a gff3 file containing annotations corresponding with the genome
    gene_id_subset -- (optional) (list of str) subset of ids to consider from the needle_fasta
    """
    
    # prepare blast database with maize v5 genome (DNA)
    prepare_blast_db(haystack_dna_fasta)


    # load maize v5 gene annotations
    gff_data = gdb.load_gene_annotations(haystack_gff)


    # prepare final results dataframe
    result_df = pd.DataFrame()

    # iterate over protein sequences related to the list of ids
    for r in get_records_for_gene_ids(needle_prot_fasta,v3_ids):
        print(r.id)

        # blast each protein sequence against maize v5 genome dna
        blast_result = run_tblastn( str(r.seq), haystack_dna_fasta )

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