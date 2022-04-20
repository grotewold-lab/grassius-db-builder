from Bio import SeqIO

def get_transcript_gene_dict( fasta_filepath ):
    """
    Extract gene IDs and transcript IDs from the annotations in the given fasta file
    
    Use the output of this function function as an argument for gdb.hmmer.categorize_all_genes
    
    return a dictionary where keys are transcript IDs, and values are gene IDs
    """
    result = {}
    with open(fasta_filepath) as fin:
        for r in SeqIO.parse(fin, "fasta"):
            gene_id = get_gene_id_from_record(r)
            transcript_id = r.name
            result[transcript_id] = gene_id
    return result


def get_gene_transcript_dict( fasta_filepath ):
    """
    Extract gene IDs and transcript IDs from the annotations in the given fasta file
    
    return a dictionary where keys are gene IDs, and values are lists of transcript IDs
    """
    result = {}
    with open(fasta_filepath) as fin:
        for r in SeqIO.parse(fin, "fasta"):
            gene_id = get_gene_id_from_record(r)
            transcript_id = r.name
            if gene_id not in result.keys():
                result[gene_id] = []
            result[gene_id].append( transcript_id )
    return result
    
    
def get_all_gene_ids( fasta_filepath ):
    """
    get a set of gene IDs present in annotations in the given fasta file
    """
    result = set()
    with open(fasta_filepath) as fin:
        for r in SeqIO.parse(fin, "fasta"):
            result.add( get_gene_id_from_record(r) )
    return result
    
    
def read_records_for_gene_ids( fasta_filepath, gene_ids ):
    """
    load a subset of the given fasta file
    yields records that are related to the given list of gene IDs
    """
    
    with open(fasta_filepath) as fin:
        for record in SeqIO.parse(fin, "fasta"):
            gene = get_gene_id_from_record(record)
            if gene in gene_ids:
                yield record
            else:
                #print( f'skipping fasta record with gene id "{gene}"' )
                pass
    
def get_gene_id_from_record( record ):
    """
    Get a gene id from the given record, from a fasta file with gene annotations
    The record should be an object of type SeqRecord
    """
    
    # special case for maize v5 records
    if record.id.startswith("Zm00001eb"):
        return record.id.split("_")[0]
    
    return get_value_from_record_description( record, 'gene' )

    
def get_related_tid_from_record( record ):
    """
    For a fasta record containing a protein sequence, get the transcript id
    of the dna sequence that it derives from.
    
    Typically this information is found in the record's description.
    """
    
    # special case fo maize v5 records
    if record.id.startswith("Zm00001eb"):
        return record.id.replace('_P','_T')
    
    return get_value_from_record_description( record, 'transcript' )
        
    
def get_value_from_record_description( record, key ):
    """
    Extract a string from the given record's description
    
    return None if the given key is not present
    
    Arguments:
    ----------
    record -- (Bio.SeqRecord.SeqRecord) one fasta entry
    key -- (str) the label for the description element in question
            e.g. "gene"
    """
    
    search_str = key + ":"
    
    try:
        desc = record.description
        parts = desc.split()
        target_part = next(p for p in parts if p.startswith(search_str))
        return target_part[len(search_str):]
    except:
        return None
    