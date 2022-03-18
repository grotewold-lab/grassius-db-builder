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
            dp = record.description.split()
            gene_parts = [p for p in dp if p.startswith("gene:")]
            if len(gene_parts) == 0:
                continue
            gene = gene_parts[0].split(":")[1]
            if gene in gene_ids:
                yield record
    

def get_gene_id_from_record( record ):
    """
    Get a gene id from the given record, from a fasta file with gene annotations
    The record should be an object of type SeqRecord
    """
    try:
        desc = record.description
        parts = desc.split()
        gene_part = next(p for p in parts if p.startswith("gene:"))
        return gene_part[5:]
    except:
        return ""
    