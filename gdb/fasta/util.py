from Bio import SeqIO

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
    
    
    
def get_records_for_gene_ids( fasta_filepath, gene_ids ):
    """
    load a subset of the given fasta file
    only return records that are related to the given list of gene IDs
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
    
    
def get_all_gene_ids( fasta_filepath ):
    """
    get a set of gene IDs present in annotations in the given fasta file
    """
    result = set()
    with open(fasta_filepath) as fin:
        for r in SeqIO.parse(fin, "fasta"):
            result.add( get_gene_id_from_record(r) )
    return result
    