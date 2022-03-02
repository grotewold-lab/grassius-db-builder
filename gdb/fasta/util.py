from Bio import SeqIO

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
    