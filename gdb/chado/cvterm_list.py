# this file contains a list of boilerplate cvterms and related entries that must be present in the database



# "cvterm" is a table in the chado schema
# required cvterm rows (cvterm_id,cv_id,dbxref_id,name,definition)
required_cvterms = [
    [844, 11,877,  "DNA", 
     "An attribute describing a sequence consisting of nucleobases bound to a repeating unit made of a 2-deoxy-D-ribose ring connected to a phosphate backbone. [RSC:cb]"],
    [534, 11,562, "polypeptide", 
     "A sequence of amino acids linked by peptide bonds which may lack appreciable tertiary structure and may not be liable to irreversible denaturation. [SO:ma]"],
    [13,  8, 14,"class",""],
    [496, 11,524,"gene_by_genome_location",""],
    [1362,11,1414,"supported_by_domain_match","An attribute to describe a feature that has been predicted using sequence similarity of a known domain. [SO:ke]"],
    [327,11,355,"derives_from",""],
]


# "dbxref" is a table in the chado schema
# required dbxref rows (dbxref_id,db_id,accession)
required_dbxrefs = [
    [562,  9, "0000104"],
    [877,  9, "0000352"],
    [14,   5, "0000002"],
    [524,  9, "0000085"],
    [1414, 9, "0000908"],
    [355,  9, "derives_from"],
]
    
    
    
# "cv" is a table in the chado schema
# required cv rows ( cv_id, name, definition )
required_cvs = [
    [11,"sequence","The Sequence Ontology"],
    [8, "taxonomic_rank","A vocabulary of taxonomic ranks (species, family, phylum, etc)"]
]