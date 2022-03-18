# this file contains a list of boilerplate cvterms and related entries that must be present in the database



# "cvterm" is a table in the chado schema
# required cvterm rows (cvterm_id,cv_id,dbxref_id,name,definition)
required_cvterms = [
    [844,11,877,  "DNA", 
     "An attribute describing a sequence consisting of nucleobases bound to a repeating unit made of a 2-deoxy-D-ribose ring connected to a phosphate backbone. [RSC:cb]"],
    [534,11,562, "polypeptide", 
     "A sequence of amino acids linked by peptide bonds which may lack appreciable tertiary structure and may not be liable to irreversible denaturation. [SO:ma]"],
]
    
    
    
# "cv" is a table in the chado schema
# required cv rows ( cv_id, name, definition )
required_cvs = [
    [11,"sequence","The Sequence Ontology"],
]


# "dbxref" is a table in the chado schema
# required dbxref rows (dbxref_id,db_id,accession)
required_dbxrefs = [
    [562, 9, "0000104"],
    [877, 9, "0000352"]
]