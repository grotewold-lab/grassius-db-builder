# this file contains a list of boilerplate cvterms and related entries that must be present in the database
# this file also contains logic to insert and validate those entries



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
    [327, 11,355,"derives_from",""],
    [435, 11,463,"clone","A piece of DNA that has been inserted in a vector so that it can be propagated in a host bacterium or some other organism. [SO:ke]"]
]


# "dbxref" is a table in the chado schema
# required dbxref rows (dbxref_id,db_id,accession)
required_dbxrefs = [
    [562,  9, "0000104"],
    [877,  9, "0000352"],
    [14,   5, "0000002"],
    [524,  9, "0000085"],
    [1414, 9, "0000908"],
    [463,  9, "0000151"],
    [355,  9, "derives_from"],
]
    
    
    
# "cv" is a table in the chado schema
# required cv rows ( cv_id, name, definition )
required_cvs = [
    [11,"sequence","The Sequence Ontology"],
    [8, "taxonomic_rank","A vocabulary of taxonomic ranks (species, family, phylum, etc)"]
]



def init_dbxrefs(cur):
    """
    Make sure the dbxref table has all necessary entries

    this is used in ChadoBuilder constructor
    """

    cur.execute("SELECT dbxref_id from dbxref")
    existing_dbxref_ids = [v[0] for v in cur.fetchall()]

    for entry in required_dbxrefs:
        dbxref_id = entry[0]
        if dbxref_id not in existing_dbxref_ids:
            cur.execute( """
                INSERT INTO dbxref 
                (dbxref_id,db_id,accession) 
                VALUES (%s,%s,%s) 
                """, entry )


def init_cvs(cur):
    """
    Make sure the cv table has all necessary entries

    this is used in ChadoBuilder constructor
    """

    cur.execute("SELECT cv_id from cv")
    existing_cv_ids = [v[0] for v in cur.fetchall()]

    for entry in required_cvs:
        cv_id = entry[0]
        if cv_id not in existing_cv_ids:
            cur.execute( """
                INSERT INTO cv 
                (cv_id,name,definition) 
                VALUES (%s,%s,%s) 
                """, entry )


def init_cvterms(cur):
    """
    Make sure the cvterm table has all necessary entries

    build a dictionary (self.cvterms) to efficiently lookup cvterm IDs

    this is used in ChadoBuilder constructor
    """

    cur.execute("SELECT cvterm_id from cvterm")
    existing_cvterm_ids = [v[0] for v in cur.fetchall()]

    cvterms = {}

    for entry in required_cvterms:
        cvterm_id = entry[0]
        if cvterm_id not in existing_cvterm_ids:
            cur.execute( """
                INSERT INTO cvterm 
                (cvterm_id,cv_id,dbxref_id,name,definition) 
                VALUES (%s,%s,%s,%s,%s) 
                """, entry )

        name = entry[3]
        cvterms[name] = cvterm_id

    return cvterms