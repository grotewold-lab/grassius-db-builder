
import pandas as pd
import gzip
import json
import pickle
import requests
import psycopg2
from Bio import SeqIO
    
from gdb import InputManager
       
def apply_patch(cb):
    """
    create/replace the table 'transcript_domains'
    
    this table is used to boost performance of the custom family query
    (query transcripts based on required/forbidden domains)
    """ 
    
    
    with psycopg2.connect(cb.conn_str) as conn:
        with conn.cursor() as cur:
            # create/replace table
            cur.execute("""
                DROP TABLE IF EXISTS transcript_domains;
                DROP TABLE IF EXISTS all_domains;
                CREATE TABLE transcript_domains (
                    dd_id SERIAL PRIMARY KEY,
                    tid text,
                    domains text
                );
            """)

            
            # get a list of maize v5 protein transcript IDs
            cur.execute("""
                SELECT f.uniquename FROM feature f
                WHERE f.type_id = 534 AND f.organism_id = 4
            """ )
            all_tids = [x[0] for x in cur.fetchall()]
            
            i = 0
            
            # iterate through protein transcripts
            fasta_filepath = InputManager()['maize_v5_proteins']
            with open(fasta_filepath) as fin:
                for r in SeqIO.parse(fin, "fasta"):
                    
                    i += 1
                    if (i%1000) == 0:
                        print( i )
                        
                    tid = r.id
                    seq = str(r.seq)
                    seqlen = len(seq)
                    

                    # get domain annotation
                    cur.execute("""
                        SELECT fp.value FROM featureprop fp
                        JOIN feature f ON f.feature_id = fp.feature_id
                        WHERE f.uniquename=%s AND fp.type_id = 61467
                    """, (tid,) )
                    anno = cur.fetchall()
                    sanno = "[" + ",".join([x[0] for x in anno]) + "]"
                    anno = json.loads(sanno)


                    # insert sequence length into domain annotations
                    for entry in anno:
                        entry["color"] = "none"
                        entry["seqlen"] = seqlen
                    sanno = json.dumps( anno )

                    # append to sql file
                    cur.execute("""
                        INSERT INTO transcript_domains (tid,domains)
                        VALUES( %s, %s );
                    """, (tid,sanno) )




