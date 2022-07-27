
import pandas as pd
import gzip
import json
import pickle
import requests
import psycopg2
        
def apply_patch( cb ):
    """
    create/replace table 'acc_list'
    
    this table is used to boost the performance of the 
    autocomplete text boxes on the "custom family" page
    """
    
    # connect to the database
    with psycopg2.connect(cb.conn_str) as conn:
        with conn.cursor() as cur:
            
            # create/replace table
            cur.execute("""
                DROP TABLE IF EXISTS acc_list;
                CREATE TABLE acc_list (
                    al_id SERIAL PRIMARY KEY,
                    accession text
                );
            """)

            # get all domain annotations
            cur.execute("""
                SELECT fp.value FROM featureprop fp
                JOIN feature f ON f.feature_id = fp.feature_id
                WHERE fp.type_id = 61467
            """ )
            anno = cur.fetchall()
            sanno = "[" + ",".join([x[0] for x in anno]) + "]"
            anno = json.loads(sanno)


            # extract set of unique accession names
            all_accessions = set()
            for entry in anno:
                acc = entry['accession'].split('.')[0]
                print(acc)
                all_accessions.add(acc)

            # insert rows
            for acc in all_accessions:
                cur.execute("""
                    INSERT INTO acc_list (accession)
                    VALUES( %s );
                """, (acc,) )




