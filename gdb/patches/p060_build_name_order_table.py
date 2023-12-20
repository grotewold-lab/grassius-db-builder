
import pandas as pd
import gzip
import json
import pickle
import requests
import psycopg2
from Bio import SeqIO
    
from gdb import InputManager
       

def parse(name):
    print(f"parsing name '{name}'")

    i = len(name)-1
    while name[i] in '0123456789':
        i -= 1
    try:
        return name,name[:i+1],int(name[i+1:])
    except ValueError as e:
        return name,name,0



def apply_patch(cb):
    """
    create/replace the table 'transcript_domains'
    
    this table is used to boost performance of the custom family query
    (query transcripts based on required/forbidden domains)
    """ 
    
    
    with psycopg2.connect(cb.conn_str) as conn:
        with conn.cursor() as cur:
        
            cur.execute("""
                select distinct(f.name)
                from feature f
                join featureprop fp
                    on fp.feature_id = f.feature_id
                    and fp.type_id = 1362;
            """)
            all_names = [r[0] for r in cur.fetchall()]
            
            all_parsed_names = []    
            all_prefixes = set()
            for name in all_names:
                if name == 'NaN':
                    continue
                name,prefix,number = parse(name)
                all_prefixes.add(prefix)
                all_parsed_names.append([name,prefix,number])
            
            all_prefixes = sorted(list(all_prefixes))
            
            svals = ''
            for name,prefix,number in all_parsed_names:
                number += 10000 * all_prefixes.index(prefix)
                svals += f"('{name}',{number}),"
            svals = svals[:-1]
                
            cur.execute(f"""
                DROP TABLE IF EXISTS name_orders;
                CREATE TABLE name_orders (
                    noid SERIAL PRIMARY KEY,
                    name text,
                    sortorder integer
                );
                CREATE INDEX name_orders_name ON name_orders(name);
                CREATE INDEX name_orders_sortorder ON name_orders(sortorder);
                INSERT INTO name_orders
                    (name,sortorder)
                    VALUES {svals};
            """)



