
import pandas as pd
import gzip
import json
import pickle
import requests
import psycopg2
import os
from Bio import SeqIO
    
from gdb import InputManager
       

def parse(name):
    i = len(name)-1
    while name[i] in '0123456789':
        i -= 1
    return name,name[:i+1],int(name[i+1:])




def apply_patch(cb):
    """
    create/replace the table 'transcript_domains'
    
    this table is used to boost performance of the custom family query
    (query transcripts based on required/forbidden domains)
    """ 
    
    # load subgenome lists
    folder = os.path.dirname(os.path.abspath(__file__))
    all_fname_formats = [
        'Maize_subgenome{0}_RefGen_AGPv4_final.txt',
        'Maize_subgenome{0}_RefGen_V3_final.txt',
    ]
    all_gids_by_sg = {}
    for sg in [1,2]:
        all_gids = []
        for fname_format in all_fname_formats:
            fname = fname_format.format(sg)
            path = os.path.join(folder,fname)
            count = 0
            with open(path) as fin:
                while True:
                    line = fin.readline()
                    if not line:
                        break
                    line = line.strip()
                    if line == '':
                        continue
                    all_gids.append(line.strip())
                    count += 1
            print( f"loaded {count} gene ids from {fname}" )
        all_gids_by_sg[sg] = all_gids
    
    
    
    with psycopg2.connect(cb.conn_str) as conn:
        with conn.cursor() as cur:

            # add new column
            cur.execute("""\
                drop table if exists subgenome;
                create table subgenome( 
                    sgid SERIAL PRIMARY KEY,
                    geneid TEXT,
                    subgenome TEXT
                );
                create index idx_subgenome_geneid on subgenome(geneid);
                create index idx_subgenomd_subenome on subgenome(subgenome);
            """)
            
            # insert data for new genes
            for sg,all_gids in all_gids_by_sg.items():
                batch_size = 100
                i = 0
                while i < len(all_gids):
                    subset = all_gids[i:i+batch_size]
                    svals = ""
                    for gid in subset:
                        svals += (f"('{gid}','{sg}'),")
                    svals = svals[:-1]
                    cur.execute(f"""
                        INSERT INTO subgenome (geneid,subgenome) VALUES
                            {svals};
                    """)
                    i += batch_size
            
            
            
            
            
            
            
    
    
    
    
    
        




