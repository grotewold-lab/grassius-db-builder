
import os
import pandas as pd
import gzip
import json
import pickle
import requests
import psycopg2

from gdb.hmmer import get_family_criteria

def apply_patch( cb ):
    """
    create/replace tables 'family_domain_colors' and 'default_domains'
    
    these tables are used to boost performance of domain annotation graphics
    """
    
    #def lookup_tid_of_interest(name):
    #    url = "http://localhost:8080/proteininfor/Maize/" + name
    #    global test
    #    resp = requests.get( url )
    #    src = resp.text
    #    se = "http://maizegdb.org/gene_center/gene/"
    #    i = src.index(se) + len(se)
    #    j = src.index( '"', i )
    #    return src[i:j]
    #
    #tids_of_interest = {}
    #df = pd.read_csv("protein_names.csv")
    #all_names = sorted(list(set(df['name'])))
    #for name in all_names:
    #    if name not in tids_of_interest.keys():
    #        tid = lookup_tid_of_interest(name)
    #        print(name,tid)
    #        tids_of_interest[name] = tid
    #with open('tids_of_interest.pi', 'wb') as f:
    #    pickle.dump( tids_of_interest, f )


    # identify transcript of interest for each protein name
    folder = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(folder,'tids_of_interest.pi')
    with open(path, 'rb') as f:
        tids_of_interest = pickle.load(f)


    #all_anno = get_domain_annotations()
    family_criteria_df = get_family_criteria()

    # build color codes based on family criteria
    all_colors = ["red","green","blue"]
    all_color_codes = {}
    df = family_criteria_df
    for row in df.index:
        family_name = df.loc[row,'GRASSIUS']
        sreq = df.loc[row,'Required']
        domains = sorted(list(set([ s.split("#")[0] for s in sreq.split(":")])))
        if len(domains) > len(all_colors):
            raise Exception('need more colors')
        family_colors = { domains[i]:all_colors[i] for i in range(len(domains)) }
        all_color_codes[family_name] = family_colors


    # connect to the database
    with psycopg2.connect(cb.conn_str) as conn:
        with conn.cursor() as cur:

            # create/replace table
            cur.execute("""
                DROP TABLE IF EXISTS family_domain_colors;
                CREATE TABLE family_domain_colors (
                    fdc_id SERIAL PRIMARY KEY,
                    family text,
                    domain text,
                    color text
                );
            """)

            # insert rows
            for family,family_colors in all_color_codes.items():
                for domain,color in family_colors.items():
                    cur.execute(f"""
                        INSERT INTO family_domain_colors (family,domain,color)
                        VALUES( %s,%s,%s );
                    """, (family,domain,color) )



            # create/replace table
            cur.execute("""
                DROP TABLE IF EXISTS default_domains;
                CREATE TABLE default_domains (
                    dd_id SERIAL PRIMARY KEY,
                    protein_name text,
                    domains text
                );
            """)

            # iterate through protein names
            for name,tid in tids_of_interest.items():
                print(name)
                
                if len(tid) == 0:
                    cur.execute("""
                        INSERT INTO default_domains (protein_name,domains)
                        VALUES( %s,'');
                    """, (name,) )
                    continue

                # get sequence length
                protein_tid = tid.replace('_T','_P')
                cur.execute("""
                    SELECT f.residues FROM feature f
                    WHERE f.uniquename=%s
                """, (protein_tid,) )
                seq = cur.fetchall()[0][0]
                seqlen = len(seq)

                # get family
                cur.execute("""
                    SELECT fp.value FROM featureprop fp
                    JOIN feature f ON f.feature_id = fp.feature_id
                    WHERE f.uniquename=%s AND fp.type_id = 1362
                """, (tid,) )
                family = cur.fetchall()[0][0]
                if family in all_color_codes.keys():
                    color_codes = all_color_codes[family]
                else:
                    color_codes = {}

                # get domain annotation
                cur.execute("""
                    SELECT fp.value FROM featureprop fp
                    JOIN feature f ON f.feature_id = fp.feature_id
                    WHERE f.uniquename=%s AND fp.type_id = 61467
                """, (tid,) )
                anno = cur.fetchall()
                sanno = "[" + ",".join([x[0] for x in anno]) + "]"
                anno = json.loads(sanno)

                # insert color info into domain annotations
                for entry in anno:
                    acc = entry['accession'].split('.')[0]
                    if acc in color_codes.keys():
                        color = color_codes[acc]
                    else:
                        color = "none"
                    entry["color"] = color
                    entry["seqlen"] = seqlen
                sanno = json.dumps( anno )

                # insert new entry into database
                cur.execute("""
                    INSERT INTO default_domains (protein_name,domains)
                    VALUES( %s,%s );
                """, (name,sanno) )
