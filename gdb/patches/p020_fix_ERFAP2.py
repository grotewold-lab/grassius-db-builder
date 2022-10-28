# run this after swap_synonyms.py
# switch remaining "ZmERFAP2_..." protein names to "ZmEREB..."


import psycopg2

# local imports
import gdb
from gdb.fasta import *
from gdb.hmmer import *
from gdb.itak import *
from gdb.grassius import *
from gdb.chado import *


def has_dd_table(cb):
    with psycopg2.connect(cb.conn_str) as conn:
        with conn.cursor() as cur:
            try:
                cur.execute("select * from default_domains limit 1")
                return True
            except:
                return False


def apply_patch(cb):
    """
    after swapping synonyms with names for the ERFAP2 family,
    switch remaining "ZmERFAP2_..." protein names to "ZmEREB..."
    """
    
    with psycopg2.connect(cb.conn_str) as conn:
        with conn.cursor() as cur:

            # find next available suffix for ZmEREB...
            cur.execute("""
                SELECT grassius_name FROM gene_name
                WHERE grassius_name LIKE 'ZmEREB%'
            """)
            init_ereb_names = [x[0] for x in cur.fetchall()]
            init_ereb_suffixes = [int(s[6:]) for s in init_ereb_names]
            next_ereb_suffix = max(init_ereb_suffixes)+1

            # find proteins that need to be renamed
            cur.execute("""
                SELECT gn_id,grassius_name FROM gene_name
                WHERE grassius_name LIKE 'ZmERFAP2_%'
            """)
            to_rename = cur.fetchall()

            to_rename = sorted(to_rename, key=lambda x: int(x[1][9:]))

            # perform renaming
            for gn_id,old_protein_name in to_rename:
                new_protein_name = f'ZmEREB{next_ereb_suffix}'
                new_synonym = old_protein_name
                next_ereb_suffix += 1
                print( f"{old_protein_name} -> {new_protein_name}" )

                # ZmERFAP2_21
                # grep -n "FROM stdin\|ZmERFAP2_21" build_db.sql

                # default_domains
                # select * from default_domains where protein_name = 'ZmERFAP2_21'
                if has_dd_table(cb):
                    cur.execute("""
                        UPDATE default_domains
                        SET protein_name = %s
                        WHERE protein_name = %s
                    """, (new_protein_name,old_protein_name) )

                # default_maize_names
                # select * from default_maize_names where name = 'ZmERFAP2_21';
                cur.execute("""
                    UPDATE default_maize_names
                    SET name = %s
                    WHERE name = %s
                """, (new_protein_name,old_protein_name) )


                # feature
                # select * from feature where name='ZmERFAP2_21';
                cur.execute("""
                    UPDATE feature
                    SET name = %s
                    WHERE name = %s
                """, (new_protein_name,old_protein_name) )

                # gene_interaction
                # select * from gene_interaction where target_name = 'ZmERFAP2_21';
                cur.execute("""
                    UPDATE gene_interaction
                    SET target_name = %s
                    WHERE target_name = %s
                """, (new_protein_name,old_protein_name) )


                # gene_name
                # select * from gene_name where grassius_name = 'ZmERFAP2_21';
                cur.execute("""
                    UPDATE gene_name
                    SET grassius_name = %s, synonym = %s
                    WHERE gn_id = %s
                """, (new_protein_name,new_synonym,gn_id) )


                # searchable_clones
                # select * from searchable_clones where name = 'ZmERFAP2_21';
                cur.execute("""
                    UPDATE searchable_clones
                    SET name = %s
                    WHERE name = %s
                """, (new_protein_name,old_protein_name) )


                # uniprot_ids
                # select * from uniprot_ids where gene_name = 'ZmERFAP2_21';
                cur.execute("""
                    UPDATE uniprot_ids
                    SET gene_name = %s
                    WHERE gene_name = %s
                """, (new_protein_name,old_protein_name) )


            # update default_maize_names.name_sort_order
            cur.execute("""
                SELECT dmn_id,name,name_sort_order 
                FROM default_maize_names
                WHERE family = 'AP2/ERF-AP2';
            """)
            result = cur.fetchall()
            base_value = min([x[2] for x in result])
            all_prefix_offsets = {
                "ZmEREB": 0,
                "ZmERFAP2_": 300
            }
            for dmn_id,name,_ in result:

                accepted = False
                for prefix,offset in all_prefix_offsets.items():
                    if name.startswith( prefix ):
                        suffix = int( name.replace(prefix,'') )
                        new_sort_order = base_value + offset + suffix
                        accepted = True
                        break

                if not accepted:
                    raise Exception('unrecognized name: ' + name )

                cur.execute("""
                    UPDATE default_maize_names
                    SET name_sort_order = %s
                    WHERE dmn_id = %s
                """, (new_sort_order,dmn_id) )
