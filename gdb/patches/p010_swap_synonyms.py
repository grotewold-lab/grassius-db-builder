# swap protein names with synonyms for "AP2/ERF-AP2" family

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
    swap protein names with synonyms for "AP2/ERF-AP2" family
    """
    
    with psycopg2.connect(cb.conn_str) as conn:
        with conn.cursor() as cur:

            cur.execute("""
                SELECT gn_id,grassius_name,synonym FROM gene_name
                WHERE grassius_name LIKE 'ZmERFAP2_%'
                AND synonym LIKE 'ZmEREB%';
            """)
            result = cur.fetchall()

            for gn_id,grassius_name,synonym in result:
                old_synonym = synonym
                new_synonym = grassius_name
                old_protein_name = grassius_name
                new_protein_name = synonym
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



            # update gene_name.name_sort_order
            cur.execute("""
                SELECT gn.gn_id,gn.grassius_name,name_sort_order 
                FROM gene_name gn
                JOIN default_maize_names dmn
                    ON dmn.name = gn.grassius_name
                WHERE dmn.family = 'AP2/ERF-AP2';
            """)
            result = cur.fetchall()
            base_value = min([x[2] for x in result])
            all_prefix_offsets = {
                "ZmEREB": 0,
                "ZmERFAP2_": 300
            }
            for gn_id,name,_ in result:

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
                    UPDATE gene_name
                    SET name_sort_order = %s
                    WHERE gn_id = %s
                """, (new_sort_order,gn_id) )

