
# insert 2023 gene interaction data 
# - includes both old data and new
# - includes distances
# - includes notes
# - discussed with john gray summer 2023

import tarfile
import psycopg2
import tempfile
import os


from gdb import InputManager


def apply_patch(cb):
    
    path = InputManager()['2023_gene_interactions']
    # ~3m lines in sql file
    
    sql = ""
    i = 0
            
    temp_dir = tempfile.TemporaryDirectory()
    print()
    # use temp_dir, and when done:
    temp_dir.cleanup()
    with tarfile.open(path, 'r') as tar:
        tar.extractall(temp_dir.name)
        
    path = os.path.join(temp_dir.name,'20231221-gene-interaction.sql')
    with open(path,'r') as fin:
        
        while True:
            line = fin.readline()
            if not line: 
                break
            i += 1
            line = line.strip()
            if line == "":
                continue
            sql += line + "\n"
            if line.startswith("INSERT"):
                break

        with psycopg2.connect(cb.conn_str) as conn:
            with conn.cursor() as cur:
                cur.execute(sql)

                while True:
                    line = fin.readline()
                    if not line: 
                        break
                    i += 1
                    if (i%1000) == 0:
                        print(f"gdb/patches/20231221-gene-interaction.sql line {i}")
                    line = line.strip()
                    if line == "":
                        continue
                    if line.startswith("--"):
                        continue
                    cur.execute(line)
                    
                    
    temp_dir.cleanup()




            
            
            
    
    
    
    
    
        




