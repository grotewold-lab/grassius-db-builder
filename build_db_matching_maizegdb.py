
# standard imports
import pandas as pd

# local imports
import gdb
from gdb import InputManager
from gdb.chado import ChadoBuilder
from gdb.fasta import get_gene_id_from_record


im = InputManager()
cb = ChadoBuilder()

# read maizegdb gene_id associations    
all_gene_ids = []
all_names = []
line_count = 0
filepath = im["maizegdb_gene_id_associations"]
with open(filepath, "r") as fin:
    
    # skip first line
    fin.readline()
    
    while True:
        
        line = fin.readline()
        if not line:
            break
        
        parts = line.split("\t")
        for gid in parts[1:]:
            gid = gid.strip()
            if len(gid) == 0:
                continue
            elif gid.startswith("B73v1_"):
                continue
            elif gid.startswith("B73v2_"):
                continue
            elif gid.startswith("B73v3_"):
                gid = gid[6:]
            elif gid.startswith("GRMZM"):
                pass
            elif gid.startswith("Zm00001d"):
                pass
            elif gid.startswith("Zm00001eb"):
                pass
            else:
                raise Exception("unrecognized prefix for gene id " + gid)
                
            all_gene_ids.append( gid )
            all_names.append( parts[0] )
            
        
# build dummy metadata 
metadata_df = pd.DataFrame(data={
    "gene_id": all_gene_ids,
    "name": all_names,
    "class": ["dummy_class"] * len(all_gene_ids),
    "family": ["dummy_family"] * len(all_gene_ids)
})
            
    
# insert sequences from fasta files
for version in ["v3","v4","v5"]:
    organism = f"Maize_{version}"
    fasta_filepath = im[f"maize_{version}_cdna"]
    cb.insert_sequences( organism, metadata_df, fasta_filepath, is_protein=False )


# save a snapshot of the database that was built
cb.write_snapshot( "build_db.sql.tar.gz" )