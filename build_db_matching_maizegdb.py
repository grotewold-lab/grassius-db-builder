
# standard imports
import pandas as pd

# local imports
import gdb
from gdb import InputManager
from gdb.chado import ChadoBuilder
from gdb.fasta import get_gene_id_from_record
from gdb.grassius import get_maizegdb_associations


im = InputManager()
#cb = ChadoBuilder()

# read maizegdb gene_id associations  
metadata_df = get_maizegdb_associations()
            
# build dummy metadata 
metadata_df["class"] = ["dummy_class"] * len(metadata_df.index)
metadata_df["family"] = ["dummy_family"] * len(metadata_df.index)

raise Exception("test")
            
    
# insert sequences from fasta files
for version in ["v3","v4","v5"]:
    organism = f"Maize_{version}"
    fasta_filepath = im[f"maize_{version}_cdna"]
    cb.insert_sequences( organism, metadata_df, fasta_filepath, is_protein=False )


# save a snapshot of the database that was built
cb.write_snapshot( "build_db.sql.tar.gz" )