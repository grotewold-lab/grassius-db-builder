# local imports
import gdb
from gdb import InputManager
from gdb.chado import ChadoBuilder


import pandas as pd


def test_chado_builder_constructor():
    cb = ChadoBuilder()

    
def test_chado_builder():
    """
    build a sample database with 3 genes
    """
    
    cb = ChadoBuilder()
    im = InputManager()
    
    metadata_df = pd.DataFrame(data={
        "gene_id": ["GRMZM2G169820","GRMZM2G153233","GRMZM2G078274"],
        "name": ["ZmARF1","ZmARF2","ZmARF3"],
        "class": ["TF"] * 3,
        "family": ["ARF"] * 3
    })
    
    organism = "Maize_v3"
    fasta_filepath = im["maize_v3_cdna"]
    
    cb.insert_sequences( organism, metadata_df, fasta_filepath, is_protein=False )
    
    # query gene IDs
    # expected response: [('GRMZM2G078274',), ('GRMZM2G153233',), ('GRMZM2G169820',)]
    response = cb.query("SELECT DISTINCT value FROM featureprop WHERE type_id=496")
    assert len(response) == 3
    
    
    