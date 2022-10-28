
import pandas as pd

from gdb import InputManager
from gdb.chado import *
       
def apply_patch(cb):
    """
    create/replace chado features for non-maize species
    (rice, sorghum, sugarcane, brachypodium)
    """ 
    
    im = InputManager()
    df = pd.read_csv( im['non_maize_metadata'] )
    cb = ChadoBuilder()


    # insert sequences from fasta files
    for organism in ['Rice','Sorghum','Sugarcane','Brachypodium']:
        for suffix in ["cdna","proteins"]:
            fasta_filepath = im[f"{organism}_{suffix}".lower()]
            cb.insert_sequences( organism + "_", df, fasta_filepath, 
                                is_protein=(suffix=='proteins') )
