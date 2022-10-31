# this file contains a list of grassius-specific species labels that must be present in the chado table "organism"


# local imports
from ..grassius import get_species_descriptions
    
    

# required organism rows (abbreviation,genus,species,common_name,infraspecific_name)
required_organisms = [
    ["Z. mays","Zea","mays","Maize","v3"],
    ["Z. mays","Zea","mays","Maize","v4"],
    ["Z. mays","Zea","mays","Maize","v5"],    
    ["O. sativa","Oryza","sativa","Rice",""],
    ["S. bicolor","Sorghum","bicolor","Sorghum",""],
    ["S. officinarum","Saccharum","officinarum","Sugarcane",""],
    ["B. distachyon","Brachypodium","distachyon","Brachypodium",""]
]


def get_all_organisms():
    """
    Get a list of chado organisms in useful readable format
    
    each organism is represented by common name and infraspecific name
    
    return a list of strings e.g. "Maize_v4"
    """
    return ["_".join(r[-2:]) for r in required_organisms]



        
def init_organisms(cur):
    """
    Make sure the organism table has all necessary entries

    build a dictionary (self.organisms) to efficiently lookup organism IDs

    this is used in constructor
    """        

    species_descriptions = get_species_descriptions()
    
    organisms = {}

    for entry in required_organisms:
        c_name, is_name = entry[-2:]

        if c_name not in species_descriptions.keys():
            print(f'WARNING missing description for organism with common_name "{c_name}"')
            description = None
        else:
            description = species_descriptions[c_name]
        
        cur.execute("""
            SELECT organism_id
            FROM organism
            WHERE (common_name = %s) and (infraspecific_name = %s);
        """, (c_name, is_name))
        result = cur.fetchone()
        
        if result is not None:
            org_id = result[0]
        else:
            cur.execute( """
                INSERT INTO organism
                (abbreviation,genus,species,common_name,infraspecific_name,comment) 
                VALUES (%s,%s,%s,%s,%s,%s) 
                RETURNING organism_id
            """, entry + [description] )
            result = cur.fetchone()
            org_id = result[0]

        organisms[c_name+"_"+is_name] = org_id

    # check consistency with get_all_organisms()
    for org in get_all_organisms():
        if org not in organisms.keys():
            raise Exception( f"internal consistency failed. id for organism {org}" )

    return organisms