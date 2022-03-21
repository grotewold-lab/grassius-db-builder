# this file contains a list of grassius-specific species labels that must be present in the chado table "organism"


    

# required organism rows (abbreviation,genus,species,common_name,infraspecific_name)
required_organisms = [
    ["Z. mays","Zea","mays","Maize","v3"],
    ["Z. mays","Zea","mays","Maize","v4"],
    ["Z. mays","Zea","mays","Maize","v5"],    
]


def get_all_organisms():
    """
    Get a list of chado organisms in useful readable format
    
    each organism is represented by common name and infraspecific name
    
    return a list of strings e.g. "Maize_v4"
    """
    return ["_".join(r[-2:]) for r in required_organisms]