# this file contains logic for building database tables that 
# are necessary for the grassius website but are not part of 
# the generic chado schema.

# These tables mainly contain redundant data for the purposes 
# of fast performance

# Some of these tables are traditional but not advantagious. 
# They should be replaced by new rows in chado tables
# (feature, featureprop, feature_relationship)

# These functions should only be called in
# ChadoBuilder.build_grassius_tables()


def build_tfome_metadata( cur, old_grassius_tfomes ):
    """
    Build table "tfome_metadata"
    
    If the table already exists it will be replaced.
    
    Arguments:
    ----------
    old_grassius_tfomes -- (DataFrame) output from
                            get_old_grassius_tfomes()
    """
    
    # create empty table
    cur.execute("DROP TABLE IF EXISTS tfome_metadata")
    cur.execute("""
        CREATE TABLE tfome_metadata (
            tfmd_id SERIAL PRIMARY KEY,
            utname               text,
            transcript_number    text,
            vector               text,
            insert_gene_bank_id  text,
            five_prime_name      text,
            five_prime_seq       text,
            five_prime_temp      text,
            three_prime_name     text,
            three_prime_seq      text,
            three_prime_temp     text,
            pcr_condition        text,
            request_info         text,
            notes                text,
            template             text
        )
    """)
    
    # insert data
    df = old_grassius_tfomes
    for row in df.index:
        
        values = df.loc[row,['utname', 
            'transcript_number', 'vector', 'insert_gene_bank_id',
            'five_prime_name', 'five_prime_seq', 'five_prime_temp',
            'three_prime_name', 'three_prime_seq', 'three_prime_temp',
            'pcr_condition', 'request_info', 'notes','template']].values
        
        cur.execute("""
            INSERT INTO tfome_metadata
            (utname,transcript_number,vector,insert_gene_bank_id,
            five_prime_name,five_prime_seq,five_prime_temp,
            three_prime_name,three_prime_seq,three_prime_temp,
            pcr_condition,request_info,notes,template)
            VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)
        """, tuple(values) )
        


def build_gene_clone( cur, old_grassius_tfomes ):
    """
    Build table "gene_clone"
    
    this is a traditional table who's data will
    correspond with some rows in the chado table
    "feature_relationship"
    
    If the table already exists it will be replaced.
    
    Arguments:
    ----------
    old_grassius_tfomes -- (Dataframe) output from 
                            get_old_grassius_tfomes()
    """
    
    # create empty table
    cur.execute("DROP TABLE IF EXISTS gene_clone")
    cur.execute("""
        CREATE TABLE gene_clone (
            gc_id SERIAL PRIMARY KEY,
            v3_id text,
            clone_name text
        )
    """)

    # insert data
    df = old_grassius_tfomes
    all_names = set()
    for row in df.index:
        gid,utn = df.loc[row,['gene_id','utname']]
        cur.execute("""
            INSERT INTO gene_clone 
            (v3_id,clone_name)
            VALUES (%s,%s)
        """, (gid,utn) )
    

        
def build_gene_interaction( cur ):
    """
    Build table "gene_interaction"
    
    this is a placeholder to make the website work.
    
    TODO: insert gene interactions
    
    If the table already exists it will be replaced
    """
    
    
    # create empty table
    cur.execute("DROP TABLE IF EXISTS gene_interaction")
    cur.execute("""
        CREATE TABLE gene_interaction (
            gi_id SERIAL PRIMARY KEY,
            gene_id text,
            target_id text,
            pubmed_id bigint,
            interaction_type text,
            experiment text,
            protein_name text,
            target_name text
        )
    """)
    
    


def build_seq_features( cur ):
    """
    Build table "seq_features"
    
    this is a placeholder to make the website work.
    
    TODO: insert seq features
    TODO: switch to using chado featureprop
    
    If the table already exists it will be replaced
    """
    
    # create empty table
    cur.execute("DROP TABLE IF EXISTS seq_features")
    cur.execute("""
        CREATE TABLE seq_features (
            sf_id SERIAL PRIMARY KEY,
            feature_id bigint,
            secondary_structures text
        )
    """)
    
    

def build_uniprot_ids( cur ):
    """
    Build table "uniprot_ids".
    
    This is a placeholder to make the website work.
    
    TODO: insert uniprot IDs
    TODO: switch to using chado featureprop
    
    If the table already exists it will be replaced
    """
    
    # create empty table
    cur.execute("DROP TABLE IF EXISTS uniprot_ids")
    cur.execute("""
        CREATE TABLE uniprot_ids (
            ui_id SERIAL PRIMARY KEY,
            gene_name text,
            uniprot_id text
        )
    """)
    

def build_searchable_clones( cur, metadata_df, old_grassius_tfomes ):
    """
    Build table "searchable_clones", based on the 
    current contents of the database

    This table is for fast performance when listing proteins 
    and searching by clone name.

    If the table already exists it will be replaced

    Arguments:
    ----------
    metadata_df -- (DataFrame) a dataframe with columns:
                        "gene_id","name","class","family"
    old_grassius_tfomes -- (dataframe) output from 
                            get_old_grassius_tfomes()
    """

    # create empty table
    cur.execute("DROP TABLE IF EXISTS searchable_clones")
    cur.execute("""
        CREATE TABLE searchable_clones (
            cs_id SERIAL PRIMARY KEY,
            name text,
            clone_list text
        )
    """)

    # insert data
    df = metadata_df
    all_gids = set(old_grassius_tfomes['gene_id'])
    for gid in all_gids:
        
        df_sub = df[df['gene_id'] == gid]
        if len(df_sub.index) == 0:
            continue
        
        name = df_sub['name'].values[0]
        
        clone_list = ' '.join( 
            old_grassius_tfomes.loc[
                old_grassius_tfomes['gene_id'] == gid,
                'utname'].values)
        
        cur.execute("""
            INSERT INTO searchable_clones 
            (name,clone_list)
            VALUES (%s,%s)
        """, (name,clone_list) )


                
        
        
def build_comment_system_urls( cur, all_family_names ):
    """
    Build table "comment_system_urls"

    This table is a placeholder to make the website work.

    This is needed because of the comment system which was
    added in Feb 2022, but is no longer used.

    If the table already exists it will be replaced

    Arguments:
    ----------
    all_family_names -- (list of str)
    """

    # create empty table
    cur.execute("DROP TABLE IF EXISTS comment_system_urls")
    cur.execute("""
        CREATE TABLE comment_system_urls (
            csu_id SERIAL PRIMARY KEY,
            type text,
            name text,
            url text
        )
    """)

    # insert one row for each family
    for family in all_family_names:
        cur.execute("""
            INSERT INTO comment_system_urls 
            (type,name,url)
            VALUES ('family',%s,'')
        """, (family,))
                    
        
        
def build_default_maize_names( cur, metadata_df, gene_versions, all_family_names ):
    """
    Build table "default_maize_names".

    For consistency and fast performance when listing proteins,
    this table contains a default gene id for each protein/version.

    For fast performance when searching for gene IDs, this table 
    contains the full set of gene IDs for each protein, concatenated
    into one string.

    If the table already exists it will be replaced

    Arguments:
    ----------
    metadata_df -- (DataFrame) a dataframe with columns:
                        "gene_id","name","class","family"
    gene_versions -- (dictionary) where keys are gene_ids,
                        values are genome versions e.g. 'v3'
    all_family_names -- (list of str) list of distinct family names
    """

    sorted_families = sorted(list(all_family_names))
        
    # create empty table
    cur.execute("DROP TABLE IF EXISTS default_maize_names")
    cur.execute("""
        CREATE TABLE default_maize_names (
            dmn_id SERIAL PRIMARY KEY,
            name text,
            name_sort_order int,
            family text,
            v3_id text,
            v4_id text,
            v5_id text,
            all_ids text
        )
    """)

    # insert one row for each distinct protein name
    df = metadata_df
    all_names = set(df["name"])
    for name in all_names:
        df_sub = df[df["name"] == name]
        all_gene_ids = df_sub["gene_id"].values
        family = df_sub["family"].values[0]
        concat_ids = " ".join(all_gene_ids)

        # pick default gene IDs (whichever appears first in metadata)
        did = {}
        for v in ['v3','v4','v5']:
            did[v] = next( (gid 
                            for gid in all_gene_ids 
                            if gene_versions[gid] == v)
                          ,"")

        # pick name_sort_order
        # (hidden number used for sorting by protein name)
        family_sort_val = sorted_families.index(family)
        name_sort_val = df_sub["suffix"].values[0]
        nso = int(family_sort_val*10000 + name_sort_val)
            
        # insert row
        cur.execute("""
            INSERT INTO default_maize_names 
            (name,name_sort_order,family,v3_id,v4_id,v5_id,all_ids)
            VALUES (%s,%s,%s,%s,%s,%s,%s)
        """, (name,nso,family,did["v3"],did["v4"],did["v5"],concat_ids))
                
                
                
def build_gene_name( cur, metadata_df, old_grassius_names ):
    """
    Build table "gene_name".
    
    If the table already exists it will be replaced.

    This is a placeholder to make the website work.
    
    Arguments:
    ----------
    metadata_df -- (DataFrame) a dataframe with columns:
                        "class","family"
    old_grassius_names -- (DataFrame) names from old grassius website
                          output from get_old_grassius_names()
    """

        
    # create empty table
    cur.execute("DROP TABLE IF EXISTS gene_name")
    cur.execute("""
        CREATE TABLE gene_name (
            gn_id SERIAL PRIMARY KEY,
            grassius_name text,
            accepted text,
            synonym text
        )
    """)

    # insert one row for each protein name
    all_names = set(metadata_df["name"])
    for name in all_names:
        accepted = 'no'
        synonym = ''
        
        # check if name was in old grassius
        old_match = old_grassius_names[old_grassius_names['name']==name]
        if len(old_match.index) > 0:
            accepted = old_match['accepted'].values[0]
            synonym = old_match['synonym'].values[0]
            
        # check if any related gene_ids were in old grassius
        else:
            all_gene_ids = metadata_df.loc[metadata_df['name']==name,'gene_id'].values
            old_names = old_grassius_names.loc[
                            old_grassius_names['v3_id'].isin(all_gene_ids),
                            'name'].values
            accepted = 'no'
            synonym = ' '.join(old_names)
            
        cur.execute("""
            INSERT INTO gene_name 
            (grassius_name,accepted,synonym)
            VALUES (%s,%s,%s)
        """, (name,accepted,synonym))
    
    
    
        
def build_family_tables( cur, metadata_df, family_desc_df ):
    """
    Build tables "family" and "class_family".

    If either table already exists it will be replaced.

    Arguments:
    ----------
    metadata_df -- (DataFrame) a dataframe with columns:
                        "class","family"
    family_desc_df -- (DataFrame) a dataframe loaded from private 
                        input "family_descriptions"
    """


    # create empty tables
    cur.execute("DROP TABLE IF EXISTS family")
    cur.execute("""
        CREATE TABLE family (
            family_id SERIAL PRIMARY KEY,
            familyname text,
            abbr text,
            description text
        )
    """)
    cur.execute("DROP TABLE IF EXISTS class_family")
    cur.execute("""
        CREATE TABLE class_family (
            cf_id SERIAL PRIMARY KEY,
            class text,
            family text
        )
    """)

    # insert values
    # iterate through metadata rows
    inserted_families = []
    for row in metadata_df.index:

        clazz,family = metadata_df.loc[row,["class","family"]]
        if family in inserted_families:
            continue

        cur.execute("""
            INSERT INTO class_family
            (class, family)
            VALUES (%s,%s)
        """, (clazz,family))

        # attempt to find corresponding row in family descriptions
        match_rows = family_desc_df[family_desc_df["familyname"]==family].index
        if len(match_rows) == 0:
            print( f'WARNING: family "{family}" is not in family descriptions' )
            abbr,desc = "",""
        else:
            abbr,desc = family_desc_df.loc[match_rows[0],["abbr","description"]]
        cur.execute("""
            INSERT INTO family 
            (familyname, abbr, description)
            VALUES (%s,%s,%s)
        """, (family,abbr,desc))

        # make sure this family will not be inserted again
        inserted_families.append(family)

