# local imports
from .docker_util import *
from ..fasta import *
from .cvterm_list import required_dbxrefs,required_cvs,required_cvterms
from .chado_organisms import required_organisms

import psycopg2

default_port = 8642


class ChadoBuilder:
    
    def __init__(self, port=None):
        """
        Construct a new instance of ChadoBuilder
        
        Find or create a docker container to host the database while it is being built, 
        and test the connection to the database
        
        There can be at most one such container.  If another instance of ChadoBuilder 
        is created with a different port, the first container will be replaced
            
        For the rest of this object's life:
            It is assumed that the docker container will continue running. 
            Otherwise exceptions will be thrown
    
        Arguments:
        ----------
        port -- (int) (optional) the host port that will be used to access the database
        """
        
        if port is None:
            port = default_port
        
        # find or create docker container
        container = get_existing_db_container()
        if (container is None) or (not container_maps_to_port( container, port )):
            container = init_db_container( port, replace=True )
        else:
            print( "connected to existing database container" )
                    
        # set instance variables
        self.port = port
        self.docker_container = container
        self.conn_str = f"dbname=postgres host=localhost port={port} user=postgres password=postgres"
        
        with psycopg2.connect(self.conn_str) as conn:
            with conn.cursor() as cur:
                
                # make sure we can connect to the database
                self._test_db_connection(cur)

                # check/initialize boilerplate tables and member vars
                self._init_dbxrefs(cur)
                self._init_cvs(cur)
                self._init_cvterms(cur) # self.cvterms
                self._init_organisms(cur) # self.organisms
                
                
    def _test_db_connection(self,cur):
        """
        raise an exception if we can't make a query
        
        this is used in constructor
        """
        try:
            cur.execute("SELECT * from cvterm limit 1")
        except:
            raise Exception("could not connect to the database")
        
        
    def _init_dbxrefs(self, cur):
        """
        Make sure the dbxref table has all necessary entries
        
        this is used in constructor
        """
        
        cur.execute("SELECT dbxref_id from dbxref")
        existing_dbxref_ids = [v[0] for v in cur.fetchall()]

        for entry in required_dbxrefs:
            dbxref_id = entry[0]
            if dbxref_id not in existing_dbxref_ids:
                cur.execute( """
                    INSERT INTO dbxref 
                    (dbxref_id,db_id,accession) 
                    VALUES (%s,%s,%s) 
                    """, entry )
        
        
    def _init_cvs(self, cur):
        """
        Make sure the cv table has all necessary entries
        
        this is used in constructor
        """
        
        cur.execute("SELECT cv_id from cv")
        existing_cv_ids = [v[0] for v in cur.fetchall()]

        for entry in required_cvs:
            cv_id = entry[0]
            if cv_id not in existing_cv_ids:
                cur.execute( """
                    INSERT INTO cv 
                    (cv_id,name,definition) 
                    VALUES (%s,%s,%s) 
                    """, entry )
                        
        
    def _init_cvterms(self, cur):
        """
        Make sure the cvterm table has all necessary entries
        
        build a dictionary (self.cvterms) to efficiently lookup cvterm IDs
        
        this is used in constructor
        """
        
        cur.execute("SELECT cvterm_id from cvterm")
        existing_cvterm_ids = [v[0] for v in cur.fetchall()]

        cvterms = {}
        
        for entry in required_cvterms:
            cvterm_id = entry[0]
            if cvterm_id not in existing_cvterm_ids:
                cur.execute( """
                    INSERT INTO cvterm 
                    (cvterm_id,cv_id,dbxref_id,name,definition) 
                    VALUES (%s,%s,%s,%s,%s) 
                    """, entry )

            name = entry[3]
            cvterms[name] = cvterm_id
        
        self.cvterms = cvterms
        
        
    def _init_organisms(self,cur):
        """
        Make sure the organism table has all necessary entries
        
        build a dictionary (self.organisms) to efficiently lookup organism IDs
        
        this is used in constructor
        """        
        
        organisms = {}
        
        for entry in required_organisms:
            c_name, is_name = entry[-2:]
            
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
                    (abbreviation,genus,species,common_name,infraspecific_name) 
                    VALUES (%s,%s,%s,%s,%s) 
                    RETURNING organism_id
                """, entry )
                result = cur.fetchone()
                org_id = result[0]
            
            organisms[(c_name,is_name)] = org_id
            
        # check consistency with get_all_organisms()
        for org in get_all_organisms():
            if org not in organisms.keys():
                raise Exception( f"internal consistency failed. id for organism {org}" )
            
        self.organisms = organisms
        
        
        
    def insert_sequences( self, organism, metadata_df, fasta_filepath, is_protein ):
        """
        Insert genetic sequences into the database
        
        Metadata must be provided for a subset of gene IDs present in the fasta file
        
        DNA sequences should be inserted before thei corresponding protein sequences. 
        When inserting protein sequences, the "transcript" field of the fasta headers
        is used to relate the new sequences with existing DNA sequences.
        
        Arguments:
        ----------
        organism -- (pair of str) organism common-name and infraspecific-name
                            e.g. ["Maize","v3"]
        metadata_df -- (DataFrame) a dataframe with the following columns:
                            "gene_id","name","class","family"
        fasta_filepath -- (str) the path to a fasta file containing sequences
        is_protein -- (bool) true if the given fasta contains protein sequences
                             false if the it contains dna sequences
        """
        
        # make sure the required metadata has been provided
        for col in ["gene_id","name","class","family"]:
            if col not in metadata_df.columns:
                raise Exception( f'metadata is missing a required column: "{col}"' )
        df = metadata_df
        all_gene_ids = df["gene_id"].value
                
        # connect to the database
        with psycopg2.connect(self.conn_str) as conn:
            with conn.cursor() as cur:
                
                # iterate through relevant fasta records
                for rec in read_records_for_gene_ids(fasta_filepath,all_gene_ids):
                    gene_id = get_gene_id_from_record(rec)
                    transcript_id = rec.id
                    name,clazz,family = df.loc[ df["gene_id"] == gene_id, ["name","class","family"] ].values[0]
                    # insert one sequence
                    _insert_sequence( cur, org_id, str(rec.seq), gene_id, transcript_id, name, clazz, family, is_protein )
            
            
            
    def _insert_sequence( cur, org_id, sequence, gene_id, transcript_id, name, clazz, family, is_protein ):
        """
        used in insert_sequences()
        """
        
        if _get_feature_id( cur, transcript_id ) is not None:
            print( f'transcript ID "{transcript_id}" already taken, skipping insert...' )
        
        if is_protein:
            type_id = self.cvterms["polypeptide"]
        else:
            type_id = self.cvterms["DNA"]
            
        # insert feature with sequence and identifiers
        cur.execute("""
            INSERT INTO chado.feature 
            (organism_id, name, uniquename, residues, seqLen, 
                md5checksum, type_id)
            VALUES (%s,%s,%s,%s,%s,%s,%s)
            RETURNING feature_id
        """, (org_id,name,transcript_id,sequence,len(sequence),
              self._checksum(dna_seq),type_id))
        (feature_id,) = cur.fetchone()
        
        # insert featureprops
        _insert_featureprop( cur, feature_id, "class", clazz )
        _insert_featureprop( cur, feature_id, "supported_by_domain_match", family )
        _insert_featureprop( cur, feature_id, "gene_by_genome_location", gene_id )
        
    def _insert_featureprop( cur, feature_id, type_name, value ):
        """
        insert a row into the featureprop table
        
        used in _insert_sequence()
        """
        type_id = self.cvterms[type_name]
        cur.execute("""
            INSERT INTO chado.featureprop (feature_id, type_id, value)
            VALUES (%s,%s,%s);
        """, (feature_id,type_id,value))
        
        
    def _checksum( seq ):
        """
        get a value to insert into table "feature", column "md5checksum"
        
        corresponds with the given sequence in column "residues"
        
        Using in _insert_sequence()
        """
        return  hashlib.md5(seq.encode('utf-8')).hexdigest()
        
        
        
    def _get_feature_id( cur, uniquename ):
        """
        used in _insert_sequence()
        """
        
        cur.execute("""
            SELECT feature_id
            FROM chado.feature
            WHERE uniquename = %s;
        """, (uniquename,))
        result = cur.fetchone()
        
        if result is None:
            return None
        return result[0]
        
                
        
        
        
        
        
    