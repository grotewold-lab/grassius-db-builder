# local imports
from .docker_util import *
from ..fasta import *
from ..grassius import *
from .chado_cvterms import init_dbxrefs,init_cvs,init_cvterms
from .chado_organisms import init_organisms

import psycopg2
import hashlib
import gzip

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
                init_dbxrefs(cur)
                init_cvs(cur)
                self.cvterms = init_cvterms(cur)
                self.organisms = init_organisms(cur) 
                
                
    def _test_db_connection(self,cur):
        """
        raise an exception if we can't make a query
        
        this is used in constructor
        """
        try:
            cur.execute("SELECT * from cvterm limit 1")
        except:
            raise Exception("could not connect to the database")
        
        
        
        
    def query( self, sql, vars=None ):
        """
        Run a query against the database that is being built
        
        Arguments:
        ----------
        sql - (str) the query to execute, which may container placeholders
        vars - (optional) (tuple) the values to insert into the query
        """
        with psycopg2.connect(self.conn_str) as conn:
            with conn.cursor() as cur:
                cur.execute(sql, vars)
                return cur.fetchall()
            
        
    def write_snapshot( self, output_path ):
        """
        Save a snapshot of the database to a .sql.tar.gz file
        """
        
        dc = self.docker_container
        
        print( "running pg_dump in docker container..." )
        dc.exec_run( "echo 'localhost:5432:postgres:postgres:postgres' > ~/.pgpass" )
        dc.exec_run( "pg_dump -h localhost -p 5432 -U postgres --clean --if-exists -f /build_db.sql" )

        print( "extracting archive from docker container..." )
        bits, stat = dc.get_archive('/build_db.sql')
        with gzip.open(output_path, 'wb', compresslevel=6) as fout:
            for chunk in bits:
                fout.write(chunk)
                
        print( "database snapshot saved to " + output_path )
        
    
        
        
    def insert_sequences( self, organism, metadata_df, fasta_filepath, is_protein ):
        """
        Insert genetic sequences into the database
        
        Metadata must be provided for a subset of gene IDs present in the fasta file
        
        DNA sequences should be inserted before thei corresponding protein sequences. 
        When inserting protein sequences, the "transcript" field of the fasta headers
        is used to relate the new sequences with existing DNA sequences.
        
        Arguments:
        ----------
        organism -- (str) organism common-name and infraspecific-name
                            e.g. "Maize_v3"
        metadata_df -- (DataFrame) a dataframe with columns:
                            "gene_id","name","class","family"
        fasta_filepath -- (str) the path to a fasta file containing sequences
        is_protein -- (bool) true if the given fasta contains protein sequences
                             false if the it contains dna sequences
        """
        
        # validate organism
        if organism not in self.organisms.keys():
            raise Exception( f"unrecognized organism \"{organism}\". expected one of: {self.organisms.keys()}")
        org_id = self.organisms[organism]
        
        # make sure the required metadata has been provided
        for col in ["gene_id","name","class","family"]:
            if col not in metadata_df.columns:
                raise Exception( f'metadata is missing a required column: "{col}"' )
        df = metadata_df
        all_gene_ids = df["gene_id"].values
                
        # connect to the database
        with psycopg2.connect(self.conn_str) as conn:
            with conn.cursor() as cur:
                
                # iterate through relevant fasta records
                for rec in read_records_for_gene_ids(fasta_filepath,all_gene_ids):
                    gene_id = get_gene_id_from_record(rec)
                    transcript_id = rec.id
                    name,clazz,family = df.loc[ df["gene_id"] == gene_id, ["name","class","family"] ].values[0]
                    # insert one sequence
                    self._insert_sequence( cur, org_id, str(rec.seq), gene_id, transcript_id, name, clazz, family, is_protein )
            
            
            
    def _insert_sequence( self, cur, org_id, sequence, gene_id, transcript_id, name, clazz, family, is_protein ):
        """
        used in insert_sequences()
        """
        
        if self._get_feature_id( cur, transcript_id ) is not None:
            #print( f'transcript ID "{transcript_id}" already taken, skipping insert...' )
            return
        
        if is_protein:
            type_id = self.cvterms["polypeptide"]
        else:
            type_id = self.cvterms["DNA"]
            
        # insert feature with sequence and identifiers
        cur.execute("""
            INSERT INTO feature 
            (organism_id, name, uniquename, residues, seqLen, 
                md5checksum, type_id)
            VALUES (%s,%s,%s,%s,%s,%s,%s)
            RETURNING feature_id
        """, (org_id,name,transcript_id,sequence,len(sequence),
              self._checksum(sequence),type_id))
        (feature_id,) = cur.fetchone()
        
        # insert featureprops
        self._insert_featureprop( cur, feature_id, "class", clazz )
        self._insert_featureprop( cur, feature_id, "supported_by_domain_match", family )
        self._insert_featureprop( cur, feature_id, "gene_by_genome_location", gene_id )
        
    def _insert_featureprop( self, cur, feature_id, type_name, value ):
        """
        insert a row into the featureprop table
        
        used in _insert_sequence()
        """
        type_id = self.cvterms[type_name]
        cur.execute("""
            INSERT INTO featureprop (feature_id, type_id, value)
            VALUES (%s,%s,%s);
        """, (feature_id,type_id,value))
        
        
    def _checksum( self, seq ):
        """
        get a value to insert into table "feature", column "md5checksum"
        
        corresponds with the given sequence in column "residues"
        
        Using in _insert_sequence()
        """
        return  hashlib.md5(seq.encode('utf-8')).hexdigest()
        
        
        
    def _get_feature_id( self, cur, uniquename ):
        """
        used in _insert_sequence()
        """
        
        cur.execute("""
            SELECT feature_id
            FROM feature
            WHERE uniquename = %s;
        """, (uniquename,))
        result = cur.fetchone()
        
        if result is None:
            return None
        return result[0]
    
    
    
    def build_grassius_tables( self, metadata_df, gene_versions, 
                              family_desc_df, old_grassius_names ):
        """
        Build tables which are necessary for grassius, but 
        are not a part of the chado schema.
        
        See gdb/grassius/grassius_tables.py for more details
        
        If any of the tables already exist, they will be replaced.
        
        Arguments:
        ----------
        metadata_df -- (DataFrame) a dataframe with columns:
                            "gene_id","name","class","family"
        gene_versions -- (dictionary) where keys are gene_ids,
                            values are genome versions e.g. 'v3'
        family_desc_df -- (DataFrame) a dataframe loaded from private 
                            input "family_descriptions"
        old_grassius_names -- (DataFrame) names from old grassius website
                              output from gdb.grassius.get_old_grassius_names()
        """
        
        all_family_names = set(metadata_df["family"])
        
        # extract numbers from protein names, if necessary
        if 'suffix' not in metadata_df.columns:
            metadata_df = parse_protein_names(metadata_df)
            
        # connect to the database
        with psycopg2.connect(self.conn_str) as conn:
            with conn.cursor() as cur:
                build_gene_interaction( cur )
                build_seq_features( cur )
                build_uniprot_ids( cur )
                build_searchable_clones( cur )
                build_comment_system_urls( cur, all_family_names )
                build_default_maize_names( cur, metadata_df, gene_versions, all_family_names )
                build_gene_name( cur, metadata_df, old_grassius_names )
                build_family_tables( cur, metadata_df, family_desc_df )
        
        
        
        
    