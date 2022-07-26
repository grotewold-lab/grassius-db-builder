# local imports
from .docker_util import *
from ..fasta import *
from ..grassius import *
from .chado_cvterms import init_dbxrefs,init_cvs,init_cvterms
from .chado_organisms import init_organisms

import psycopg2
import hashlib
import gzip
from Bio import SeqIO

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
        
    
    def insert_tfomes( self, old_grassius_tfomes ):
        """
        Insert tfome sequences into the chado tabls 'feature' and 
        'feature_relationship'
        
        also insert data into the grassius-specific 'tfome_metadata' table
        
        Tfomes will be related to This function should be called after all other 
        sequences have been inserted using insert_sequences()
        
        Arguments:
        ----------
        old_grassius_tfomes -- (DataFrame) output from 
                                gdb.grassius.get_old_grassius_tfomes
        """
        
        org_id = self.organisms["Maize_v3"]
        
                
        # connect to the database
        with psycopg2.connect(self.conn_str) as conn:
            with conn.cursor() as cur:
                
                # iterate through tfomes
                for row in old_grassius_tfomes.index:
                    utname,dna_seq,prot_seq,gene_id,tn = old_grassius_tfomes.loc[
                        row,['utname','sequence','translation','gene_id','transcript_number']]

                    # insert dna seq
                    tfome_dna_fid = self._insert_tfome_sequence( cur, org_id, utname, dna_seq, False)

                    # insert protein seq
                    tfome_prot_fid = self._insert_tfome_sequence( cur, org_id, utname, prot_seq, True)

                    # insert relationship
                    # protein -> (derives from) -> dna
                    self._insert_frel( cur, tfome_prot_fid, "derives_from", tfome_dna_fid )
                    
                    # insert relationship
                    # tfome dna -> (clone) -> genomic dna transcript
                    self._insert_tfome_frel( cur, tfome_dna_fid,utname, gene_id,tn )
            
            
    def _insert_tfome_frel( self, cur, tfome_dna_fid,utname, gene_id,tn ):
        """
        If necessary, insert one new feature_relationship
        to relate a tfome to some genomic dna
        
        tfome dna -> (clone) -> genomic dna transcript
        
        used in insert_tfomes()
        """
        if len(gene_id.strip()) == 0:
            return
        if "_FG" in gene_id:
            gt_name = gene_id.split('_')[0] + "_FG" + tn
        else:
            gt_name = f'{gene_id}_{tn}'
        gt_fid = self._get_feature_id( cur, gt_name )
        if gt_fid is None:
            #print(f'could not find dna transcript "{gt_name}" for tfome "{utname}". Skipping adding feature_relationship...')
            pass
        else:
            self._insert_frel( cur, tfome_dna_fid, "clone", gt_fid )
            
            
        
            
    def _insert_tfome_sequence( self, cur, org_id, utname, sequence, is_protein ):
        """
        If necessary, insert one new chado feature with the given tfome sequence.
        
        return the relevant feature_id
        
        used in insert_tfomes()
        """
        
        if is_protein:
            type_id = self.cvterms["polypeptide"]
            uniquename = utname + "_P"
        else:
            type_id = self.cvterms["DNA"]
            uniquename = utname
        
        existing_id = self._get_feature_id( cur, uniquename )
        if existing_id is not None:
            return existing_id
            
        # insert feature with sequence and identifiers
        cur.execute("""
            INSERT INTO feature 
            (organism_id, name, uniquename, residues, seqLen, 
                md5checksum, type_id)
            VALUES (%s,%s,%s,%s,%s,%s,%s)
            RETURNING feature_id
        """, (org_id,utname,uniquename,sequence,len(sequence),
              self._checksum(sequence),type_id))
        (feature_id,) = cur.fetchone()
        
        return feature_id
            
        
        
    def insert_sequences( self, organism, metadata_df, fasta_filepath, is_protein ):
        """
        Insert genetic sequences into the database
        
        This should not be used for tfomes' sequences. Instead use insert_tfomes().
        
        Metadata must be provided for a subset of gene IDs present in the fasta file
        
        DNA sequences should be inserted before their corresponding protein sequences. 
        When inserting protein sequences, the "transcript" field of each fasta entry's
        description is used to relate the new sequences with existing DNA sequences.
        
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
                #for rec in read_records_for_gene_ids(fasta_filepath,all_gene_ids):
                with open(fasta_filepath) as fin:
                    for rec in SeqIO.parse(fin, "fasta"):
                        gene_id = get_gene_id_from_record(rec)
                        transcript_id = rec.id
                        
                        # find corresponding metadata
                        match = df[df['gene_id']==gene_id]
                        if len(match.index) > 0:
                            name,clazz,family = df.loc[ match.index[0], ["name","class","family"] ].values
                        else:
                            name,clazz,family = transcript_id,None,None

                        # insert one sequence
                        fid = self._insert_sequence( cur, org_id, str(rec.seq), gene_id, 
                                                    transcript_id, name, clazz, family, is_protein )

                        # insert feature_relationship for protein sequences
                        # protein -> (derives from) -> dna
                        if is_protein:
                            related_tid = get_related_tid_from_record( rec )
                            if related_tid is not None:
                                related_fid = self._get_feature_id( cur, related_tid )
                                if related_fid is None:
                                    raise Exception(f'missing feature with uniquename "{related_tid}"')
                                self._insert_frel( cur, fid, "derives_from", related_fid )
                        
                        
            
    def _insert_sequence( self, cur, org_id, sequence, gene_id, transcript_id, name, clazz, family, is_protein ):
        """
        if necessary, insert one new chado feature with the given sequence
        
        also insert featureprops if necessary
        
        return the relevant feature_id
        
        used in insert_sequences()
        """
        
        existing_id = self._get_feature_id( cur, transcript_id )
        if existing_id is not None:
            return existing_id
        
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
        self._insert_fprop( cur, feature_id, "gene_by_genome_location", gene_id )
        if clazz is not None:
            self._insert_fprop( cur, feature_id, "class", clazz )
        if family is not None:
            self._insert_fprop( cur, feature_id, "supported_by_domain_match", family )
            
        return feature_id
        
            
    def _insert_frel( self, cur, subject_fid, type_name, object_fid ):
        """
        if necessary, insert a row into the feature_relationship table
        
        used in insert_sequences() and insert_tfomes()
        """
        type_id = self.cvterms[type_name]
        
        # check for exiting frel
        cur.execute("""
            SELECT feature_relationship_id 
            FROM feature_relationship
            WHERE (subject_id = %s) AND (type_id = %s) AND (object_id = %s)
        """, (subject_fid,type_id,object_fid) )
        existing = cur.fetchone()
        
        # insert new frel if necessary
        if existing is None:
            cur.execute("""
                INSERT INTO feature_relationship (subject_id, type_id, object_id)
                VALUES (%s,%s,%s)
            """, (subject_fid,type_id,object_fid))
        
                                      
    def _insert_fprop( self, cur, feature_id, type_name, value ):
        """
        insert a row into the featureprop table
        
        used in _insert_sequence()
        """
        type_id = self.cvterms[type_name]
        
        cur.execute("""
            INSERT INTO featureprop (feature_id, type_id, value)
            VALUES (%s,%s,%s)
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
        return the feature ID of an existing feature with the 
        given uniquename, or None if it does not exist
        
        used in _insert_sequence(), _insert_tfome_sequence(), 
        insert_domain_annotations()
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
                              family_desc_df, old_grassius_names, 
                              old_grassius_tfomes, gene_interactions, 
                              domain_descriptions ):
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
        old_grassius_tfomes -- (DataFrame) output from 
                                gdb.grassius.get_old_grassius_tfomes()
        gene_interactions -- (DataFrame) loaded from private
                                input "gene_interactions"
        domain_descriptions -- (DataFrame) output from
                                gdb.grassius.get_domain_descriptions()
        """
        
        all_family_names = set(metadata_df["family"])
        
        # extract numbers from protein names, if necessary
        if 'suffix' not in metadata_df.columns:
            metadata_df = parse_all_protein_names(metadata_df)
            
        # for performance, extract a dictionary from metadata
        protein_name_dict = get_protein_name_dict(metadata_df)
            
        # connect to the database
        with psycopg2.connect(self.conn_str) as conn:
            with conn.cursor() as cur:
                
                # build tables
                build_domain_descriptions( cur, domain_descriptions )
                build_tfome_metadata( cur, old_grassius_tfomes )
                build_gene_clone( cur, old_grassius_tfomes )
                build_searchable_clones( cur, protein_name_dict, old_grassius_tfomes )
                build_gene_interaction( cur, protein_name_dict, gene_interactions )
                build_seq_features( cur )
                build_uniprot_ids( cur, protein_name_dict )
                build_comment_system_urls( cur, all_family_names )
                build_default_maize_names( cur, metadata_df, gene_versions, all_family_names, old_grassius_tfomes )
                build_gene_name( cur, metadata_df, old_grassius_names )
                build_family_tables( cur, metadata_df, family_desc_df )
        
                
    def insert_domain_annotations( self, all_domain_annos ):
        """
        Insert domain annotations into the "featureprop" table.
        
        This should be called after inserting protein sequences using
        insert_sequences()
        
        Any existing domain annotations will be deleted and replaced.
        
        Arguments:
        ----------
        all_domain_annos: (dict) json data output from
                                 gdb.grassius.get_domain_annotations()
        """
        
            
        # connect to the database
        with psycopg2.connect(self.conn_str) as conn:
            with conn.cursor() as cur:
        
                # remove existing annotations from database
                cur.execute("""
                    DELETE FROM featureprop
                    WHERE type_id = 61467;
                """)

                # iterate through the given annotations
                for tid,anno in all_domain_annos.items():
                    all_uns = [tid]
                    if '_P' in tid:
                        all_uns.append( tid.replace( '_P', '_T' ) )
                    i = 0
                    
                    for un in all_uns:

                        # find related feature
                        fid = self._get_feature_id( cur, un )
                        if fid is None:
                            continue

                        # iterate through annotations for this transcript
                        if not isinstance(anno,list):
                            anno = [anno]
                        for entry in anno:
                            name = entry['@name']
                            acc = entry['@acc']#.split(".")[0]
                            all_doms = entry['domains']
                            if not isinstance(all_doms,list):
                                all_doms = [all_doms]
                            for d in all_doms:                        

                                # insert one featureprop entry
                                dom = {
                                    "name": name, 
                                    "accession": acc, 
                                    "start": d['@alisqfrom'], 
                                    "end": d['@alisqto']
                                }
                                sdom = str(dom).replace("'",'"')
                                cur.execute("""
                                    INSERT INTO featureprop ( feature_id, type_id, value, rank ) 
                                    VALUES ( %s, 61467, %s, %s );
                                """, (fid,sdom,i))

                                i += 1
            
    def insert_secondary_structures(self):
        """
        Build into the grassius-specific table "seq_features"
        
        This should be called after inserting protein sequences using
        insert_sequences()
        
        The table will be populated with Jan2022 secondary structure
        data in html form.
        
        This is a band-aid to make the website work for now. Eventually 
        the website should be updated and secondary structures
        should be stored in json form similar to domains, within the 
        chado table "featureprop".
        """
        
        with psycopg2.connect(self.conn_str) as conn:
            with conn.cursor() as cur:
                
                df = get_secondary_structures()
                
                # iterate through secondary structure data
                for row in df.index:
                    tid,struct = df.loc[row,['transcript_id','structure']]
                    
                    # find corresponding feature id
                    fid = self._get_feature_id( cur, tid )
                    if fid is None:
                        continue
                        
                    # insert a row
                    cur.execute("""
                        INSERT INTO seq_features 
                        (feature_id,secondary_structures)
                        VALUES (%s,%s)
                    """, (fid,struct) )
                    
                    
        
        
        
        
        
        
        
        
        
    