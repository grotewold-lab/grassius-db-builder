# local imports
from .docker_util import *

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
            
        # connect to the database
        conn_str = f"dbname=postgres host=localhost port={port} user=postgres password=postgres"
        with psycopg2.connect(conn_str) as conn:
            with conn.cursor() as cur:
                try:
                    cur.execute("SELECT * from cvterm limit 10")
                except:
                    raise Exception("could not connect to the database")
                
        # set instance variables
        self.port = port
        self.docker_container = container
        self.conn_str = conn_str
        
    