# this file contains utility functions related to chado

# local imports
from ..input_manager import InputManager

import docker
import tempfile
import shutil
import time
import numpy as np

image_name = "gdb-chado-image"
container_name = "gdb-chado-container"
    
    
def container_maps_to_port( container, port ):
    """
    Return true if the given docker container is mapped to the given port on the host.
    
    This should only be used for a container that was created by init_db_container()
    
    Arguments:
    ----------
    container -- (docker.models.containers.Container) the docker container in question
    replace -- (int) the expected port number on the host machine
    """
    
    sport = str(port)
    
    for entry in container.ports['5432/tcp']:
        if entry['HostPort'] != sport:
            return False
        
    return True
        
    
def get_existing_db_container():
    """
    Get the relevant docker container if it exists
    otherwise return None
    """
    try:
        return docker.from_env().containers.get(container_name)
    except docker.errors.NotFound:
        return None

    
    
def init_db_container( port, replace=False ):
    """
    Start running an empty chado database locally in docker
    
    Python must be running in administrator mode
    
    Arguments:
    ----------
    port -- (int) the host port that will be used to access the database
    replace -- (bool) (optional) True if an existing container from previous runs should be replaced
                                otherwise an exception will be thrown if the container already exists
    """
    
    client = docker.from_env()
    
    # check if container is already running
    ec = get_existing_db_container()
    if ec is not None:
        if replace:
            if ec.status != "exited":
                ec.kill()
            client.containers.prune()
        else:
            raise Exception("container already exists")
            
    # build image if necessary
    if not image_exists():
        init_template_db_image()
        
    # start new container and keep it running
    container = client.containers.run(image_name,
                                 name=container_name, 
                                 ports={'5432/tcp':port},
                                 environment={"POSTGRES_USER":"postgres",
                                              "POSTGRES_PASSWORD":"postgres"},
                                 detach=True, 
                                 tty=True)
    print( "created new docker container for database." )
    
    # wait for initialization process to finish
    _moniter_initialization( container )
    
    return container


def _moniter_initialization( container ):
    """
    Wait for the given container to finish initializing
    Report progress while waiting
    """
    linecount = 0
    expected_lines = 5860
    
    # pick line counts to report progress
    report_points = { 
        int(pct*expected_lines) : str(int(pct*100))+"% complete" 
        for pct in np.arange( .05,1,.05 )
    }
    
    print( "initializing database..." )
    
    # start recieving log entries from the container
    buffer = bytearray()
    for b in container.logs(stream=True):
        if not isinstance(b,bytes):
            print( "got unexpected response format from logs. Skipping monitering of database initialization. You may need to manually wait for initialization and try again" )
            return 
            
        buffer.extend(b)
        if b.endswith(b'\n'): # hit line break
            line = buffer.decode().strip()
            buffer = bytearray()
            linecount += 1
            
            # report progess
            if linecount in report_points.keys():
                print( report_points[linecount] )   
            
            
            # check if finished
            if (linecount > 100) and (line == "LOG:  database system is ready to accept connections"):
                break
    print( "database is ready!" )

    
    
def image_exists():
    """
    Return true if the relevant docker image already exists
    
    This should only be used in init_template_db_container
    """
    try:
        docker.from_env().images.get(image_name)
        return True
    except docker.errors.ImageNotFound:
        return False

    
def init_template_db_image():
    """
    Make sure there is a docker image that will be used in init_template_db_container.
    
    This should only be used in init_template_db_container
    """
    
    print( f'building docker image "{image_name}"...' )
    
    # make a temporary folder to contain docker build context
    folder = tempfile.mkdtemp()
    dockerfile = folder + "/Dockerfile"
    sql_file = folder + "/init.sql"
    
    # paste sql file
    shutil.copy( InputManager()["chado_template"], sql_file )
    
    # build dockerfile
    with open(dockerfile,"w") as fout:
        fout.write("FROM postgres:9\n")
        fout.write("COPY init.sql /docker-entrypoint-initdb.d/\n")
        
        # workaround error when running image
        fout.write("RUN sed -i 's/CREATE INDEX bingroup_boxrange/--/g' /docker-entrypoint-initdb.d/init.sql\n")
        fout.write("RUN sed -i 's/CREATE INDEX binloc_boxrange/--/g' /docker-entrypoint-initdb.d/init.sql\n")
    
    # build image
    docker.from_env().images.build( path=folder, tag=image_name )
    
    # delete temporary folder
    shutil.rmtree(folder)
