# this file contains utility functions related to chado

# local imports
from ..input_manager import InputManager

import docker
import tempfile
import shutil

image_name = "gdb-chado-image"
container_name = "gdb-chado-container"
default_port = 8642
    
    
def get_existing_template_db_container():
    """
    Get the relevant docker container if it exists
    otherwise return None
    """
    try:
        return docker.from_env().containers.get(container_name)
    except docker.errors.NotFound:
        return None

    
    
def init_template_db_container( replace=False, port=None ):
    """
    Start running an empty chado database locally in docker
    
    Python must be running in administrator mode
    
    Arguments:
    ----------
    replace -- (bool) (optional) True if an existing container from previous runs should be replaced
                                otherwise an exception will be thrown if the container already exists
    port -- (int) (optional) the host port that will be usable for access to the database
    """
    
    if port is None:
        port = default_port
    
    client = docker.from_env()
    
    # check if container is already running
    ec = get_existing_template_db_container()
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
    return client.containers.run(image_name,
                                 name=container_name, 
                                 ports={'5432/tcp':8642},
                                 detach=True, 
                                 tty=True)


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
    
    # build image
    docker.from_env().images.build( path=folder, tag=image_name )
    
    #debug
    print( folder )
    #shutil.rmtree(folder)
