# this file contains utility functions related to chado

import docker

container_name = "gdb_chado_container"
    
    
def get_existing_container():
    """
    Get the relevant docker container if it exists
    otherwise return None
    """
    for c in docker.from_env().containers.list(all=True):
        if c.name == container_name:
            return c
    return None

    
def init_template_db( replace=False ):
    """
    Start running an empty chado database locally in docker
    
    Python must be running in administrator mode
    
    Arguments:
    ----------
    replace -- (bool) (optional) True if an existing container from previous runs should be replaced
                                otherwise an exception will be thrown if the container already exists
    """
    
    
    client = docker.from_env()
    
    # check if container is already running
    ec = get_existing_container()
    if ec is not None:
        if replace:
            if ec.status != "exited":
                ec.kill()
            client.containers.prune()
        else:
            raise Exception("container already exists")
            
        
    # start new container
    return client.containers.run("postgres:9", detach=True, name=container_name)