# locate, download, and verify input files
import os
import pandas as pd
import gzip
import hashlib
import shutil
import urllib.request as request
from contextlib import closing

class InputManager:

    
    def __init__(self):
        
        # locate local folder and load all metadata
        self.input_dir = os.path.join( os.path.dirname(__file__), "../inputs" )
        path = os.path.join( self.input_dir, "input_list.csv" )
        self.df = pd.read_csv(path).fillna("")
        
        # prepare headers to emulate a browser when downloading files
        self.dl_headers={
            'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_9_3) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/35.0.1916.47 Safari/537.36'
        }
        
    
    
    def prepare_all_inputs(self):
        """
        make sure all inputs are ready for use
        download and extract public inputs if necessary
        verify integrity of all inputs
        """
        for name in self.list_input_names():
            self.get_input_filepath(name)
            
        backup_folder = os.path.abspath(self.input_dir)
        print( f'all inputs verified. Recommend backing up this folder:\n\t{backup_folder}' )
            
    
    def list_input_names(self):
        """
        get a list of short names that correspond with input files
        """
        return list(self.df["name"])
    
    
    def get_input_filepath(self, name):
        """
        download and extract the file if necessary
        verify integrity
        return the local path for the given input file
        """
        
        # get metadata
        if name not in self.df["name"].values:
            raise Exception(f'unrecognized input name "{name}"')
        filename,url,md5sum = self.df.loc[
            self.df["name"] == name,
            ["filename","public_url","md5sum"]].values[0]
        
        # pick the expected local file path
        is_public = (len(url) > 0)
        subdir_name = "public_inputs" if is_public else "private_inputs"
        filepath = os.path.join( self.input_dir, subdir_name, filename )
        
        # if the local file doesn't exist...
        if not os.path.exists( filepath ):
            if not is_public:
                raise Exception(f'missing private input "{name}".\n\texpected local path: {filepath}')
            
            # download the public file
            compressed_path = filepath + ".gz"
            print( f'downloading public input "{name}"...\n\turl: {url}' )
            req = request.Request( url, data=None, headers=self.dl_headers )
            with closing(request.urlopen(req)) as rin:
                with open(compressed_path, 'wb') as fout:
                    shutil.copyfileobj(rin, fout)
            
            # decompress the downloaded file
            print( f'extracting downloaded archive for "{name}"...' )
            with gzip.open(compressed_path, 'rb') as fin:
                with open(filepath, 'wb') as fout:
                    shutil.copyfileobj(fin,fout)
            
        # check integrity
        print( f'checking integrity of input "{name}"...' )
        with open(filepath, "rb") as f:
            file_hash = hashlib.md5()
            while chunk := f.read(8192):
                file_hash.update(chunk)
            if file_hash.hexdigest() != md5sum:
                raise Exception(f'integrity check failed for input "{name}"\n\tlocal path: {filepath}')
       
        return filepath
        
        
        
        
        
        