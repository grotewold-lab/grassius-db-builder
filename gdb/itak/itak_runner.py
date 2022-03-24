# this file contains logic to run iTAK

# a pre-built iTAK executable comes from public inputs
#  rules will be replaced (iTAK/database/TF_Rule.txt)
#  hmm files will be replaced (iTAK/database/*.hmm)

# this file contains logic to replace those files


# local imports
from ..input_manager import InputManager

import zipfile
import os
import subprocess
from pathlib import Path


class ItakRunner:
    
    def __init__(self):
        
        # locate or create local iTAK folder
        im = InputManager()
        zip_path = im["itak_git_repo"]
        self.folder = str(Path(zip_path).parent.absolute()) + "/iTAK"
        if os.path.exists( self.folder ):
            print( f"using existing itak folder:\n\t{self.folder}" )
        else:
            print( f"creating local itak folder:\n\t{self.folder}" )
            with zipfile.ZipFile(zip_path, 'r') as zip_ref:
                zip_ref.extractall(self.folder)
        
        # verify that the iTAK executable exists
        exe_path = self.folder + "/iTAK-master/iTAK.pl"
        if not os.path.exists( exe_path ):
            raise Exception( f"missing iTAK executable:\n\t{exe_path}" )
            
            
    def run_test(self):
        """
        Run iTAK with included sample data, 
        as described on the github README
        """
        
        wd = self.folder + "/iTAK-master"
        
        command = [ "perl", "iTAK.pl", "test_seq" ]
        p = subprocess.Popen(command, cwd=wd)
        p.wait()
        
        