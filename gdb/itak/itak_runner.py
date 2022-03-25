# this file contains logic to run iTAK

# a pre-built iTAK executable comes from public inputs
#  rules will be replaced (iTAK/database/TF_Rule.txt)
#  hmm files will be replaced (iTAK/database/*.hmm)

# this file contains logic to replace those files


# local imports
from ..input_manager import InputManager
from .itak_util import read_itak_output, build_rules_file
from ..hmmer import concatenate_hmms

import zipfile
import os,sys,stat
import subprocess
import shutil
import tempfile
from pathlib import Path



class ItakRunner:
    
    def __init__(self, reset=False):
        """
        create a new instance of ItakRunner
        
        Arguments:
        ----------
        reset -- (optional) (bool) if true, the iTAK database will be reset. 
                      after resetting, the database will contain sample data 
                      from the github repo
        """
        
        # locate or create local iTAK folder
        im = InputManager()
        zip_path = im["itak_git_repo"]
        self.folder = str(Path(zip_path).parent.absolute()) + "/iTAK"
        if os.path.exists( self.folder ):
            if reset:
                print( f"removing existing itak folder:\n\t{self.folder}" )
                shutil.rmtree( self.folder )
            else:
                print( f"using existing itak folder:\n\t{self.folder}" )
        if not os.path.exists(self.folder):
            print( f"creating local itak folder:\n\t{self.folder}" )
            with zipfile.ZipFile(zip_path, 'r') as zip_ref:
                zip_ref.extractall(self.folder)
        
        # verify that the iTAK executable exists
        exe_path = self.folder + "/iTAK-master/iTAK.pl"
        if not os.path.exists( exe_path ):
            raise Exception( f"missing iTAK executable:\n\t{exe_path}" )
            
        # make sure binaries have execute permission
        bin_folder = self.folder + "/iTAK-master/bin/"
        os.chmod(bin_folder + "/hmmbuild", stat.S_IXUSR)
        os.chmod(bin_folder + "/hmmpress", stat.S_IXUSR)
        os.chmod(bin_folder + "/hmmscan", stat.S_IXUSR)
        
            
    def run_test(self):
        """
        Run iTAK with included sample data, 
        as described on the github README
        """
        
        wd = self.folder + "/iTAK-master"
        
        command = [ "perl", "iTAK.pl", "test_seq" ]
        p = subprocess.Popen(command, cwd=wd)
        p.wait()
        
        return read_itak_output( wd + "/test_seq_output" )
    
        
    def set_database(self, hmm_filepath, rules_df):
        """
        Replace the iTAK database, so future runs will apply the 
        given hm models and criteria
        
        Arguments:
        ----------
        hmm_filepath -- (str) the path to an hmm file
        rules_df -- (DataFrame) the criteria for classifying transcripts
                           output from gdb.hmmer.read_family_criteria()
                           must contain columns "Required", "Forbidden", and "GRASSIUS"
        """
        
    
        database_folder = self.folder + "/iTAK-master/database"
        
        # append hmm file to the original hmm file in the database
        temp_folder = tempfile.mkdtemp()
        original_hmm_path = database_folder + "/Tfam_domain.hmm"
        combined_hmm_path = temp_folder + "/combined.hmm"
        concatenate_hmms( [original_hmm_path,hmm_filepath], combined_hmm_path )
        shutil.copyfile(combined_hmm_path, original_hmm_path)
        shutil.rmtree(temp_folder)
        
        # rebuild rules
        build_rules_file(database_folder, rules_df)
        
        
        
    def run_itak(self, fasta_filepath):
        """
        classify transcripts from the given fasta file
        
        database should be prepared ahead of time using set_database()
        """
        
        # copy fasta file to working dir
        wd = self.folder + "/iTAK-master"
        fname = os.path.basename( fasta_filepath )
        shutil.copyfile( fasta_filepath, os.path.join( wd, fname ) )
        
        # run iTAK
        command = [ "perl", "iTAK.pl", fname ]
        p = subprocess.Popen(command, cwd=wd)
        p.wait()
        
        return read_itak_output( wd + f"/{fname}_output" )
        
        
        
        
        
        
        
        
        
        