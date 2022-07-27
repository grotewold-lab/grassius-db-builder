import os
import importlib.util
import sys

def apply_all_patches(cb):
    """
    Apply all patches in order
    
    Arguments:
    ----------
    cb - a ChadoBuilder instance, connected to the 
            database that should be patched.
    """
    folder = os.path.dirname(os.path.abspath(__file__))
    all_patch_files = _get_patch_filenames(folder)
    for fname in all_patch_files:
        path = os.path.join(folder,fname)
        print(f"applying patch:\n\t{path}\n")
        
        spec = importlib.util.spec_from_file_location("patch", path)
        patch = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(patch)
        patch.apply_patch(cb)
        
            
def _get_patch_filenames(folder):
    result = []
    for fname in os.listdir(folder):
        if fname.startswith('p') and fname.endswith('.py'):
            result.append(fname)
    return sorted(result)
    