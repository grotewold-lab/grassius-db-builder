"""
gdb.patches module

Patches are one-off procedures that should be applied 
at the end of the database builder pipeline. They may 
also be applied to an existing database. Re-applying 
patches may be necessary if the database is modified, 
otherwise they will have no effect.
"""

from .apply_patches import apply_all_patches