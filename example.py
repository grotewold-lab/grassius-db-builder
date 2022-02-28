import gdb
from gdb.blast import read_blast_output


# download and/or check integrity of all inputs
im = gdb.InputManager()
im.prepare_all_inputs()


# parse blast output example
filepath = im.get_input_filepath("blast_output_example")
br = read_blast_output(filepath)
print( br.data )


