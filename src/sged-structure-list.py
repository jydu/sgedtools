#! /usr/bin/python

""" Created on 09/03/20 by jdutheil

    List all amino-acid residues in a structure
"""

import getopt, sys, glob

from Bio.PDB import *

parser = PDBParser()

cmd_args = sys.argv
arg_list = cmd_args[1:]

unix_opt = "p:o:a:c"
full_opt = ["pdb=", "output=", "chain=", "csv"]
try:
    arguments, values = getopt.getopt(arg_list, unix_opt, full_opt)
except getopt.error as err:
    print(str(err))
    sys.exit(2)

tabsep = True  # TSV by default
group_col = "PDB"
for arg, val in arguments:
    if arg in ("-p", "--pdb"):
        pdb_file = val
        print("PDB file: %s" % val)
    elif arg in ("-o", "--output"):
        output_file = val
        print("Output SGED file: %s" % output_file)
    elif arg in ("-a", "--chain"):
        chain_sel = val
        print("PDB chain to use: %s" % chain_sel)
    elif arg in ("-c", "--csv"):
        tabsep = False

if tabsep:
    print("SGED file is in TSV format")
    delim = "\t"
else:
    print("SGED file is in CSV format")
    delim = ","


print("Parsing PDB file %s..." % pdb_file)
structure = parser.get_structure("STRUCT", pdb_file)

# First we need to check that there is only one model:
if len(structure) > 1:
    print(
        "Warning, %s models in PDB file %s. Using the first one."
        % (len(structure), pdb_file)
    )

model = structure[0]
chain = model[chain_sel]

# Then we write down all residues (TODO: eventually with some basic info):
def res_to_str(id):
    s = str(id[1])
    if id[2] != " ":
        s = s + id[2]
    return s


with open(output_file, "w") as handle:
    handle.write("# SGED input PDB = %s\n" % pdb_file)
    handle.write("# SGED input PDB chain = %s\n" % chain_sel)
    handle.write("%s\n" % group_col)
    for residue in chain:
        if is_aa(residue):
            handle.write(
                "[%s%s]\n" % (residue.get_resname(), res_to_str(residue.get_id()))
            )

print("Done.")
