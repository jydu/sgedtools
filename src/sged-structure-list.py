#!/usr/bin/env python3

""" Created on 09/03/20 by jdutheil

    List all amino-acid residues in a structure.
    This program is part of the SgedTools package.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

import getopt, sys, glob

from Bio.PDB import *

cmd_args = sys.argv
arg_list = cmd_args[1:]

unix_opt = "p:i:f:o:ch"
full_opt = ["pdb=", "pdb-id=", "pdb-format=", "output=", "csv", "help"]

def usage() :
    print(
"""
sged-structure-list

    List all amino-acid residues in a structure

Available arguments:
    --pdb (-p): Input protein data bank file (required).
    --pdb-format (-f): Format of the protein data bank file (default: PDB).
        Either PDB or mmCif is supported. In addition, remote:PDB or remote:mmCif
        allow to directly download the structure file from the Protein Data Bank.
        In this case, --pdb-id indicates the PDB id.
    --pdb-id (-i): Specify the id of the PDB file to retrieve remotely.
    --output (-o): Output SGED file (required).
    --csv (-c): Input SGED file is with comas instead of tabs (default: tabs).
    --help (-h): Print this message.
"""
    )
    sys.exit()


try:
    arguments, values = getopt.getopt(arg_list, unix_opt, full_opt)
except getopt.error as err:
    print(str(err))
    sys.exit(2)

pdb_format = "PDB"
tabsep = True  # TSV by default
for arg, val in arguments:
    if arg in ("-p", "--pdb"):
        pdb_file = val
        print("PDB file: %s" % pdb_file)
    elif arg in ("-i", "--pdb-id"):
        pdb_id = val
        print("PDB id: %s" % val)
    elif arg in ("-f", "--pdb-format"):
        pdb_format = val
        if val != "PDB" and val != "mmCif" and val[0:7] != "remote:":
            print(
                "Structure format should be either PDB or mmCif, or remote:PDB, remote:mmCIF, etc. if you would like to retrieve the file from RCSB"
            )
            exit(-1)
        print("PDB format: %s" % pdb_format)
    elif arg in ("-o", "--output"):
        output_file = val
        print("Output SGED file: %s" % output_file)
    elif arg in ("-c", "--csv"):
        tabsep = False
    elif arg in ("-h", "--help"):
        usage()

if tabsep:
    print("SGED file is in TSV format")
    delim = "\t"
else:
    print("SGED file is in CSV format")
    delim = ","

# Check options:

if not 'pdb_file' and not 'pdb_id' in globals():
    print("Error: a structure file or id should be provided.")
    usage()

if not 'output_file' in globals():
    print("Error: an output file should be provided.")
    usage()

# Parse the PDB and compute

if pdb_format.startswith("remote:"):
    remote_format = pdb_format[7:]
    pdb_server = PDBList(
        server="ftp://ftp.wwpdb.org", pdb = None, obsolete_pdb = False, verbose = True
    )
    pdb_format = remote_format

if pdb_format.upper() == "PDB":
    pdb_format = "PDB" #Case needs to be respected for remote access
    parser = PDBParser()
elif pdb_format.upper() == "MMCIF":
    pdb_format == "mmCif" #Case needs to be respected for remote access
    parser = MMCIFParser()
else:
    print("ERROR!!! Unsupported structure format: %s" % pdb_format)
    exit(-1)

if "pdb_server" in locals():
    pdb_file = pdb_server.retrieve_pdb_file(
        pdb_code = pdb_id,
        obsolete = False,
        pdir = ".",
        file_format = remote_format,
        overwrite = False,
    )
    print("Downloaded PDB file %s..." % pdb_file)

structure = parser.get_structure("STRUCT", pdb_file)

# First we need to check that there is only one model:
if len(structure) > 1:
    print(
        "Warning, %s models in PDB file %s. Using the first one."
        % (len(structure), pdb_file)
    )

model = structure[0]

# Then we write down all residues:

def res_to_str(id):
    s = str(id[1])
    if id[2] != " ":
        s = s + id[2]
    return s


with open(output_file, "w") as handle:
    handle.write("# SGED input PDB = %s\n" % pdb_file)
    handle.write("Group%sChain\n" % delim)
    for chain in model:
        for residue in chain:
            if is_aa(residue):
                handle.write(
                        "[%s:%s%s]%s%s\n" % (chain.id, residue.get_resname(), res_to_str(residue.get_id()), delim, chain.id)
                )

print("Done.")
