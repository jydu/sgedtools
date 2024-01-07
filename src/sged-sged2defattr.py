#!/usr/bin/env python3

import pandas
import getopt, sys

""" Created on 05/01/24 by jdutheil

    Converts data from a SGED file to attributes readable by Chimera.
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

cmd_args = sys.argv
arg_list = cmd_args[1:]

unix_opt = "s:o:g:v:ch"
full_opt = ["sged=", "output=", "group=", "variable=", "csv", "help"]

def usage() :
    print(
"""
sged-sged2defattr

    Converts data from a SGED file to attributes readable by Chimera.

Available arguments:
    --sged (-s): Input SGED file (required).
    --group (-g): Column where group coordinates are stored (default: Group).
    --output (-o): Output attribute file (required).
    --variable (-v): Column to export as attributes
    --csv (-c): Input SGED file is with comas instead of tabs (default)
    --help (-h): Print this message.
"""
    )
    sys.exit()

try:
    arguments, values = getopt.getopt(arg_list, unix_opt, full_opt)
except getopt.error as err:
    print(str(err))
    sys.exit(2)

tabseq = True
group_col = "Group"
min_post_prob = 0.
for arg, val in arguments:
    if arg in ("-s", "--sged"):
        sged_file = val
        print("Input SGED file: %s" % sged_file)
    elif arg in ("-g", "--group"):
        group_col = val
        print("Group coordinates are in column: %s" % group_col)
    elif arg in ("-o", "--output"):
        output_file = val
        print("Output attribute file: %s" % output_file)
    elif arg in ("-v", "--variable"):
        var_col = val
        print("Variable to export: %s" % var_col)
    elif arg in ("-c", "--csv"):
        tabseq = False
    elif arg in ("-h", "--help"):
        usage()
        
if tabseq:
    print("SGED file is in TSV format.")
    delim = "\t"
else:
    print("SGED file is in CSV format.")
    delim = ","

# Check required arguments

if not 'sged_file' in globals():
    print("Error: a SGED input file should be specified.")
    usage()
if not 'output_file' in globals():
    print("Error: an ouput file should be specified.")
    usage()
if not 'var_col' in globals():
    print("Error: a variable columnn should be specified.")
    usage()

# Convert residue selection:
def sged2chimera(selection):
    tmp = selection[1 : (len(selection) - 1)]
    tmp = tmp.replace(" ", "")
    tmp = tmp.split(";")
    sorted_res = dict()
    for desc in tmp:
        if desc == 'NA' :
            pass
            #  ignore missing data
        else :
            m = desc.split(":")
            if len(m) == 2 :
                chain = m[0]
                res   = m[1]
                if not chain in sorted_res :
                    sorted_res[chain] = []
                if len(res) <= 3 :
                    raise NameError("Unvalid PDB residue %s, must be of the form ALA123." % res)
                else :
                    sorted_res[chain].append(res[3:])
            else :
                raise NameError("Unvalid PDB residue %s, must be of the form chain:residue." % desc)
    result = ""
    for chain, residues in sorted_res.items() :
        result = result + "/" + chain + ":" + ",".join(residues) + " "
    return result


# Start parsing
with open(sged_file) as csv_file:
    df = pandas.read_csv(
        csv_file, sep=delim, dtype=str, comment='#'
    )
    groups = df[group_col]
    with open(output_file, "w") as handle:
        handle.write(
                "# Created by SgedTools\nattribute: %s\nmatch mode: any\nrecipient: residues\n" % var_col
        )
        for i, g in enumerate(groups):
            sel = sged2chimera(g)
            if sel != "" :
                handle.write("\t%s\t%s\n" % (sel, df[var_col].iloc[i]))

print("Done.")
