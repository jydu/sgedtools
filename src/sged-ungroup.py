#!/usr/bin/env python3

""" Created on 13/02/20 by jdutheil

    Convert multi-sites groups into single sites groups. 
    Allow to specify which column to replicate.
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

import getopt, sys
import pandas

cmd_args = sys.argv
arg_list = cmd_args[1:]

unix_opt = "s:o:d:g:ch"
full_opt = ["sged=", "output=", "data=", "group=", "csv", "help"]

def usage() :
    print(
"""
sged-ungroup

    Convert multi-sites groups into single sites groups. 
    Allow to specify which column to replicate.

Available arguments:
    --sged (-s): Input SGED file (required).
    --output (-o): Output SGED file (required).
    --group (-g): Column where group coordinates are stored (default: Group).
    --data (-d): Column selection (default: empty selection). 
        Indicates which column should be added to the output file.
        Entries in the input file for the selected columns will be 
        duplicated in each entry in the output file.
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

tabsep = True  # TSV by default
selected_cols = []
group_col = "Group"
for arg, val in arguments:
    if arg in ("-s", "--sged"):
        sged_file = val
        print("SGED file: %s" % sged_file)
    elif arg in ("-o", "--output"):
        output_file = val
        print("Output ungrouped file: %s" % output_file)
    elif arg in ("-d", "--data"):
        selected_cols = val.split(",")
    elif arg in ("-g", "--group"):
        group_col = val
        print("Group coordinates are in column: %s" % group_col)
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

# Check required arguments
if not 'sged_file' in globals():
    usage()
if not 'output_file' in globals():
    usage()

# Start parsing
with open(sged_file) as csv_file:
    df = pandas.read_csv(csv_file, sep=delim, dtype=str, comment='#')
    groups = df[group_col]
    with open(output_file, "w") as handle:
        handle.write("Group%s%s\n" % (delim, delim.join(df[selected_cols].columns)))
        for i, g in enumerate(groups):
            tmp = g[1 : (len(g) - 1)]
            tmp = tmp.replace(" ", "")
            positions = tmp.split(";")
            for j in positions:
                handle.write(
                    "[%s]%s%s\n" % (j, delim, delim.join(df[selected_cols].iloc[i]))
                )

print("Done.")
