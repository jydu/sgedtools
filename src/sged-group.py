#!/usr/bin/env python3

""" Created on 15/09/23 by jdutheil

    Group individual sites or groups into (super) groups.
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

unix_opt = "s:b:o:d:g:ch"
full_opt = ["sged=", "by=", "output=", "data=", "group=", "csv", "help"]

def usage() :
    print(
"""
sged-group

    Group individual sites or groups into (super) groups.
    A column can be specified to group sites.
    Other columns are discarded.

Available arguments:
    --sged (-s): Input SGED file (required).
    --by (-b): Group according to the specified column (default: group everything).
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
    elif arg in ("-b", "--by"):
        group_by = val
        print("Group sites according to column: %s" % group_by)
    elif arg in ("-o", "--output"):
        output_file = val
        print("Output grouped file: %s" % output_file)
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
    df = pandas.read_csv(csv_file, sep = delim, dtype = str, comment = '#')
    if not 'group_by' in globals():
        df['%by%'] = [True] * len(df.index)
        group_by = '%by%'
    with open(output_file, "w") as handle:
        handle.write("Group\n")
        grouped = df.groupby(group_by)
        for subset in grouped:
            groups = subset[1][group_col]
            supergroup = []
            for g in groups:
                tmp = g[1 : (len(g) - 1)]
                tmp = tmp.replace(" ", "")
                positions = tmp.split(";")
                supergroup =  supergroup + positions
            handle.write("[%s]\n" % (";".join(supergroup)))

print("Done.")
