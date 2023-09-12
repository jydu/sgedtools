#! /usr/bin/python

""" Created on 19/08/20 by jdutheil

    Take the groups in a SGED files and combine them in pairs.
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

unix_opt = "s:t:o:g:h:j:ch"
full_opt = ["sged=", "output=", "group=", "csv", "help"]

def usage() :
    print(
"""
sged-get-all-pairs

    Take the groups in a SGED files and combine them in pairs.

Available arguments:
    --sged (-s): Input SGED file (required).
    --group (-g): Column where group coordinates are stored (default: Group).
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

tabsep = True  # TSV by default
group_col = "Group"
for arg, val in arguments:
    if arg in ("-s", "--sged"):
        sged_file = val
        print("SGED file: %s" % sged_file)
    elif arg in ("-o", "--output"):
        output_file = val
        print("Output pairs file: %s" % output_file)
    elif arg in ("-g", "--group"):
        group_col = val
        print("Coordinates are in column: %s" % group_col)
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
    print("Error: a SGED input file should be specified.")
    usage()
if not 'output_file' in globals():
    print("Error: an ouput file should be specified.")
    usage()

# Start parsing

with open(sged_file) as csv_file:
    df = pandas.read_csv(csv_file, sep = delim, dtype = str, comment = "#")

newgroups = []
for i in range(len(df) - 1):
    group1 = df.loc[i, group_col][1:-1]  # removes []
    for j in range(i + 1, len(df)):
        group2 = df.loc[j, group_col][1:-1]  # removes []
        newgroups.append("[%s;%s]" % (group1, group2))

newdf = pandas.DataFrame({group_col: newgroups})

# Write results:
newdf.to_csv(output_file, sep = delim, na_rep = "NA", index = False)

print("Done.")
