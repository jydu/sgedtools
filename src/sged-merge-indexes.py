#!/usr/bin/env python3

""" Created on 07/09/23 by jdutheil

    Concatenate several compatible indexes into a single one.
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

unix_opt = "i:l:o:"
full_opt = ["index=", "index-list=", "output="]

def usage() :
    print(
"""
sged-liftover-index

    Concatenate several compatible indexes into a single one.
    Each index must specify non-overlapping positions.

Available arguments:
    --index (-i): Index to merge (must be used at least twice)
    --index-lst (-l): Alternatively, a file specifying the list of indexes (one per line)
    --output (-o): Output index file (required).
    --help (-h): Print this message.
"""
    )
    sys.exit()

try:
    arguments, values = getopt.getopt(arg_list, unix_opt, full_opt)
except getopt.error as err:
    print(str(err))
    sys.exit(2)

index_list = []
for arg, val in arguments:
    if arg in ("-i", "--index"):
        index_file = val
        print("Index file: %s" % index_file)
        index_list.append(index_file)
    elif arg in ("-l", "--index-list"):
        with open(val, 'r') as handle:
            index_file = handle.readline().strip()
            print("Index file: %s" % index_file)
            index_list.append(index_file)
    elif arg in ("-o", "--output"):
        output_file = val
        print("Output index file: %s" % output_file)
    elif arg in ("-h", "--help"):
        usage()

# Check arguments:
if len(index_list) < 2:
    print("Error: at least two index files should be provided.")
    usage()

if not 'output_file' in globals():
    print("Error: an output file should be provided.")
    usage()

# Get the indexes:
pd_list = []
for index_file in index_list:
    with open(index_file, 'r') as handle:
        index = pandas.read_csv(handle, sep = ",", comment = "#", dtype = str, index_col = 0)
        pd_list.append(index)

# Now merge them:
index = pandas.concat(pd_list, verify_integrity = True) #Check that there are no duplicate

# Write the output :
with open(output_file, "w") as handle:
    handle.write("# SGED index file version 1.00\n")
    handle.write("# SGED index start\n")
    index.to_csv(handle, na_rep = "NA")

print("Done.")
