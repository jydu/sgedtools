#!/usr/bin/env python3

""" Created on 07/09/23 by jdutheil

    Combine two indexes into a new one.
    If the first index translates A into B and the second one B into C,
    the combined index will translate A into C.
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

unix_opt = "1:2:o:"
full_opt = ["index1=", "index2=", "output="]

def usage() :
    print(
"""
sged-liftover-index

    Combine two indexes into a new one.
    If the first index translates A into B and the second one B into C,
    the combined index will translate A into C.

Available arguments:
    --index1 (-1): First input index file (required).
    --index2 (-2): Second input index file (required).
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

for arg, val in arguments:
    if arg in ("-1", "--index1"):
        index1_file = val
        print("First index file: %s" % index1_file)
    elif arg in ("-2", "--index2"):
        index2_file = val
        print("Second index file: %s" % index2_file)
    elif arg in ("-o", "--output"):
        output_file = val
        print("Output index file: %s" % output_file)
    elif arg in ("-h", "--help"):
        usage()

# Check arguments:
if not 'index1_file' in globals():
    print("Error: an index file should be provided.")
    usage()

if not 'index2_file' in globals():
    print("Error: a second index file should be provided.")
    usage()

if not 'output_file' in globals():
    print("Error: an output file should be provided.")
    usage()

# Get indexes:
index1 = pandas.read_csv(open(index1_file), sep = ",", comment = "#", dtype = str)
index1.index = index1.index.map(str)
index1.dropna(inplace=True)

index2 = pandas.read_csv(open(index2_file), sep = ",", comment = "#", dtype = str, index_col = 0)
index2.index = index2.index.map(str)
index2.dropna(inplace=True)

# Parse and write the output index on the go:
with open(output_file, "w") as handle:
    handle.write("# SGED index file version 1.00\n")
    handle.write("# SGED index start\n")
    handle.write("%s,%s\n" % (index1.columns.values[0], index2.columns.values[0]))
   
    for i in range(len(index1.index)):
        first = index1.iloc[i, 0]
        x = index1.iloc[i, 1]
        if x in index2.index:
            second = index2.loc[x].iloc[0]
        else:
            second = "NA"

        handle.write(
            "%s,%s\n" % (first, second)
        )

print("Done.")
