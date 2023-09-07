#! /usr/bin/python

""" Created on 07/09/23 by jdutheil

    Concatenate several compatible indexes into a single one.
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
