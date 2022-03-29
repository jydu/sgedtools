#! /usr/bin/python

""" Created on 13/02/20 by jdutheil

    Convert multi-sites groups into single sites groups. 
    Allow to specify which column to replicate.
"""

import getopt, sys
import pandas

cmd_args = sys.argv
arg_list = cmd_args[1:]

unix_opt = "s:o:d:c"
full_opt = ["sged=", "output=", "data=", "csv"]
try:
    arguments, values = getopt.getopt(arg_list, unix_opt, full_opt)
except getopt.error as err:
    print(str(err))
    sys.exit(2)

tabsep = True  # TSV by default
selected_cols = []
for arg, val in arguments:
    if arg in ("-s", "--sged"):
        sged_file = val
        print("SGED file: %s" % sged_file)
    elif arg in ("-o", "--output"):
        output_file = val
        print("Output ungrouped file: %s" % output_file)
    elif arg in ("-d", "--data"):
        selected_cols = val.split(",")
    elif arg in ("-c", "--csv"):
        tabsep = False

if tabsep:
    print("SGED file is in TSV format")
    delim = "\t"
else:
    print("SGED file is in CSV format")
    delim = ","

# Start parsing
with open(sged_file) as csv_file:
    df = pandas.read_csv(csv_file, sep=delim, dtype=str)
    groups = df["Group"]
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
