#! /usr/bin/python

""" Created on 19/08/20 by jdutheil

    Take the groups in a SGED files and combine them in pairs.
"""

import getopt, sys
import pandas

cmd_args = sys.argv
arg_list = cmd_args[1:]

unix_opt = "s:t:o:g:h:j:c"
full_opt = ["sged=", "output=", "group=", "csv"]
try:
    arguments, values = getopt.getopt(arg_list, unix_opt, full_opt)
except getopt.error as err:
    print(str(err))
    sys.exit(2)

tabsep = True  # TSV by default
group_col = "Group"
join_type = "outer"
for arg, val in arguments:
    if arg in ("-s", "--sged"):
        sged_file = val
        print("SGED file: %s" % sged_file)
    elif arg in ("-o", "--output"):
        output_file = val
        print("Output merged file: %s" % output_file)
    elif arg in ("-g", "--group"):
        group_col = val
        print("Coordinates are in column: %s" % group_col)
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
    df = pandas.read_csv(csv_file, sep=delim, dtype=str, comment="#")

newgroups = []
for i in range(len(df) - 1):
    group1 = df.loc[i, group_col][1:-1]  # removes []
    for j in range(i + 1, len(df)):
        group2 = df.loc[j, group_col][1:-1]  # removes []
        newgroups.append("[%s;%s]" % (group1, group2))

newdf = pandas.DataFrame({group_col: newgroups})

# Write results:
newdf.to_csv(output_file, sep=delim, na_rep="NA", index=False)

print("Done.")
