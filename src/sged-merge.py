#! /usr/bin/python

""" Created on 24/02/20 by jdutheil

    Merge two SGED files.
"""

import getopt, sys
import pandas

cmd_args = sys.argv
arg_list = cmd_args[1:]

unix_opt = "s:t:o:g:h:j:c"
full_opt = [
    "sged1=",
    "sged2=",
    "output=",
    "group=",
    "group1=",
    "group2=",
    "join=",
    "csv",
]
try:
    arguments, values = getopt.getopt(arg_list, unix_opt, full_opt)
except getopt.error as err:
    print(str(err))
    sys.exit(2)

tabsep = True  # TSV by default
group_col1 = "Group"
group_col2 = "Group"
join_type = "outer"
for arg, val in arguments:
    if arg in ("-s", "--sged1"):
        sged_file1 = val
        print("First SGED file: %s" % sged_file1)
    elif arg in ("-t", "--sged2"):
        sged_file2 = val
        print("Second SGED file: %s" % sged_file2)
    elif arg in ("-o", "--output"):
        output_file = val
        print("Output merged file: %s" % output_file)
    elif arg in (
        "-g",
        "--group",
        "--group1",
    ):  # Note: if only this arg is passed, group col name is assumed to be identical in both files
        group_col1 = val.split(",")
        group_col2 = val.split(",")
        print("Coordinates are in column: %s" % group_col1)
    elif arg in ("-h", "--group2"):
        group_col2 = val.split(",")
        print("Coordinates for second file are in column: %s" % group_col2)
    elif arg in ("-j", "--join"):
        join_type = val
        print("Join type: %s" % join_type)
    elif arg in ("-c", "--csv"):
        tabsep = False

if tabsep:
    print("SGED file is in TSV format")
    delim = "\t"
else:
    print("SGED file is in CSV format")
    delim = ","

# Start parsing
with open(sged_file1) as csv_file1:
    df1 = pandas.read_csv(csv_file1, sep=delim, dtype=str, comment="#")
with open(sged_file2) as csv_file2:
    df2 = pandas.read_csv(csv_file2, sep=delim, dtype=str, comment="#")

frames = [df1, df2]
df = pandas.merge(df1, df2, how=join_type, left_on=[x for x in group_col1], right_on=[x for x in group_col2])

# Write results:
df.to_csv(output_file, sep=delim, na_rep="NA", index=False)

print("Done.")
