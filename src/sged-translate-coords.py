#! /usr/bin/python

""" Created on 14/02/20 by jdutheil

    Translate coordinates according to an index.
"""

import getopt, sys
import pandas

cmd_args = sys.argv
arg_list = cmd_args[1:]

unix_opt = "s:o:i:n:r:c"
full_opt = ["sged=", "output=", "index=", "name=", "reverse", "csv"]
try:
    arguments, values = getopt.getopt(arg_list, unix_opt, full_opt)
except getopt.error as err:
    print(str(err))
    sys.exit(2)

tabsep = True  # TSV by default
reversed_index = False  # Tells if the index is inverted (second column is the key, the first one being the value). Works only if the index has 1:1 matches.
for arg, val in arguments:
    if arg in ("-s", "--sged"):
        sged_file = val
        print("SGED file: %s" % sged_file)
    elif arg in ("-o", "--output"):
        output_file = val
        print("Output translated file: %s" % output_file)
    elif arg in ("-i", "--index"):
        index_file = val
    elif arg in ("-n", "--name"):
        tln_name = val
    elif arg in ("-r", "--reverse"):
        reversed_index = True
    elif arg in ("-c", "--csv"):
        tabsep = False


if tabsep:
    print("SGED file is in TSV format.")
    delim = "\t"
else:
    print("SGED file is in CSV format.")
    delim = ","

index_col = 0
if reversed_index:
    print("Reverse index will be used.")
    index_col = 1

# Get index:
index = pandas.read_csv(
    open(index_file), sep=",", comment="#", dtype=str, index_col=index_col
)
index.index = index.index.map(str)
index.dropna(inplace=True)

# Start parsing
with open(sged_file) as csv_file:
    df = pandas.read_csv(
        csv_file, sep=delim, dtype=str, keep_default_na=False
    )  # NA in columns ignored
    groups = df["Group"]
    df.drop("Group", axis=1, inplace=True)
    with open(output_file, "w") as handle:
        handle.write(
            "Group%s%s%s%s\n" % (delim, tln_name, delim, delim.join(df.columns))
        )
        for i, g in enumerate(groups):
            tmp = g[1 : (len(g) - 1)]
            tmp = tmp.replace(" ", "")
            positions = tmp.split(";")
            translations = [
                index.loc[x][0] if x in index.index else "NA" for x in positions
            ]
            tln_group = "[%s]" % ";".join(translations)
            handle.write(
                "%s%s%s%s%s\n" % (g, delim, tln_group, delim, delim.join(df.iloc[i]))
            )

print("Done.")
