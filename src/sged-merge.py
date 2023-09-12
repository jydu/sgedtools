#! /usr/bin/python

""" Created on 24/02/20 by jdutheil

    Merge two SGED files.
    his program is part of the SgedTools package.

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

unix_opt = "s:t:o:g:k:j:ch"
full_opt = [
    "sged1=",
    "sged2=",
    "output=",
    "group=",
    "group1=",
    "group2=",
    "join=",
    "csv",
    "help"
]

def usage() :
    print(
"""
sged-merge

    Merge two SGED files.

Available arguments:
    --sged1 (-s): First SGED file (required).
    --sged2 (-t): Second SGED file (required).
    --group/--group1 (-g): List of columns to use for merging 
        (one or several values separated by comas, default: Group).
    --group2 (-k): List of columns to use for merging for the second file,
        if different from the first file (default: Group).
    --output (-o): Output SGED file (required).
    --join (-j): Join type to use for merging (default: outer).
        One of left, right, inner, outer, or cross.
        See https://pandas.pydata.org/docs/reference/api/pandas.merge.html for details.
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
        print("Columns to use for merging: %s" % group_col1)
    elif arg in ("-k", "--group2"):
        group_col2 = val.split(",")
        print("Columns to use for merging for second file: %s" % group_col2)
    elif arg in ("-j", "--join"):
        join_type = val
        print("Join type: %s" % join_type)
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

if not 'sged_file1' in globals():
    print("Error: a first SGED input file should be specified.")
    usage()
if not 'sged_file2' in globals():
    print("Error: a second SGED input file should be specified.")
    usage()
if not 'output_file' in globals():
    print("Error: an ouput file should be specified.")
    usage()

# Start parsing

with open(sged_file1) as csv_file1:
    df1 = pandas.read_csv(csv_file1, sep=delim, dtype=str, comment="#")
with open(sged_file2) as csv_file2:
    df2 = pandas.read_csv(csv_file2, sep=delim, dtype=str, comment="#")

frames = [df1, df2]
df = pandas.merge(df1, df2, 
                  how = join_type,
                  left_on = [x for x in group_col1],
                  right_on = [x for x in group_col2])

# Write results:
df.to_csv(output_file, sep = delim, na_rep = "NA", index = False)

print("Done.")
