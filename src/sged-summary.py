#! /usr/bin/python

""" Created on 10/09/23 by jdutheil

    Summarize site statistics per group.
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
import numpy

cmd_args = sys.argv
arg_list = cmd_args[1:]

unix_opt = "s:t:o:g:k:p:f:n:ch"
full_opt = [
    "sged-groups=",
    "sged-sites=",
    "output=",
    "group=",
    "group1=",
    "group2=",
    "property=",
    "function=",
    "name=",
    "csv",
    "help"
]

def usage() :
    print(
"""
sged-merge

    Summarize site statistics per group.

Available arguments:
    --sged-groups (-s): The input groups SGED file (required).
    --sged-sites (-t): The input sites SGED file (required).
    --group/--group1 (-g): Column where group coordinates are stored (default: Group).
    --group2 (-k): Column where group coordinates are stored in the second file,
        if different from the first file (default: Group).
    --property (-p): Column in the Site file to be summarized.
        Should be specified at least once, can be used multiple times.
    --function (-f): The function to apply to all site properties for each group.
        One of min, max, mean, median.
    --output (-o): Output SGED file (required).
    --name (-n): Column name for the summary results (default: Summary).
        Should be called as many times as the --property argument.
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
summary_columns = []
properties = []
for arg, val in arguments:
    if arg in ("-s", "--sged-groups"):
        sged_groups = val
        print("Groups SGED file: %s" % sged_groups)
    elif arg in ("-t", "--sged-sites"):
        sged_sites = val
        print("Sites SGED file: %s" % sged_sites)
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
    elif arg in ("-k", "--group2"):
        group_col2 = val.split(",")
        print("Coordinates for second file are in column: %s" % group_col2)
    elif arg in ("-p", "--property"):
        properties.append(val)
        print("Summarise property: %s" % val)
    elif arg in ("-f", "--function"):
        if val == "min":
            summary_function = numpy.min
        elif val == "max":
            summary_function = numpy.max
        elif val == "median":
            summary_function = numpy.median
        elif val == "mean":
            summary_function = numpy.mean
        else:
            print("Error: unknow function.")
            usage()
        print("Summary function: %s" % val)
    elif arg in ("-n", "--name"):
        summary_columns.append(val)
        print("Summary results in: %s" % val)
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

if not 'sged_groups' in globals():
    print("Error: a Group SGED input file should be specified.")
    usage()
if not 'sged_sites' in globals():
    print("Error: a Site SGED input file should be specified.")
    usage()
if not 'output_file' in globals():
    print("Error: an ouput file should be specified.")
    usage()
if not 'summary_function' in globals():
    print("Error: a summary function should be specified.")
    usage()
if not 'summary_columns' in globals():
    print("Error: a column name for the summary should be specified.")
    usage()
if len(properties) != len(summary_columns):
    print("Error: one column name should be specified for each property.")
    usage()
if len(properties) == 0:
    print("Error: at least one property should be specified.")
    usage()


# Start parsing

with (
    open(sged_groups) as csv_file1,
    open(sged_sites) as csv_file2
):
    df_groups = pandas.read_csv(
        csv_file1, sep = delim, dtype = str, comment = "#")
    df_sites = pandas.read_csv(
        csv_file2, sep = delim, dtype = str, comment = "#",
        index_col = group_col2)
    for j,ppt in enumerate(properties):
        summary = ["NA"] * len(df_groups.index) 
        for i,g in enumerate(df_groups[group_col1]):
            sites = ["[%s]" % x for x in g[1 : (len(g) - 1)].split(";")]
            vec = [float(x) for x in df_sites.loc[sites, ppt]]
            summary[i] = summary_function(vec)
        df_groups[summary_columns[j]] = summary

# Write results:

df_groups.to_csv(output_file, sep = delim, na_rep = "NA", index = False)

print("Done.")
