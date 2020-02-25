#! /usr/bin/python

""" Created on 24/02/20 by jdutheil

    Merge two SGED files.
"""

import getopt, sys
import pandas

cmd_args = sys.argv
arg_list = cmd_args[1:]

unix_opt = "s:t:o:g:c"
full_opt = ["sged1=", "sged2=", "output=", "groups=", "chain=", "csv"]
try:
  arguments, values = getopt.getopt(arg_list, unix_opt, full_opt)
except getopt.error as err:
  print (str(err))
  sys.exit(2)

tabsep = True # TSV by default
group_col = "Group"
measures = []
for arg, val in arguments:
  if arg in ("-s", "--sged1"):
    sged_file1 = val
    print "First SGED file: %s" % sged_file1
  elif arg in ("-t", "--sged2"):
    sged_file2 = val
    print "Second SGED file: %s" % sged_file2
  elif arg in ("-o", "--output"):
    output_file = val
    print "Output ungrouped file: %s" % output_file
  elif arg in ("-g", "--groups"):
    group_col = val
    print "PDB coordinates are in column: %s" % group_col
  elif arg in ("-c", "--csv"):
    tabsep = False

if tabsep:
  print "SGED file is in TSV format"
  delim = '\t'
else:
  print "SGED file is in CSV format"
  delim = ','

# Start parsing
with open(sged_file1) as csv_file1:
  df1 = pandas.read_csv(csv_file1, sep = delim, dtype = str)
with open(sged_file2) as csv_file2:
  df2 = pandas.read_csv(csv_file2, sep = delim, dtype = str)

frames = [df1, df2]
df = pandas.merge(df1, df2, how = 'outer', on = group_col)

# Write results:
df.to_csv(output_file, sep = delim, na_rep = 'NA', index = False)

print "Done."

