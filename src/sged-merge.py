#! /usr/bin/python

""" Created on 24/02/20 by jdutheil

    Merge two SGED files.
"""

import getopt, sys
import pandas

cmd_args = sys.argv
arg_list = cmd_args[1:]

unix_opt = "s:t:o:g:j:c"
full_opt = ["sged1=", "sged2=", "output=", "groups=", "join=", "csv"]
try:
  arguments, values = getopt.getopt(arg_list, unix_opt, full_opt)
except getopt.error as err:
  print (str(err))
  sys.exit(2)

tabsep = True # TSV by default
group_col = "Group"
measures = []
join_type = "outer"
for arg, val in arguments:
  if arg in ("-s", "--sged1"):
    sged_file1 = val
    print "First SGED file: %s" % sged_file1
  elif arg in ("-t", "--sged2"):
    sged_file2 = val
    print "Second SGED file: %s" % sged_file2
  elif arg in ("-o", "--output"):
    output_file = val
    print "Output merged file: %s" % output_file
  elif arg in ("-g", "--groups"):
    group_col = val
    print "PDB coordinates are in column: %s" % group_col
  elif arg in ("-j", "--join"):
    join_type =  val
    print "Join type: %s" % join_type
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
  df1 = pandas.read_csv(csv_file1, sep = delim, dtype = str, comment = '#')
with open(sged_file2) as csv_file2:
  df2 = pandas.read_csv(csv_file2, sep = delim, dtype = str, comment = '#')

frames = [df1, df2]
df = pandas.merge(df1, df2, how = join_type, on = group_col)

# Write results:
df.to_csv(output_file, sep = delim, na_rep = 'NA', index = False)

print "Done."

