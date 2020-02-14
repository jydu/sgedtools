#! /usr/bin/python

""" Created on 14/02/20 by jdutheil

    Get structural info from a PDB file for groups specified in a SGED file.
"""

import getopt, sys
import pandas
import numpy
from Bio.PDB import *

cmd_args = sys.argv
arg_list = cmd_args[1:]

unix_opt = "s:p:o:m:g:a:c"
full_opt = ["sged=", "pdb=", "output=", "measures=", "groups=", "chain=", "csv"]
try:
  arguments, values = getopt.getopt(arg_list, unix_opt, full_opt)
except getopt.error as err:
  print (str(err))
  sys.exit(2)

tabsep = True # TSV by default
group_col = "Group"
measures = []
for arg, val in arguments:
  if arg in ("-s", "--sged"):
    sged_file = val
    print "SGED file: %s" % sged_file
  elif arg in ("-p", "--pdb"):
    pdb_file = val
    print "PDB file: %s" % pdb_file
  elif arg in ("-o", "--output"):
    output_file = val
    print "Output ungrouped file: %s" % output_file
  elif arg in ("-m", "--measures"):
    measures = val.split(',')
  elif arg in ("-g", "--groups"):
    group_col = val
    print "PDB coordinates are in column: %s" % group_col
  elif arg in ("-a", "--chain"):
    chain_sel = val
    print "PDB chain to use: %s" % chain_sel
  elif arg in ("-c", "--csv"):
    tabsep = False

if tabsep:
  print "SGED file is in TSV format"
  delim = '\t'
else:
  print "SGED file is in CSV format"
  delim = ','

# Parse the PDB and compute 
parser = PDBParser()
structure = parser.get_structure('STRUCT', pdb_file)

if (len(structure) != 1):
  print( "Error, more than one model in PDB file. Exiting.")
  exit
model = structure[0]
chain = model[chain_sel]

# Start parsing
with open(sged_file) as csv_file:
  df = pandas.read_csv(csv_file, sep = delim, dtype = str)
  groups = df[group_col]
  for measure in measures:
    results = [numpy.nan for x in groups]
    if measure == "AlphaDist":
      for i, g in enumerate(groups):
        tmp = g[1:(len(g)-1)]
        tmp = tmp.replace(' ', '')
        res_sel = tmp.split(";")
        # TODO need to account for NA: either discard the fulll group, or simply discard the site
        positions = [ int(x[3:]) for x in res_sel]
        states    = [ x[0:3] for x in res_sel]
        calphas   = []
        for j, pos in enumerate(positions):
          if chain[pos].resname == states[j]:
            calphas.append(chain[pos]['CA'])
          else:
            print "ERROR! There is no residue %s." % res_sel[j]
            exit(-2)
        # Compute all pairwise distances between residues CA:
        distances = []
        for j in range(1, len(calphas)):
          for k in range(j):
            distances.append(calphas[j] - calphas[k])
        results[i] = numpy.max(distances)
      df["AlphaDist"] = results

  # Write results:
  df.to_csv(output_file, sep = delim, na_rep = 'NA', index = False)

print "Done."

