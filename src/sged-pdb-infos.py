#! /usr/bin/python

""" Created on 14/02/20 by jdutheil

    Get structural info from a PDB file for groups specified in a SGED file.
"""

import getopt, sys
import pandas
import numpy
from Bio.PDB import *
from Bio.PDB.DSSP import DSSP
from Bio.SeqUtils import *

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
    results_max = [numpy.nan for x in groups]
    results_min = [numpy.nan for x in groups]
    results_med = [numpy.nan for x in groups]
    results_mea = [numpy.nan for x in groups]

    if measure == "AlphaDist":
      for i, g in enumerate(groups):
        tmp = g[1:(len(g)-1)]
        tmp = tmp.replace(' ', '')
        res_sel = tmp.split(";")
        # Ignore missing data:
        res_sel_cleaned = [x for x in res_sel if x != "NA"]
        positions = [ int(x[3:]) for x in res_sel_cleaned]
        states    = [ x[0:3] for x in res_sel_cleaned]
        calphas   = []
        for j, pos in enumerate(positions):
          if chain[pos].resname == states[j]:
            calphas.append(chain[pos]['CA'])
          else:
            print "ERROR! There is no residue %s in PDB file." % res_sel[j]
            exit(-2)
        # Compute all pairwise distances between residues CA:
        distances = []
        for j in range(1, len(calphas)):
          for k in range(j):
            distances.append(calphas[j] - calphas[k])
        results_max[i] = numpy.max(distances) if len(distances) > 0 else numpy.nan
        results_min[i] = numpy.min(distances) if len(distances) > 0 else numpy.nan
        results_med[i] = numpy.median(distances) if len(distances) > 0 else numpy.nan
        results_mea[i] = numpy.mean(distances) if len(distances) > 0 else numpy.nan
      df["AlphaDistMax"]    = results_max
      df["AlphaDistMin"]    = results_min
      df["AlphaDistMedian"] = results_med
      df["AlphaDistMean"]   = results_mea

    elif measure == "DSSPsum":
      dssp = DSSP(model, pdb_file)
      results_str = [numpy.nan for x in groups]
      results_rsa_max = [numpy.nan for x in groups]
      results_rsa_min = [numpy.nan for x in groups]
      results_rsa_med = [numpy.nan for x in groups]
      results_rsa_mea = [numpy.nan for x in groups]
      
      for i, g in enumerate(groups):
        tmp = g[1:(len(g)-1)]
        tmp = tmp.replace(' ', '')
        res_sel = tmp.split(";")
        # Ignore missing data:
        res_sel_cleaned = [x for x in res_sel if x != "NA"]
        positions = [int(x[3:]) for x in res_sel_cleaned]
        states    = [x[0:3] for x in res_sel_cleaned]
        motifs    = [numpy.nan for x in positions]
        rsa       = [numpy.nan for x in positions]
        for j, pos in enumerate(positions):
          res = dssp[(chain_sel, pos)] 
          letter = IUPACData.protein_letters_3to1[states[j].title()]
          if res[1] == letter:
             motifs[j] = res[2]
             rsa[j] = res[3]
          else:
            print "ERROR! There is no residue %s in DSSP file." % res_sel[j]
            exit(-2)
        results_str[i] = "".join(motifs)
        results_rsa_max[i] = numpy.max(rsa) if len(rsa) > 0 else numpy.nan
        results_rsa_min[i] = numpy.min(rsa) if len(rsa) > 0 else numpy.nan
        results_rsa_med[i] = numpy.median(rsa) if len(rsa) > 0 else numpy.nan
        results_rsa_mea[i] = numpy.mean(rsa) if len(rsa) > 0 else numpy.nan
      df["RsaMax"]    = results_rsa_max
      df["RsaMin"]    = results_rsa_min
      df["RsaMedian"] = results_rsa_med
      df["RsaMean"]   = results_rsa_mea
      df["SecondaryStructure"] = results_str

    elif measure == "DSSP": #Best for single sites:
      dssp = DSSP(model, pdb_file)
      results_str = [numpy.nan for x in groups]
      results_rsa = [numpy.nan for x in groups]
      
      for i, g in enumerate(groups):
        tmp = g[1:(len(g)-1)]
        tmp = tmp.replace(' ', '')
        res_sel = tmp.split(";")
        # Ignore missing data:
        res_sel_cleaned = [x for x in res_sel if x != "NA"]
        positions = [int(x[3:]) for x in res_sel_cleaned]
        states    = [x[0:3] for x in res_sel_cleaned]
        motifs    = [numpy.nan for x in positions]
        rsa       = [numpy.nan for x in positions]
        for j, pos in enumerate(positions):
          res = dssp[(chain_sel, pos)] 
          letter = IUPACData.protein_letters_3to1[states[j].title()]
          if res[1] == letter:
             motifs[j] = res[2]
             rsa[j] = res[3]
          else:
            print "ERROR! There is no residue %s in DSSP file." % res_sel[j]
            exit(-2)
        results_str[i] = "".join(motifs) if len(motifs) > 0 else numpy.nan
        results_rsa[i] = numpy.max(rsa) if len(rsa) > 0 else numpy.nan
      df["Rsa"]                = results_rsa
      df["SecondaryStructure"] = results_str




  # Write results:
  df.to_csv(output_file, sep = delim, na_rep = 'NA', index = False)

print "Done."
