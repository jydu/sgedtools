#! /usr/bin/python

""" Created on 14/02/20 by jdutheil

    Get structural info from a PDB file for groups specified in a SGED file.
"""

import getopt, sys, os.path
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

if (len(structure) > 1):
  print("Warning, %s models in PDB file %s. Using the first one." % (len(structure), pdb_file))

model = structure[0]
chain = model[chain_sel]

class ModelSelect(Select):
  def accept_model(self, model):
    if model.get_id() == 0:
      return 1
    else:
      return 0

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
      #DSSP cannot handle multiple models, we get the first one only
      pdb_file2 = pdb_file
      if len(structure) > 1:
        io = PDBIO()
        io.set_structure(structure)
        ext = pdb_file[-4:]
        if ext.lower() == ".pdb":
          pdb_file2 = pdb_file[:-4] + "_model0" + ext
        else:
          pdb_file2 = pdb_file + "_model0"
        if not os.path.isfile(pdb_file2):
          io.save(pdb_file2, ModelSelect())

      dssp = DSSP(model, pdb_file2)
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
        positions = [x[3: ] for x in res_sel_cleaned]
        states    = [x[0:3] for x in res_sel_cleaned]
        motifs    = [numpy.nan for x in positions]
        rsa       = [numpy.nan for x in positions]
        for j, pos in enumerate(positions):
          insert_code = ' '
          if len(pos) > 3:
            insert_code = pos[3:]
            pos = pos[:3]
          if (chain_sel, (' ', int(pos), insert_code)) in dssp:
            res = dssp[(chain_sel, (' ', int(pos), insert_code))] 
            states_res = states[j].title()
            if states_res in IUPACData.protein_letters_3to1:
              letter = IUPACData.protein_letters_3to1[states_res]
            else:
              letter = "X"
            if res[1] == letter:
               motifs[j] = res[2]
               if res[3] == 'NA':
                 rsa[j] = numpy.nan
               else:
                 rsa[j] = res[3]
            else:
              print "ERROR! There is no residue %s in DSSP file." % res_sel[j]
              exit(-2)
          else:
            motifs[j] = " "
            rsa[j] = numpy.nan
        results_str[i] = "".join(motifs)
        results_rsa_max[i] = numpy.nanmax(rsa) if len(rsa) > 0 and not all(numpy.isnan(rsa)) else numpy.nan
        results_rsa_min[i] = numpy.nanmin(rsa) if len(rsa) > 0 and not all(numpy.isnan(rsa)) else numpy.nan
        results_rsa_med[i] = numpy.nanmedian(rsa) if len(rsa) > 0 and not all(numpy.isnan(rsa)) else numpy.nan
        results_rsa_mea[i] = numpy.nanmean(rsa) if len(rsa) > 0 and not all(numpy.isnan(rsa)) else numpy.nan
      df["RsaMax"]    = results_rsa_max
      df["RsaMin"]    = results_rsa_min
      df["RsaMedian"] = results_rsa_med
      df["RsaMean"]   = results_rsa_mea
      df["SecondaryStructure"] = results_str

    elif measure == "DSSP": #Best for single sites:
      #DSSP cannot handle multiple models, we get the first one only
      pdb_file2 = pdb_file
      if len(structure) > 1:
        io = PDBIO()
        io.set_structure(structure)
        ext = pdb_file[-4:]
        if ext.lower() == ".pdb":
          pdb_file2 = pdb_file[:-4] + "_model0" + ext
        else:
          pdb_file2 = pdb_file + "_model0"
        if not os.path.isfile(pdb_file2):
          io.save(pdb_file2, ModelSelect())

      dssp = DSSP(model, pdb_file2)
      results_str = [numpy.nan for x in groups]
      results_rsa = [numpy.nan for x in groups]
      
      for i, g in enumerate(groups):
        tmp = g[1:(len(g)-1)]
        tmp = tmp.replace(' ', '')
        res_sel = tmp.split(";")
        # Ignore missing data:
        res_sel_cleaned = [x for x in res_sel if x != "NA"]
        positions = [x[3: ] for x in res_sel_cleaned]
        states    = [x[0:3] for x in res_sel_cleaned]
        motifs    = [numpy.nan for x in positions]
        rsa       = [numpy.nan for x in positions]
        for j, pos in enumerate(positions):
          insert_code = ' '
          if len(pos) > 3:
            insert_code = pos[3:]
            pos = pos[:3]
          if (chain_sel, (' ', int(pos), insert_code)) in dssp:
            res = dssp[(chain_sel, (' ', int(pos), insert_code))]
            states_res = states[j].title()
            if states_res in IUPACData.protein_letters_3to1:
              letter = IUPACData.protein_letters_3to1[states_res]
            else:
              letter = "X"
            if res[1] == letter:
               motifs[j] = res[2]
               if res[3] == 'NA':
                 rsa[j] = numpy.nan
               else:
                 rsa[j] = res[3]
            else:
              print "ERROR! There is no residue %s in DSSP file." % res_sel[j]
              exit(-2)
          else:
             motifs[j] = " "
             rsa[j] = numpy.nan
        results_str[i] = "".join(motifs) if len(motifs) > 0 else numpy.nan
        results_rsa[i] = numpy.nanmax(rsa) if len(rsa) > 0 and not all(numpy.isnan(rsa)) else numpy.nan
      df["Rsa"]                = results_rsa
      df["SecondaryStructure"] = results_str

    elif measure == "ResidueDepth":
      rd = ResidueDepth(model)

      results_res_depth = [numpy.nan for x in groups]
      results_ca_depth  = [numpy.nan for x in groups]
      
      for i, g in enumerate(groups):
        tmp = g[1:(len(g)-1)]
        tmp = tmp.replace(' ', '')
        res_sel = tmp.split(";")
        # Ignore missing data:
        res_sel_cleaned = [x for x in res_sel if x != "NA"]
        positions = [x[3:] for x in res_sel_cleaned]
        res_depth = [numpy.nan for x in res_sel_cleaned]
        ca_depth  = [numpy.nan for x in res_sel_cleaned]
        for j, pos in enumerate(positions):
          insert_code = ' '
          if len(pos) > 3:
            insert_code = pos[3:]
            pos = pos[:3]
          if (chain_sel, (' ', int(pos), insert_code)) in rd:
            (res_depth[j], ca_depth[j]) = rd[(chain_sel, (' ', int(pos), insert_code))]
        results_res_depth[i] = numpy.nanmean(res_depth) if len(res_depth) > 0 and not all(numpy.isnan(res_depth)) else numpy.nan
        results_ca_depth[i]  = numpy.nanmean( ca_depth) if len( ca_depth) > 0 and not all(numpy.isnan( ca_depth)) else numpy.nan
      df["ResidueDepth"] = results_res_depth
      df[ "CalphaDepth"] = results_ca_depth




  # Write results:
  df.to_csv(output_file, sep = delim, na_rep = 'NA', index = False)

print "Done."

