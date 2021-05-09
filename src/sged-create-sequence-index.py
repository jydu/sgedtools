#! /usr/bin/python

""" Created on 31/07/20 by jdutheil

    Convert alignment coordinates to species-specific coordinates
"""

import getopt, sys, glob

from Bio.SeqUtils import *
from Bio import SeqIO

cmd_args = sys.argv
arg_list = cmd_args[1:]

unix_opt = "a:r:o:f:"
full_opt = ["alignment=", "reference=", "output=", "format="]
try:
  arguments, values = getopt.getopt(arg_list, unix_opt, full_opt)
except getopt.error as err:
  print (str(err))
  sys.exit(2)

aln_format = "ig"
for arg, val in arguments:
  if arg in ("-a", "--alignment"):
    aln_file = val
    print("Alignment file: %s" % aln_file)
  elif arg in ("-r", "--reference"):
    ref_seq = val
    print("Output reference sequence: %s" % ref_seq)
  elif arg in ("-f", "--format"):
    aln_format = val
    print("Input alignment format: %s" % aln_format)
  elif arg in ("-o", "--output"):
    output_file = val
    print("Output index file: %s" % output_file)
#TODO: check that all args are provided, that file exist, and eventually allow for various alignment formats.

print("Parsing structure(s)...")

# We retrieve the original sequence from the alignment:
with open(aln_file, "r") as handle:
    aln_seqs = SeqIO.to_dict(SeqIO.parse(handle, aln_format))
aln_seq = aln_seqs[ref_seq]

print("Build the index...")

# Build the index of the sequence:
aln_index = dict()
pos = 0
for i, c in enumerate(aln_seq):
  if c != '-':
    pos = pos + 1
    aln_index[pos] = i

print("Write the results...")

with open(output_file, "w") as handle:
  handle.write("# SGED index file version 0.99\n")
  handle.write("# SGED input alignment = %s\n" % aln_file)
  handle.write("# SGED input alignment sequence = %s\n" % ref_seq)
  handle.write("# SGED index start\n")
  handle.write("AlnPos,RefRes\n")
  for seq_pos, aln_pos in aln_index.items():
      handle.write("%s,%s\n" % (aln_pos + 1, seq_pos))

print("Done.")

