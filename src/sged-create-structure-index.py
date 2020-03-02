#! /usr/bin/python

""" Created on 12/02/20 by jdutheil

    Convert alignment coordinates to species-specific and structure coordinates
"""

import getopt, sys, glob

from Bio.PDB import *
from Bio.SeqUtils import *
from Bio import SeqIO
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
blosum62 = matlist.blosum62
parser = PDBParser()

cmd_args = sys.argv
arg_list = cmd_args[1:]

unix_opt = "p:a:o:"
full_opt = ["pdb=", "alignment=", "output="]
try:
  arguments, values = getopt.getopt(arg_list, unix_opt, full_opt)
except getopt.error as err:
  print (str(err))
  sys.exit(2)

pdb_files = []
for arg, val in arguments:
  if arg in ("-p", "--pdb"):
    pdb_files = pdb_files + glob.glob(val)
    print "PDB file: %s" % val
  elif arg in ("-a", "--alignment"):
    aln_file = val
    print "Alignment file: %s" % aln_file
  elif arg in ("-o", "--output"):
    output_file = val
    print "Output index file: %s" % output_file
#TODO: check that all args are provided, that file exist, and eventually allow for various alignment formats.

print "Parsing structure(s)..."

pdb_seqs = dict()
for pdb_file in pdb_files:
  print("Parsing PDB file %s..." % pdb_file)
  structure = parser.get_structure('STRUCT', pdb_file)

  # First we need to check that there is only one model:
  if (len(structure) > 1):
    print("Warning, %s models in PDB file %s. Using the first one." % (len(structure), pdb_file))

  model = structure[0]

  # Check how many chain there are in the model:
  nb_chains = len(model)

  # We retrieve the sequences for each chain:
  for chain in model:
    chain_id = chain.get_id()
    chain_seq = ""
    for residue in chain:
      res = residue.get_resname().title()
      letter = ""
      if is_aa(residue, standard = True):
        letter = IUPACData.protein_letters_3to1[res]
      elif res == "Mse":
        letter = "S"
      chain_seq = chain_seq + letter
    if len(chain_seq) > 0:
      pdb_seqs[pdb_file + "|" + chain_id] = chain_seq

print "Compare structure(s) and alignment..."

# We retrieve the original sequence from the alignment:
with open(aln_file, "rU") as handle:
  aln_seqs = SeqIO.to_dict(SeqIO.parse(handle, "ig")) 

# Align each PDB sequence with each alignment sequence and get the best score
best_pdb = ""
best_aln = ""
best_score = 0
for pdb_id, pdb_seq in pdb_seqs.items():
  for aln_id, aln_seq in aln_seqs.items():
    score = pairwise2.align.globaldx(str(aln_seq.seq).replace('-', ''), pdb_seq, blosum62, score_only = True)
    if score > best_score:
      best_score = score
      best_pdb = pdb_id
      best_aln = aln_id

print "Best match between sequence %s and chain %s, with a score of %s." % (best_aln, best_pdb, best_score) 
aln_seq = aln_seqs[best_aln]
pdb_seq = pdb_seqs[best_pdb]

print "Build the index..."
(best_pdb_file, best_pdb_chain) = best_pdb.split("|")
structure = parser.get_structure('STRUCT', best_pdb_file)
model = structure[0]
nb_chains = len(model)

# Build the index of the sequence:
aln_index = dict()
pos = 0
for i, c in enumerate(aln_seq):
  if c != '-':
    pos = pos + 1
    aln_index[pos] = i

# Build the index for the PDB sequence:
chain = model[best_pdb_chain]
pdb_index = dict()
pos = 0
for residue in chain:
  res = residue.get_resname().title()
  if is_aa(residue, standard = True):
    pos = pos + 1
    letter = IUPACData.protein_letters_3to1[res]
    pdb_index[pos] = "%s%s" % (residue.get_resname(), residue.get_id()[1])
  elif res == "Mse":
    pos = pos + 1
    letter = "S"
    pdb_index[pos] = "%s%s" % (residue.get_resname(), residue.get_id()[1])

# Get the best alignment:
pairwise_aln = pairwise2.align.globaldx(str(aln_seq.seq).replace('-', ''), pdb_seq, blosum62)
print(pairwise2.format_alignment(*pairwise_aln[0]))

# Get the alignment index. If several alignments are provided, only consistent positions are kept:
def build_aln_index(aln) :
  seq1 = aln[0]
  seq2 = aln[1]
  n = len(seq1)
  pos1 = 0
  pos2 = 0
  index = dict()
  for i in range(0, n):
    if seq1[i] != '-':
      pos1 = pos1 + 1
    if seq2[i] != '-':
      pos2 = pos2 + 1
    if seq1[i] != '-' and seq2[i] != '-':
      index[pos1] = pos2
  return index

# Get index for each alignment:
indexes = dict()
count = 0
for aln in pairwise_aln:
  indexes[count] = build_aln_index(aln)
  count = count + 1

# Now get the consensus:
seq_index = dict()
for k, j in indexes[0].items():
  test = True
  for i in range(1,len(indexes)):
    if k in indexes[i]:
      if indexes[i][k] != j:
        test = False
        print "Position %s is ambiguous (2)." % k
    else:
      test = False
      print "Position %s is ambiguous (1)." % k
  if test:
    seq_index[k] = j

print "Write the results..."

# Now convert each alignment position into a PDB position and write the result to a file:
with open(output_file, "w") as handle:
  handle.write("# SGED index file version 0.99\n")
  handle.write("# SGED input alignment = %s\n" % aln_file)
  handle.write("# SGED input alignment sequence = %s\n" % best_aln)
  handle.write("# SGED input PDB = %s\n" % best_pdb_file)
  handle.write("# SGED input PDB chain = %s\n" % best_pdb_chain)
  handle.write("# SGED index start\n")
  handle.write("AlnPos,PdbRes\n")
  for seq_pos, aln_pos in aln_index.items():
    if seq_pos in seq_index :
      pdb_pos = seq_index[seq_pos]
      handle.write("%s,%s\n" % (aln_pos, pdb_index[pdb_pos]))
    else:
      handle.write("%s,NA\n" % aln_pos)

print "Done."


