#!/usr/bin/env python3

""" Created on 12/02/20 by jdutheil
    Modified on 24/06/25 by lorenzopenone

    Convert alignment coordinates to species-specific and structure coordinates
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

import getopt, sys, glob

from Bio.PDB import *
from Bio.SeqUtils import *
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Align import substitution_matrices
from Bio.Data import PDBData
from pathlib import Path

cmd_args = sys.argv
arg_list = cmd_args[1:]

# Added -s / --pid option
unix_opt = "p:l:i:j:f:a:g:o:u:v:c:s:xh"
full_opt = [
    "pdb=",
    "pdb-list=",
    "pdb-id=",
    "pdb-id-list=",
    "pdb-format=",
    "alignment=",
    "alignment-format=",
    "output=",
    "gap-open=",
    "gap-extend=",
    "coverage=",
    "pid=",
    "exclude-incomplete",
    "help"
]

def usage():
    print(
"""
sged-create-structure-index

    Create a structure index for an alignment. Align each sequence to all chains of one 
    or more input structures and find the best match.

Available arguments:
    --pdb (-p): Input protein data bank file (required).
        Can be used multiple times to selected several entries.
        File globs can be used to select multiple files.
    --pdb-list (-l): File with list of PDB files, one per line.
        Required if --pdb not specified.
    --pdb-format (-f): Force one format for all input structures. PDB or mmCif.
        If you omit this option, the script autodetects the format from 
        each file name/extension, so you can freely pass a single list 
        containing a mix of .pdb and .cif files in the same run.
        In addition, remote:PDB or remote:mmCif  allow to directly download the 
        structure file from the Protein Data Bank. In this case, --pdb-id indicates the PDB id.
    --pdb-id (-i): Specify the id of the PDB file to retrieve remotely.
        Can be used multiple times to selected several entries.
    --pdb-id-list (-j): File with list of PDB ids to retrieve remotely.
    --alignment (-a): Input alignment file (required);
    --alignment-format (-g): Input alignment format (default: fasta).
        Any format recognized by Bio::AlignIO (see https://biopython.org/wiki/AlignIO)
    --gap-open (-u): Gap opening penalty in pairwise alignment (default: -2).
    --gap-extend (-v): Gap extension penalty in pairwise alignment (default : 0).
    --coverage (-c): Threshold used with the "free coverage" rule: max(pdb_cov, aln_cov) >= threshold.
        pdb_cov = overlap / len(pdb_seq_without_gaps)
        aln_cov = overlap / len(aln_seq_without_gaps)
    --pid (-s): Minimum percent identity (0-1) required to keep a match (default: 0.0).
        PID = identical positions / aligned positions without gaps.
    --output (-o): Output index file (required).
    --exclude-incomplete (-x): Exclude incomplete chains from scan (default: false).
    --help (-h): Print this message.
"""
    )
    sys.exit()

try:
    arguments, values = getopt.getopt(arg_list, unix_opt, full_opt)
except getopt.error as err:
    print(str(err))
    sys.exit(2)

# Helper

def load_structure(filename):
    """
    Return a Bio.PDB.Structure object regardless of file format.
    Recognises:
        *.pdb  / *.ent  → PDBParser
        *.cif / *.mmcif → MMCIFParser
    """
    suffix = Path(filename).suffix.lower()
    if suffix in {".pdb", ".ent"}:
        parser = PDBParser(QUIET=True)
    elif suffix in {".cif", ".mmcif"}:
        parser = MMCIFParser(QUIET=True)
    else:
        raise ValueError(f"Unrecognised extension for {filename}")
    # Use the file stem as the structure ID
    return parser.get_structure(Path(filename).stem, filename)


aligner = PairwiseAligner(mode='global')
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
aligner.open_gap_score = 0
aligner.extend_gap_score = 0
coverage_threshold = 0.0  # free-coverage threshold
min_pid = 0.0             # minimum percent identity (0–1)

pdb_files = []
pdb_ids = []
pdb_format = "PDB"
exclude_incomplete = False
aln_format = "fasta"
for arg, val in arguments:
    if arg in ("-p", "--pdb"):
        pdb_files = pdb_files + glob.glob(val)
        print("PDB file: %s" % val)
    elif arg in ("-l", "--pdb-list"):
        with open(val, "r") as f:
            files = [line.strip() for line in f]
        for x in files:
            print("PDB file: %s" % x)
        pdb_files = pdb_files + files
    elif arg in ("-i", "--pdb-id"):
        pdb_ids.append(val)
        print("PDB id: %s" % val)
    elif arg in ("-j", "--pdb-id-list"):
        with open(val, "r") as f:
            ids = [line.strip() for line in f]
        for x in ids:
            print("PDB id: %s" % x)
        pdb_ids = pdb_ids + ids
    elif arg in ("-f", "--pdb-format"):
        pdb_format = val
        if val.upper() != "PDB" and val.upper() != "MMCIF" and val[0:7] != "remote:":
            print(
                "Structure format should be either PDB or mmCif, or remote:PDB, remote:mmCif, etc. if you would like to retrieve the file from RCSB"
            )
            sys.exit(-1)
        print("PDB format: %s" % pdb_format)
    elif arg in ("-a", "--alignment"):
        aln_file = val
        print("Alignment file: %s" % aln_file)
    elif arg in ("-g", "--alignment-format"):
        aln_format = val
        print("Alignment format: %s" % aln_format)
    elif arg in ("-o", "--output"):
        output_file = val
        print("Output index file: %s" % output_file)
    elif arg in ("-u", "--gap-open"):
        aligner.open_gap_score = float(val)
    elif arg in ("-v", "--gap-extend"):
        aligner.extend_gap_score = float(val)
    elif arg in ("-c", "--coverage"):
        coverage_threshold = float(val)
        print(f"Coverage threshold (free coverage): {coverage_threshold}")
    elif arg in ("-s", "--pid"):
        min_pid = float(val)
        print(f"Minimum PID threshold: {min_pid}")
    elif arg in ("-x", "--exclude-incomplete"):
        exclude_incomplete = True
    elif arg in ("-h", "--help"):
        usage()

# Check options:
if len(pdb_files) == 0 and len(pdb_ids) == 0:
    print("Error: at least one structure file/id should be provided.")
    usage()

if 'aln_file' not in globals():
    print("Error: an alignment file should be provided.")
    usage()

if 'output_file' not in globals():
    print("Error: an output file should be provided.")
    usage()

print("Using alignment open gap score: %f" % aligner.open_gap_score)
print("Using alignment extend gap score: %f" % aligner.extend_gap_score)

print("Parsing structure(s)...")

if pdb_format.startswith("remote:"):
    remote_format = pdb_format[7:]
    fmt_up = remote_format.upper()
    if fmt_up == "PDB":
        remote_format = "pdb"
    elif fmt_up == "MMCIF":
        remote_format = "mmCif"
    else:
        print(f"ERROR: unsupported remote format '{remote_format}'.")
        sys.exit(1)

    pdb_server = PDBList(
        server="http://files.wwpdb.org",
        pdb=None,
        obsolete_pdb=False,
        verbose=True,
    )
    pdb_format = remote_format

if "pdb_server" in locals():
    for pdb_id in pdb_ids:
        print("Retrieving structure from remote server...")
        pdb_file = pdb_server.retrieve_pdb_file(
            pdb_code=pdb_id,
            obsolete=False,
            pdir=".",
            file_format=pdb_format,
            overwrite=False,
        )
        pdb_files.append(pdb_file)
        print("Downloaded PDB file %s..." % pdb_file)

pdb_seqs = dict()
prop_incomplete = dict()
for pdb_file in pdb_files:
    print(f"Parsing structure file {pdb_file}…")
    structure = load_structure(pdb_file)

    # First we need to check that there is only one model:
    if len(structure) > 1:
        print(
            "Warning, %s models in PDB file %s. Using the first one." % (len(structure), pdb_file)
        )

    model = structure[0]

    # Check how many chains there are in the model:
    nb_chains = len(model)

    # Retrieve sequences for each chain:
    for chain in model:
        chain_id = chain.get_id()
        chain_seq = ""
        nb_incomplete = 0
        for residue in chain:
            if len(residue) < 4:
                # Incomplete residue, most likely only Ca
                nb_incomplete += 1

            if is_aa(residue):
                res = residue.get_resname().upper()
                letter = PDBData.protein_letters_3to1_extended[res]
                if letter not in Polypeptide.d1_to_index:
                    letter = "X"
                chain_seq += letter

        if len(chain_seq) > 0:
            pdb_seqs[pdb_file + "|" + chain_id] = chain_seq
            prop_incomplete[pdb_file + "|" + chain_id] = nb_incomplete / len(chain_seq)

if exclude_incomplete and prop_incomplete:
    # Look at the proportion of incomplete data. Keep only chains with the lowest proportion.
    min_prop_incomplete = min(prop_incomplete.values())
    print("Minimum proportion of incomplete data: %s" % min_prop_incomplete)
    for seq, prop in list(prop_incomplete.items()):
        if prop > min_prop_incomplete:
            print(
                "Sequence %s has a proportion of incomplete residues equal to %s and is discarded." % (seq, prop)
            )
            del pdb_seqs[seq]

print("Compare structure(s) and alignment...")

# Retrieve the original sequence(s) from the alignment:
with open(aln_file, "r") as handle:
    aln_seqs = SeqIO.to_dict(SeqIO.parse(handle, aln_format))

# Align each PDB sequence with each alignment sequence and get the best score
best_pdb = ""
best_aln = ""
best_score = 0
best_coverage = 0.0
best_pdb_cov = 0.0
best_aln_cov = 0.0
best_pid = 0.0

for pdb_id, pdb_seq in pdb_seqs.items():
    for aln_id, aln_seq in aln_seqs.items():
        # Do the pairwise alignment
        aln = aligner.align(str(aln_seq.seq).replace("-", ""), pdb_seq)[0]
        s1, s2 = aln[0], aln[1]

        # Overlap: positions with non-gap characters in both sequences
        overlap = sum((aa1 != "-") and (aa2 != "-") for aa1, aa2 in zip(s1, s2))

        # Length of sequences without gaps
        aln_len_nogap = len(str(aln_seq.seq).replace("-", ""))
        pdb_cov = overlap / len(pdb_seq) if len(pdb_seq) > 0 else 0
        aln_cov = overlap / aln_len_nogap if aln_len_nogap > 0 else 0

        # Free coverage rule: keep the max of the two coverages
        coverage = max(pdb_cov, aln_cov)

        # Compute PID (exclude gap columns from denominator)
        aligned_cols = sum((aa1 != "-") and (aa2 != "-") for aa1, aa2 in zip(s1, s2))
        ident = sum((aa1 == aa2) and (aa1 != "-") for aa1, aa2 in zip(s1, s2))
        pid = ident / aligned_cols if aligned_cols else 0.0

        # Filters
        if coverage < coverage_threshold:
            continue
        if pid < min_pid:
            continue

        score = aln.score
        if score > best_score:
            best_score = score
            best_pdb = pdb_id
            best_aln = aln_id
            best_coverage = coverage
            best_pdb_cov = pdb_cov
            best_aln_cov = aln_cov
            best_pid = pid

if best_aln and best_pdb:
    aln_seq = aln_seqs[best_aln]
    pdb_seq = pdb_seqs[best_pdb]
    print(
        f"Best match between sequence {best_aln} and chain {best_pdb}, "
        f"score {best_score:.1f}, coverage_free {best_coverage:.1%} "
        f"(pdb_cov={best_pdb_cov:.1%}, aln_cov={best_aln_cov:.1%}), "
        f"PID={best_pid:.1%}."
    )
else:
    print(
        f"No structure satisfied the thresholds "
        f"(coverage ≥ {coverage_threshold:.1%}, PID ≥ {min_pid:.1%})."
    )
    sys.exit(1)

print("Build the index...")
(best_pdb_file, best_pdb_chain) = best_pdb.split("|")
structure = load_structure(best_pdb_file)
model = structure[0]
nb_chains = len(model)

# Build the index of the alignment sequence (alignment positions → original sequence positions):
aln_index = dict()
pos = 0
for i, c in enumerate(aln_seq):
    if c != "-":
        pos += 1
        aln_index[pos] = i

# Build the index for the PDB sequence (PDB sequence position → PDB residue string):
chain = model[best_pdb_chain]
pdb_index = dict()
pos = 0

def res_to_str(res_id):
    """Convert a residue ID tuple (hetfield, resseq, icode) to a human-readable string."""
    s = str(res_id[1])
    if res_id[2] != " ":
        s = s + res_id[2]
    return s

for residue in chain:
    if is_aa(residue):
        res = residue.get_resname().upper()
        pos += 1
        letter = PDBData.protein_letters_3to1_extended.get(res, "X")
        if letter not in Polypeptide.d1_to_index:
            letter = "X"
        pdb_index[pos] = f"{residue.get_resname()}{res_to_str(residue.get_id())}"

# Get the best alignment (BioPython may return multiple optimal alignments):
pairwise_aln = aligner.align(str(aln_seq.seq).replace("-", ""), pdb_seq)
if len(pairwise_aln) > 10:
    print(
        f"{len(pairwise_aln)} alignments returned, keeping only the 10 first ones.\n"
        "This is usually due to a too low gap penalty.\n"
        "Try rerunning with --gap-open -2 for better results.\n"
    )
    pairwise_aln = [pairwise_aln[i] for i in range(9)]

print(pairwise_aln[0])

# Build alignment index for a given pairwise alignment:
def build_aln_index(aln):
    seq1 = aln[0]
    seq2 = aln[1]
    n = len(seq1)
    pos1 = 0
    pos2 = 0
    index = dict()
    for i in range(n):
        if seq1[i] != "-":
            pos1 += 1
        if seq2[i] != "-":
            pos2 += 1
        if seq1[i] != "-" and seq2[i] != "-":
            index[pos1] = pos2
    return index

# Build index for each returned alignment:
indexes = dict()
count = 0
for aln in pairwise_aln:
    indexes[count] = build_aln_index(aln)
    count += 1

# Consensus index: keep only positions consistent across all alignments
seq_index = dict()
for k, j in indexes[0].items():
    test = True
    for i in range(1, len(indexes)):
        if k in indexes[i]:
            if indexes[i][k] != j:
                test = False
                print(f"Position {k} is ambiguous (2).")
        else:
            test = False
            print(f"Position {k} is ambiguous (1).")
    if test:
        seq_index[k] = j

print("Write the results...")

# Convert each alignment position into a PDB position and write the result to a file:
with open(output_file, "w") as handle:
    handle.write("# SGED index file version 1.00\n")
    handle.write("# SGED input alignment = %s\n" % aln_file)
    handle.write("# SGED input alignment sequence = %s\n" % best_aln)
    handle.write("# SGED input PDB = %s\n" % best_pdb_file)
    handle.write("# SGED input PDB chain = %s\n" % best_pdb_chain)
    handle.write("# SGED index start\n")
    handle.write("AlnPos,PdbRes\n")
    for seq_pos, aln_pos in aln_index.items():
        if seq_pos in seq_index:
            pdb_pos = seq_index[seq_pos]
            handle.write(f"{aln_pos + 1},{best_pdb_chain}:{pdb_index[pdb_pos]}\n")
        else:
            handle.write(f"{aln_pos + 1},NA\n")

print("Done.")
