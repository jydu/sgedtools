# Libraries
import sys, getopt
import pandas as pd
from Bio import SeqIO

""" Created on 23/05/23 by BasilePajot

    Converts the output of RASER rate shift detection to SGED.
"""

# Variable input.

cmd_args = sys.argv
arg_list = cmd_args[1:]

unix_opt = "r:a:f:o:ch"
full_opt = ["raser=",
            "alignment=",
            "alignment-format=",
            "output=", 
            "csv",
            "help"
            ]

def usage() :
    print(
"""
sged-raser2sged

    Converts the output of RASER rate shift detection to SGED.
    Output coordinates according to the original alignment position.

Available arguments:
    --raser (-t): RASER output file (required).
    --alignment (-a): RASER input alignment file (required).
    --alignment-format (-f): Input alignment format (default: fasta).
        Any format recognized by Bio::AlignIO (see https://biopython.org/wiki/AlignIO).
    --output (-o): Output SGED file (required).
    --csv (-c): Input SGED file is with comas instead of tabs (default)
    --help (-h): Print this message.
"""
    )
    sys.exit()


try:
    arguments, values = getopt.getopt(arg_list, unix_opt, full_opt)
except getopt.error as err:
    print(str(err))
    sys.exit(2)

# TSV by default
tabseq = True
aln_format = "fasta"
for arg, val in arguments:
    if arg in ("-r", "--raser"):
        raser_file = val
        print("RASER input file: %s" % raser_file)
    elif arg in ("-a", "--alignment"):
        aln_file = val
        print("Alignment file: %s" % aln_file)
    elif arg in ("-f", "--alignment-format"):
        aln_format = val
        print("Alignment format: %s" % aln_format)
    elif arg in ("-o", "--output"):
        output_file = val
        print("Output SGED file: %s" % output_file)
    elif arg in ("-c", "--csv"):
        tabseq = False
    elif arg in ("-h", "--help"):
        usage()

if tabseq:
    print("SGED file is in TSV format.")
    delim = "\t"
else:
    print("SGED file is in CSV format.")
    delim = ","

# Check required arguments

if not 'raser_file' in globals():
    print("Error: an input RASER file must be specified.")
    usage()
if not 'aln_file' in globals():
    print("Error: an input alignment file must be specified.")
    usage()
if not 'output_file' in globals():
    print("Error: an output SGED file must be specified.")
    usage()


# Read the raser file
with open(raser_file, "r") as f:
    contents = f.readlines()

# Select the name of the sequence used as reference for the raser analysis and the part of the raser file we want to use for the next steps
ref_seq_id = contents[10].split(" ")[3].replace("\n", "")

# Select the reference sequence in the alignment file
with open(aln_file, "r") as handle:
    for record in SeqIO.parse(handle, aln_format):
        if record.id == ref_seq_id :
            ref_seq = record.seq
# Transform the raser file into a sged-like dataframe
contents = contents[13:]
positive_sites = []
for line in contents:
    site_info = line.split("\t")
    positive_sites.append(site_info)
    
number_cols = 0
for line in contents:
    if line.count("\t") > number_cols:
        number_cols = line.count("\t")

number_cols += 1  

if number_cols == 3:
    raser_df = pd.DataFrame(positive_sites, columns=["position", "amino_acid", "probability"])  
else:
    raser_df = pd.DataFrame(positive_sites, columns=["position", "amino_acid", "probability", "Proba > 0.95"])
raser_df.replace(to_replace=[None], value="", inplace=True) 

if number_cols == 3:
    for counter, param in enumerate(raser_df["probability"]):
        raser_df["probability"][counter] = param[0].replace("\n", "")
else:
    for counter, param in enumerate(zip(raser_df["probability"], raser_df["Proba > 0.95"])):
        raser_df["probability"][counter] = param[0].replace("\n", "")
        raser_df["Proba > 0.95"][counter] = param[1].replace("\n", "")

# Align the sequences from the alignment file and the raser output.

raser_aa = raser_df["amino_acid"]
raser_df.drop(["amino_acid", "position"], axis=1, inplace=True)

to_write = []
counter = 0
for i, aa in enumerate(ref_seq):
     if aa == "-":
             to_write.append(["["+str(i)+"]", aa]+["NA" for _ in range(raser_df.shape[1])])
     else:
             to_write.append(["["+str(i)+"]", aa]+[raser_df.iloc[counter, x].replace("\n", "") for x in range(raser_df.shape[1])])
             counter+=1

# Put the output in a dataframe

if number_cols == 3:
    df = pd.DataFrame(to_write, columns = ["Group", "amino_acid", "probability"])
else:
    df = pd.DataFrame(to_write, columns = ["Group", "amino_acid", "probability", "Proba>0.95"])


# Export to csv

df.to_csv(output_file, index = False, sep = delim)

print("Done.")

