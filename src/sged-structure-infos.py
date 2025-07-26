#!/usr/bin/env python3

""" Created on 14/02/20 by jdutheil

    Get structural info from a PDB file for groups specified in a SGED file.
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

import getopt, sys, os.path
import pandas
import numpy
import scipy.cluster.hierarchy
from Bio.PDB import *
from Bio.PDB.DSSP import DSSP
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.SeqUtils import *
from Bio.Data import PDBData
from progress.bar import Bar

cmd_args = sys.argv
arg_list = cmd_args[1:]

unix_opt = "s:p:i:f:o:m:t:g:d:ch"
full_opt = [
    "sged=",
    "pdb=",
    "pdb-id=",
    "pdb-format=",
    "output=",
    "measure=",
    "threshold=",
    "group=",
    "dssp=",
    "csv",
    "help"
]

def usage() :
    print(
"""
sged-structure-infos

    Get structural info from a PDB file for groups specified in a SGED file.

Available arguments:
    --sged (-s): Input SGED file (required).
    --pdb (-p): Input protein data bank file (required).
    --pdb-format (-f): Format of the protein data bank file (default: PDB).
        Either PDB or mmCif is supported. In addition, remote:PDB or remote:mmCif
        allow to directly download the structure file from the Protein Data Bank.
        In this case, --pdb-id indicates the PDB id.
    --pdb-id (-i): Specify the id of the PDB file to retrieve remotely.
    --group (-g): Column where group coordinates are stored (default: Group).
    --measure (-m): Structural statistic to compute. Available measures are:
        * Chain: list the chains present in each group.
        * AlphaDist: 3D distance between alpha carbons of residues.
                     Min, Max, Mean and Median for all pair within each group are reported.
        * ContactSubgraphs: number of residues clusters based on a 3D distance threshold
        * ContactMap: compute the number of contacts for each residue, based on a 3D distance threshold.
        * DSSP:  DSSP measures for each site (RSA, secondary structure)
        * DSSPsum: summary for each group of DSSP measures (RSA, secondary structure)
        * ResidueDepth: distance to the surface of the protein.
        * SecondaryStructureLabel: individual label of each secondary structure unit
            (requires mmCif format as input).
    --threshold (-t): Threshold to consider for residues to be in contact (default: 8A°).
    --output (-o): Output SGED file (required).
    --result (-r): Column where to store test results (required).
    --dssp (-d): Name of the DSSP executable (default: mkdssp).
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

pdb_format = "PDB"
tabsep = True  # TSV by default
group_col = "Group"
measures = []
dssp_exe = "mkdssp"
threshold = 8
for arg, val in arguments:
    if arg in ("-s", "--sged"):
        sged_file = val
        print("SGED file: %s" % sged_file)
    elif arg in ("-p", "--pdb"):
        pdb_file = val
        print("PDB file: %s" % pdb_file)
    elif arg in ("-i", "--pdb-id"):
        pdb_id = val
        print("PDB id: %s" % val)
    elif arg in ("-f", "--pdb-format"):
        pdb_format = val
        if val != "PDB" and val != "mmCif" and val[0:7] != "remote:":
            print(
                "Structure format should be either PDB or mmCif, or remote:PDB, remote:mmCIF, etc. if you would like to retrieve the file from RCSB"
            )
            exit(-1)
        print("PDB format: %s" % pdb_format)
    elif arg in ("-o", "--output"):
        output_file = val
        print("Output info file: %s" % output_file)
    elif arg in ("-m", "--measure"):
        measures.append(val)
        print("Measure to compute: %s" % val)
    elif arg in ("-g", "--group"):
        group_col = val
        print("PDB coordinates are in column: %s" % group_col)
    elif arg in ("-t", "--threshold"):
        threshold = float(val)
        print("Contact threshold to use: %s" % threshold)
    elif arg in ("-d", "--dssp"):
        dssp_exe = val
        print("DSSP command to use: %s" % dssp_exe)
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

# Check options:

if not 'sged_file' in globals():
    print("Error: an input SGED file should be provided.")
    usage()

if not 'pdb_file' and not 'pdb_id' in globals():
    print("Error: a structure file or id should be provided.")
    usage()

if not 'output_file' in globals():
    print("Error: an output file should be provided.")
    usage()

# Parse the PDB and compute

if pdb_format.startswith("remote:"):
    remote_format = pdb_format[7:]
    pdb_server = PDBList(
        server="ftp://ftp.wwpdb.org", pdb = None, obsolete_pdb = False, verbose = True
    )
    pdb_format = remote_format

if pdb_format.upper() == "PDB":
    pdb_format = "PDB" #Case needs to be respected for remote access
    parser = PDBParser()
elif pdb_format.upper() == "MMCIF":
    pdb_format = "mmCif" #Case needs to be respected for remote access
    parser = MMCIFParser()
else:
    print("ERROR!!! Unsupported structure format: %s" % pdb_format)
    exit(-1)

if "pdb_server" in locals():
    pdb_file = pdb_server.retrieve_pdb_file(
        pdb_code = pdb_id,
        obsolete = False,
        pdir = ".",
        file_format = remote_format,
        overwrite = False,
    )
    print("Downloaded PDB file %s..." % pdb_file)

structure = parser.get_structure("STRUCT", pdb_file)


if len(structure) > 1:
    print(
        "Warning, %s models in PDB file %s. Using the first one."
        % (len(structure), pdb_file)
    )

model = structure[0]


class ModelSelect(Select):
    def accept_model(self, model):
        if model.get_id() == 0:
            return 1
        else:
            return 0


def res_to_str(id):
    s = str(id[1])
    if id[2] != " ":
        s = s + id[2]
    return s


# Start parsing
with open(sged_file) as csv_file:
    df = pandas.read_csv(csv_file, sep=delim, dtype=str, comment="#")
    groups = df[group_col]
    for measure in measures:
        pbar = Bar(measure, max = len(groups))

        if measure == "AlphaDist":

            """
            Compute Ca-Ca distances between all group members and report group statistics
            """

            results_max = [numpy.nan for x in groups]
            results_min = [numpy.nan for x in groups]
            results_med = [numpy.nan for x in groups]
            results_mea = [numpy.nan for x in groups]
            for i, g in enumerate(groups):
                if (isinstance(g, float) and numpy.isnan(g)) or str(g).strip().upper() == "NA":
                    pbar.next()
                    continue
                # strip brackets, spaces, split on “;”
                tmp = g[1:-1].replace(" ", "")
                res_sel = tmp.split(";")
                # skip row if *any* selector is NA
                if any(tok.upper() == "NA" for tok in res_sel):
                    pbar.next()
                    continue
                # keep only non‑empty, non‑NA selectors
                res_sel_cleaned = [tok for tok in res_sel if tok and tok.upper() != "NA"]
                calphas = []
                incomplete = False
                for j, pos in enumerate(res_sel_cleaned):
                    s = pos.split(":") # Assuming format Chain:Residue
                    if len(s) != 2:
                        print("ERROR! PDB coordinate should be of the form chain:residue, for instance, A:LEU123.")
                        exit(-1)

                    chain_sel = s[0]
                    pos = s[1][3:]
                    state = s[1][:3]
                    chain = model[chain_sel]
                    insert_code = " "
                    try:
                        int(pos)
                    except:
                        n = len(pos)
                        insert_code = pos[
                            (n - 1) :
                        ]  # Assuming insertion code is one character only
                        pos = pos[: (n - 1)]
                    res_id = (" ", int(pos), insert_code)
                    if not res_id in chain:
                        res_id = (
                            "H_%s" % state,
                            int(pos),
                            insert_code,
                        )  # Try with HETATM

                    if chain[res_id].resname == state:
                        if "CA" in chain[res_id]:
                            calphas.append(chain[res_id]["CA"])
                        else:
                            incomplete = True
                            print(
                                "WARNING! Residue %s has no CA. Distances cannot be computed for group %s."
                                % (chain[res_id].resname, i + 1)
                            )
                    else:
                        print("ERROR! There is no residue %s in PDB file." % res_sel[j])
                        exit(-2)
                # Compute all pairwise distances between residues CA:
                if not incomplete:
                    distances = []
                    for j in range(1, len(calphas)):
                        for k in range(j):
                            distances.append(calphas[j] - calphas[k])
                    results_max[i] = (
                        numpy.max(distances) if len(distances) > 0 else numpy.nan
                    )
                    results_min[i] = (
                        numpy.min(distances) if len(distances) > 0 else numpy.nan
                    )
                    results_med[i] = (
                        numpy.median(distances) if len(distances) > 0 else numpy.nan
                    )
                    results_mea[i] = (
                        numpy.mean(distances) if len(distances) > 0 else numpy.nan
                    )
                pbar.next()
            df["AlphaDistMax"] = results_max
            df["AlphaDistMin"] = results_min
            df["AlphaDistMedian"] = results_med
            df["AlphaDistMean"] = results_mea

        if measure == "ContactSubgraphs":

            """
            Compute the number of subgraphs based on a certain threshold distance
            """

            results_nb_subs = [numpy.nan for x in groups]
            results_nb_mapped = [numpy.nan for x in groups]
            for i, g in enumerate(groups):
                if (isinstance(g, float) and numpy.isnan(g)) or str(g).strip().upper() == "NA":
                    pbar.next()
                    continue
                # strip brackets, spaces, split on “;”
                tmp = g[1:-1].replace(" ", "")
                res_sel = tmp.split(";")
                # skip row if *any* selector is NA
                if any(tok.upper() == "NA" for tok in res_sel):
                    pbar.next()
                    continue
                # keep only non‑empty, non‑NA selectors
                res_sel_cleaned = [tok for tok in res_sel if tok and tok.upper() != "NA"]
                calphas = []
                incomplete = False
                for j, pos in enumerate(res_sel_cleaned):
                    s = pos.split(":") # Assuming format Chain:Residue
                    if len(s) != 2:
                        print("ERROR! PDB coordinate should be of the form chain:residue, for instance, A:LEU123.")
                        exit(-1)
                    chain_sel = s[0]
                    pos = s[1][3:]
                    state = s[1][:3]
                    chain = model[chain_sel]
                    insert_code = " "
                    try:
                        int(pos)
                    except:
                        n = len(pos)
                        insert_code = pos[
                            (n - 1) :
                        ]  # Assuming insertion code is one character only
                        pos = pos[: (n - 1)]
                    res_id = (" ", int(pos), insert_code)
                    if not res_id in chain:
                        res_id = (
                            "H_%s" % state,
                            int(pos),
                            insert_code,
                        )  # Try with HETATM

                    if chain[res_id].resname == state:
                        if "CA" in chain[res_id]:
                            calphas.append(chain[res_id]["CA"])
                        else:
                            incomplete = True
                            print(
                                "WARNING! Residue %s has no CA. Distances cannot be computed for group %s."
                                % (chain[res_id].resname, i + 1)
                            )
                    else:
                        print("ERROR! There is no residue %s in PDB file." % res_sel[j])
                        exit(-2)
                # Compute all pairwise distances between residues CA:
                if not incomplete and len(calphas) > 1:
                    distances = []
                    for j in range(0, len(calphas) - 1):
                        for k in range(j + 1, len(calphas)):
                            distances.append(calphas[j] - calphas[k])
                    tree = scipy.cluster.hierarchy.single(distances)
                    clusters = scipy.cluster.hierarchy.fcluster(
                        tree, threshold, criterion="distance"
                    )
                    results_nb_subs[i] = len(numpy.unique(clusters))
                    results_nb_mapped[i] = len(clusters)
                pbar.next()

            df["ContactSubgraphs.NbSubgraphs"] = results_nb_subs
            df["ContactSubgraphs.NbMapped"] = results_nb_mapped

        if measure == "ContactMap":

            """
            This provides the mean number of contacts per residue at a given thresholds
            """

            results_contact = [numpy.nan for x in groups]
            for i, g in enumerate(groups):
                if (isinstance(g, float) and numpy.isnan(g)) or str(g).strip().upper() == "NA":
                    pbar.next()
                    continue
                # strip brackets, spaces, split on “;”
                tmp = g[1:-1].replace(" ", "")
                res_sel = tmp.split(";")
                # skip row if *any* selector is NA
                if any(tok.upper() == "NA" for tok in res_sel):
                    pbar.next()
                    continue
                # keep only non‑empty, non‑NA selectors
                res_sel_cleaned = [tok for tok in res_sel if tok and tok.upper() != "NA"]
                num_contact = [0 for x in res_sel_cleaned]
                for j, pos in enumerate(res_sel_cleaned):
                    s = pos.split(":") # Assuming format Chain:Residue
                    if len(s) != 2:
                        print("ERROR! PDB coordinate should be of the form chain:residue, for instance, A:LEU123.")
                        exit(-1)
                    chain_sel = s[0]
                    pos = s[1][3:]
                    state = s[1][:3]
                    chain = model[chain_sel]
                    insert_code = " "
                    try:
                        int(pos)
                    except:
                        n = len(pos)
                        insert_code = pos[
                            (n - 1) :
                        ]  # Assuming insertion code is one character only
                        pos = pos[: (n - 1)]
                    res_id = (" ", int(pos), insert_code)
                    if not res_id in chain:
                        res_id = (
                            "H_%s" % state,
                            int(pos),
                            insert_code,
                        )  # Try with HETATM

                    if chain[res_id].resname == state:
                        # Compute distance of this residue with all others in the structure:
                        for search_chain in model:
                            for search_res in search_chain:
                                if (
                                    is_aa(search_res)
                                    and search_res.get_id() != res_id
                                    and "CA" in search_res
                                    and "CA" in chain[res_id]
                                ):
                                    distance = search_res["CA"] - chain[res_id]["CA"]
                                    if distance <= threshold:
                                        num_contact[j] = num_contact[j] + 1
                    else:
                        print("ERROR! There is no residue %s in PDB file." % res_sel[j])
                        exit(-2)
                results_contact[i] = (
                    numpy.mean(num_contact) if len(num_contact) > 0 else numpy.nan
                )
                pbar.next()

            df["NbContact"] = results_contact

        elif measure == "DSSPsum":
            # DSSP cannot handle multiple models, we get the first one only
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

            results_str = [numpy.nan for x in groups]
            results_rsa_max = [numpy.nan for x in groups]
            results_rsa_min = [numpy.nan for x in groups]
            results_rsa_med = [numpy.nan for x in groups]
            results_rsa_mea = [numpy.nan for x in groups]

            try:
                dssp = DSSP(model, pdb_file2, dssp = dssp_exe, file_type = pdb_format)

                for i, g in enumerate(groups):
                    if (isinstance(g, float) and numpy.isnan(g)) or str(g).strip().upper() == "NA":
                        pbar.next()
                        continue
                    # strip brackets, spaces, split on “;”
                    tmp = g[1:-1].replace(" ", "")
                    res_sel = tmp.split(";")
                    # skip row if *any* selector is NA
                    if any(tok.upper() == "NA" for tok in res_sel):
                        pbar.next()
                        continue
                    # keep only non‑empty, non‑NA selectors
                    res_sel_cleaned = [tok for tok in res_sel if tok and tok.upper() != "NA"]
                    motifs = [numpy.nan for x in res_sel_cleaned]
                    rsa = [numpy.nan for x in res_sel_cleaned]
                    for j, pos in enumerate(res_sel_cleaned):
                        s = pos.split(":") # Assuming format Chain:Residue
                        if len(s) != 2:
                            print("ERROR! PDB coordinate should be of the form chain:residue, for instance, A:LEU123.")
                            exit(-1)
                        chain_sel = s[0]
                        pos = s[1][3:]
                        state = s[1][:3]
                        chain = model[chain_sel]
                        insert_code = " "
                        try:
                            int(pos)
                        except:
                            n = len(pos)
                            insert_code = pos[
                                (n - 1) :
                            ]  # Assuming insertion code is one character only
                            pos = pos[: (n - 1)]
                        if (chain_sel, (" ", int(pos), insert_code)) in dssp:
                            res = dssp[(chain_sel, (" ", int(pos), insert_code))]
                            states_res = state.title()
                            letter = IUPACData.protein_letters_3to1.get(states_res, "X")
                            if res[1] == letter:
                                motifs[j] = res[2]
                                if res[3] == "NA":
                                    rsa[j] = numpy.nan
                                else:
                                    rsa[j] = res[3]
                            else:
                                print(
                                    "ERROR! There is no residue %s in DSSP file."
                                    % res_sel[j]
                                )
                                exit(-2)
                        else:
                            motifs[j] = " "
                            rsa[j] = numpy.nan
                    results_str[i] = "".join(motifs)
                    results_rsa_max[i] = (
                        numpy.nanmax(rsa)
                        if len(rsa) > 0 and not all(numpy.isnan(rsa))
                        else numpy.nan
                    )
                    results_rsa_min[i] = (
                        numpy.nanmin(rsa)
                        if len(rsa) > 0 and not all(numpy.isnan(rsa))
                        else numpy.nan
                    )
                    results_rsa_med[i] = (
                        numpy.nanmedian(rsa)
                        if len(rsa) > 0 and not all(numpy.isnan(rsa))
                        else numpy.nan
                    )
                    results_rsa_mea[i] = (
                        numpy.nanmean(rsa)
                        if len(rsa) > 0 and not all(numpy.isnan(rsa))
                        else numpy.nan
                    )
                    pbar.next()
            except Exception as e:
                print("ERROR! DSSP computation failed. Outputing 'nan'.")
                print(str(e))

            df["RsaMax"] = results_rsa_max
            df["RsaMin"] = results_rsa_min
            df["RsaMedian"] = results_rsa_med
            df["RsaMean"] = results_rsa_mea
            df["SecondaryStructure"] = results_str

        elif measure == "DSSP":  # Best for single sites:
            # DSSP cannot handle multiple models, we get the first one only
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

            results_str = [numpy.nan for x in groups]
            results_rsa = [numpy.nan for x in groups]

            try:
                dssp = DSSP(model, pdb_file2, dssp = dssp_exe, file_type = pdb_format)

                for i, g in enumerate(groups):
                    if (isinstance(g, float) and numpy.isnan(g)) or str(g).strip().upper() == "NA":
                        pbar.next()
                        continue
                    # strip brackets, spaces, split on “;”
                    tmp = g[1:-1].replace(" ", "")
                    res_sel = tmp.split(";")
                    # skip row if *any* selector is NA
                    if any(tok.upper() == "NA" for tok in res_sel):
                        pbar.next()
                        continue
                    # keep only non‑empty, non‑NA selectors
                    res_sel_cleaned = [tok for tok in res_sel if tok and tok.upper() != "NA"]
                    motifs = [numpy.nan for x in res_sel_cleaned]
                    rsa = [numpy.nan for x in res_sel_cleaned]
                    for j, pos in enumerate(res_sel_cleaned):
                        s = pos.split(":") # Assuming format Chain:Residue
                        if len(s) != 2:
                            print("ERROR! PDB coordinate should be of the form chain:residue, for instance, A:LEU123.")
                            exit(-1)
                        chain_sel = s[0]
                        pos = s[1][3:]
                        state = s[1][:3]
                        chain = model[chain_sel]
                        insert_code = " "
                        try:
                            int(pos)
                        except:
                            n = len(pos)
                            insert_code = pos[
                                (n - 1) :
                            ]  # Assuming insertion code is one character only
                            pos = pos[: (n - 1)]
                        if (chain_sel, (" ", int(pos), insert_code)) in dssp:
                            res = dssp[(chain_sel, (" ", int(pos), insert_code))]
                            states_res = state.title()
                            letter = IUPACData.protein_letters_3to1.get(states_res, "X")
                            if res[1] == letter:
                                motifs[j] = res[2]
                                if res[3] == "NA":
                                    rsa[j] = numpy.nan
                                else:
                                    rsa[j] = res[3]
                            else:
                                print(
                                    "ERROR! There is no residue %s in DSSP file."
                                    % res_sel[j]
                                )
                                exit(-2)
                        else:
                            motifs[j] = " "
                            rsa[j] = numpy.nan
                    results_str[i] = "".join(motifs) if len(motifs) > 0 else numpy.nan
                    results_rsa[i] = (
                        numpy.nanmax(rsa)
                        if len(rsa) > 0 and not all(numpy.isnan(rsa))
                        else numpy.nan
                    )
                    pbar.next()
            except Exception as e:
                print("ERROR! DSSP computation failed. Outputing 'nan'.")
                print(str(e))

            df["Rsa"] = results_rsa
            df["SecondaryStructure"] = results_str

        elif measure == "ResidueDepth":
            results_res_depth = [numpy.nan for x in groups]
            results_ca_depth = [numpy.nan for x in groups]

            try:
                rd = ResidueDepth(model)
                for i, g in enumerate(groups):
                    if (isinstance(g, float) and numpy.isnan(g)) or str(g).strip().upper() == "NA":
                        pbar.next()
                        continue
                    # strip brackets, spaces, split on “;”
                    tmp = g[1:-1].replace(" ", "")
                    res_sel = tmp.split(";")
                    # skip row if *any* selector is NA
                    if any(tok.upper() == "NA" for tok in res_sel):
                        pbar.next()
                        continue
                    # keep only non‑empty, non‑NA selectors
                    res_sel_cleaned = [tok for tok in res_sel if tok and tok.upper() != "NA"]
                    res_depth = [numpy.nan for x in res_sel_cleaned]
                    ca_depth = [numpy.nan for x in res_sel_cleaned]
                    for j, pos in enumerate(res_sel_cleaned):
                        s = pos.split(":") # Assuming format Chain:Residue
                        if len(s) != 2:
                            print("ERROR! PDB coordinate should be of the form chain:residue, for instance, A:LEU123.")
                            exit(-1)
                        chain_sel = s[0]
                        pos = s[1][3:]
                        chain = model[chain_sel]
                        insert_code = " "
                        try:
                            int(pos)
                        except:
                            n = len(pos)
                            insert_code = pos[
                                (n - 1) :
                            ]  # Assuming insertion code is one character only
                            pos = pos[: (n - 1)]
                        if (chain_sel, (" ", int(pos), insert_code)) in rd:
                            (res_depth[j], ca_depth[j]) = rd[
                                (chain_sel, (" ", int(pos), insert_code))
                            ]
                    results_res_depth[i] = (
                        numpy.nanmean(res_depth)
                        if len(res_depth) > 0 and not all(numpy.isnan(res_depth))
                        else numpy.nan
                    )
                    results_ca_depth[i] = (
                        numpy.nanmean(ca_depth)
                        if len(ca_depth) > 0 and not all(numpy.isnan(ca_depth))
                        else numpy.nan
                    )
                    pbar.next()
            except:
                print(
                    "ERROR! Computation of molecular surface failed. Outputing 'nan'."
                )

            df["ResidueDepth"] = results_res_depth
            df["CalphaDepth"] = results_ca_depth

        elif measure == "SecondaryStructureLabel":

            """
            Provide a list of secondary structure types and labels. Only works with mmCIF input files.
            """
            if pdb_format.upper() != "MMCIF":
                print("SecondaryStructureLabel only works with mmCIF input files.")
                exit(-1)

            # Load dictionary:
            mmcif_dict = MMCIF2Dict(pdb_file)
            struct_index = dict()

            # Beta-sheets
            if "_struct_sheet_range.sheet_id" in mmcif_dict:
                sheet_id = mmcif_dict["_struct_sheet_range.sheet_id"]
                range_id = mmcif_dict["_struct_sheet_range.id"]
                sta_res = mmcif_dict["_struct_sheet_range.beg_auth_comp_id"]
                sta_cha = mmcif_dict["_struct_sheet_range.beg_auth_asym_id"]
                sta_pos = mmcif_dict["_struct_sheet_range.beg_auth_seq_id"]
                end_res = mmcif_dict["_struct_sheet_range.end_auth_comp_id"]
                end_cha = mmcif_dict["_struct_sheet_range.end_auth_asym_id"]
                end_pos = mmcif_dict["_struct_sheet_range.end_auth_seq_id"]

                # Create index. First do some checks:
                nb_elts = len(sta_res)
                if len(sta_cha) != nb_elts:
                    print(
                        "ERROR! Invalid mmCIF file, wrong number of _struct_sheet_range.beg_auth_asym_id elements."
                    )
                    exit(-1)
                if len(sta_pos) != nb_elts:
                    print(
                        "ERROR! Invalid mmCIF file, wrong number of _struct_sheet_range.beg_auth_seq_id elements."
                    )
                    exit(-1)
                if len(end_res) != nb_elts:
                    print(
                        "ERROR! Invalid mmCIF file, wrong number of _struct_sheet_range.end_auth_comp_id elements."
                    )
                    exit(-1)
                if len(end_cha) != nb_elts:
                    print(
                        "ERROR! Invalid mmCIF file, wrong number of _struct_sheet_range.end_auth_asym_id elements."
                    )
                    exit(-1)
                if len(end_pos) != nb_elts:
                    print(
                        "ERROR! Invalid mmCIF file, wrong number of _struct_sheet_range.end_auth_seq_id elements."
                    )
                    exit(-1)

                for i in range(nb_elts):
                    if sta_cha[i] != end_cha[i]:
                        print(
                            "ERROR! Chain attribute of start and end positions should be identical in beta-sheet element %i. This strand annotation will be ignored."
                            % i
                        )
                    else:
                        chain_sel = sta_cha[i]
                        chain = model[chain_sel]
                        if not chain_sel in struct_index:
                            struct_index[chain_sel] = dict()
                        s = [
                            chain[x] if x in chain else None
                            for x in range(int(sta_pos[i]), int(end_pos[i]) + 1)
                        ]
                        s = list(filter(None, s))
                        for residue in s:
                            res = residue.get_resname().upper()
                            letter = PDBData.protein_letters_3to1_extended.get(res, "X")
                            if not letter in Polypeptide.d1_to_index:
                                letter = "X"
                            struct_index[chain_sel][
                                "%s%s"
                                % (residue.get_resname(), res_to_str(residue.get_id()))
                            ] = "%s-%s" % (sheet_id[i], range_id[i])

            # Helices
            if "_struct_conf.conf_type_id" in mmcif_dict:
                helix_id = mmcif_dict["_struct_conf.conf_type_id"]
                sconf_id = mmcif_dict["_struct_conf.id"]
                sta_res = mmcif_dict["_struct_conf.beg_auth_comp_id"]
                sta_cha = mmcif_dict["_struct_conf.beg_auth_asym_id"]
                sta_pos = mmcif_dict["_struct_conf.beg_auth_seq_id"]
                end_res = mmcif_dict["_struct_conf.end_auth_comp_id"]
                end_cha = mmcif_dict["_struct_conf.end_auth_asym_id"]
                end_pos = mmcif_dict["_struct_conf.end_auth_seq_id"]

                # Some checks:
                nb_elts = len(sta_res)
                if len(sta_cha) != nb_elts:
                    print(
                        "ERROR! Invalid mmCIF file, wrong number of _struct_conf.beg_auth_asym_id elements."
                    )
                    exit(-1)
                if len(sta_pos) != nb_elts:
                    print(
                        "ERROR! Invalid mmCIF file, wrong number of _struct_conf.beg_auth_seq_id elements."
                    )
                    exit(-1)
                if len(end_res) != nb_elts:
                    print(
                        "ERROR! Invalid mmCIF file, wrong number of _struct_conf.end_auth_comp_id elements."
                    )
                    exit(-1)
                if len(end_cha) != nb_elts:
                    print(
                        "ERROR! Invalid mmCIF file, wrong number of _struct_conf.end_auth_asym_id elements."
                    )
                    exit(-1)
                if len(end_pos) != nb_elts:
                    print(
                        "ERROR! Invalid mmCIF file, wrong number of _struct_conf.end_auth_seq_id elements."
                    )
                    exit(-1)

                for i in range(nb_elts):
                    if sta_cha[i] != end_cha[i]:
                        print(
                            "ERROR! Chain attribute of start and end positions should be identical in helix element %i. This helix annotation is ignored."
                            % i
                        )
                    else:
                        chain_sel = sta_cha[i]
                        chain = model[chain_sel]
                        if not chain_sel in struct_index:
                            struct_index[chain_sel] = dict()
                        s = [
                            chain[x] if x in chain else None
                            for x in range(int(sta_pos[i]), int(end_pos[i]) + 1)
                        ]
                        s = list(filter(None, s))
                        for residue in s:
                            res = residue.get_resname().upper()
                            letter = PDBData.protein_letters_3to1_extended.get(res, "X")
                            if not letter in Polypeptide.d1_to_index:
                                letter = "X"
                            struct_index[chain_sel][
                                "%s%s"
                                % (residue.get_resname(), res_to_str(residue.get_id()))
                            ] = "%s-%s" % (helix_id[i], sconf_id[i])

            # Add labels to groups:
            results_labels = [numpy.nan for x in groups]
            for i, g in enumerate(groups):
                if (isinstance(g, float) and numpy.isnan(g)) or str(g).strip().upper() == "NA":
                    pbar.next()
                    continue
                tmp = g[1:-1].replace(" ", "")
                res_sel = tmp.split(";")

                if any(tok.upper() == "NA" for tok in res_sel):
                    pbar.next()
                    continue
                res_sel = [tok for tok in res_sel if tok and tok.upper() != "NA"]
                group_labels = ["NA" for _ in res_sel]
                for j, site in enumerate(res_sel):
                    s = site.split(":") # Assuming format Chain:Residue
                    if len(s) != 2:
                        print("ERROR! PDB coordinate should be of the form chain:residue, for instance, A:LEU123.")
                        exit(-1)
                    chain_sel = s[0]
                    site = s[1]
                    chain = model[chain_sel]
                    if chain_sel in struct_index:
                        if site in struct_index[chain_sel]:
                            group_labels[j] = struct_index[chain_sel][site]
                results_labels[i] = "[%s]" % (";".join(group_labels))

                pbar.next()
                
            df["SecondaryStructureLabels"] = results_labels

        if measure == "Chain":

            """
            List the chain id of each residue.
            """

            results = [numpy.nan for x in groups]
            for i, g in enumerate(groups):
                if (isinstance(g, float) and numpy.isnan(g)) or str(g).strip().upper() == "NA":
                    pbar.next()
                    continue
                # strip brackets, spaces, split on “;”
                tmp = g[1:-1].replace(" ", "")
                res_sel = tmp.split(";")
                # skip row if *any* selector is NA
                if any(tok.upper() == "NA" for tok in res_sel):
                    pbar.next()
                    continue
                # keep only non‑empty, non‑NA selectors
                res_sel_cleaned = [tok for tok in res_sel if tok and tok.upper() != "NA"]
                chain = dict()
                for pos in res_sel_cleaned:
                    s = pos.split(":") # Assuming format Chain:Residue
                    if len(s) != 2:
                        print("ERROR! PDB coordinate should be of the form chain:residue, for instance, A:LEU123.")
                        exit(-1)

                    chain[s[0]] = 1
                if len(chain.keys()) == 0:
                    chain = "NA"
                else:
                    chain = "".join(chain.keys())
                results[i] = chain
                pbar.next()
                
            df["Chain"] = results

        print("\r")

    # Write results:
    df.to_csv(output_file, sep = delim, na_rep = "NA", index = False)

print("Done.")
