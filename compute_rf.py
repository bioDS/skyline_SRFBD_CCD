#!/usr/bin/env python
"""RF distances between the CCD, MCC and true trees of one simulation folder.

Usage:
    python compute_rf.py SIM_DIR CCD_FILE MCC_FILE HELPER_DIR [TRUE_TREES_FILE TRUE_LINE]

  SIM_DIR          simulation folder, e.g. .../four_interval_simulation/inf/1
  CCD_FILE         CCD tree in SIM_DIR       (ccd_1_map_sRanges_42.tre)
  MCC_FILE         translated MCC tree       (mcc_sRanges_42_translated.tre)
  HELPER_DIR       folder holding renumber_trees.py and annotate_tree.py
  TRUE_TREES_FILE  optional: file of true trees, one Newick string per line
  TRUE_LINE        optional: line (1-based) of the true tree for this simulation

Three pairs are compared: ccd_mcc, true_ccd and true_mcc (the true-tree pairs
only when TRUE_TREES_FILE is given). The true tree is written to
SIM_DIR/true_tree.tre as Nexus so it goes through the same renumber/annotate
steps as the other two trees.

Writes SIM_DIR/rf_distance.csv (one row per pair) and prints one
tab-separated line per pair:
    RF_RESULT <pair> <status> <rf> <rf_norm> <ext_rf_clades> <ext_rf_splits>
status is "ok" when that pair was computed; missing values are NA.
"""
import csv
import os
import sys
import traceback
from copy import deepcopy
from pathlib import Path

from brokilon.core import read_nexus_trees
from brokilon.metrics import robinson_foulds_ccd, extended_robinson_foulds
from brokilon.ccd.domain.topology.ccd import get_clades

PAIRS = ("ccd_mcc", "true_ccd", "true_mcc")   # <first tree>_<second tree>
FIELDS = ["pair", "status", "rf", "rf_norm", "ext_rf_clades", "ext_rf_splits"]
TRUE_NEXUS = "true_tree.tre"


def clean(value):
    """One-line, tab-free text, so the RF_RESULT line stays parseable."""
    return " ".join(str(value).split())


def read_true_newick(path, line_no):
    with open(path) as fh:
        for n, line in enumerate(fh, start=1):
            if n == line_no:
                newick = line.strip()
                if not newick:
                    raise ValueError(f"line {line_no} of {path} is empty")
                return newick if newick.endswith(";") else newick + ";"
    raise ValueError(f"{path} has fewer than {line_no} lines")


def write_nexus(newick, out_file):
    with open(out_file, "w") as fh:
        fh.write("#NEXUS\n\nBegin trees;\n")
        fh.write(f"\ttree TRUE = {newick}\n")
        fh.write("End;\n")


def compare(tree_1, tree_2, taxon_map):
    # Work on copies: the clade/RF helpers can relabel nodes in place.
    c1, c2 = get_clades(deepcopy(tree_1)), get_clades(deepcopy(tree_2))
    rf = robinson_foulds_ccd(deepcopy(tree_1), deepcopy(tree_2))
    n_clades = len(c1) + len(c2)
    return {
        "rf": rf,
        "rf_norm": rf / n_clades if n_clades else float("nan"),
        "ext_rf_clades": extended_robinson_foulds(
            deepcopy(tree_1), deepcopy(tree_2),
            kind="sranges", taxon_map=taxon_map, use_clades=True),
        "ext_rf_splits": extended_robinson_foulds(
            deepcopy(tree_1), deepcopy(tree_2),
            kind="sranges", taxon_map=taxon_map, use_clades=False),
    }


def main():
    args = sys.argv[1:]
    if len(args) not in (4, 6):
        sys.exit(__doc__)
    sim_dir, ccd_file, mcc_file, helper_dir = args[:4]
    true_file = os.path.abspath(args[4]) if len(args) == 6 else None
    true_line = int(args[5]) if len(args) == 6 else None

    helper_dir = os.path.abspath(helper_dir)
    os.chdir(sim_dir)
    # helpers from HELPER_DIR first, then any copies in the simulation folder
    sys.path.insert(0, ".")
    sys.path.insert(0, helper_dir)
    from renumber_trees import build_taxon_map, renumber
    from annotate_tree import annotate

    # One shared numbering for ALL trees (the files order their taxa differently).
    # renumber() also repairs single-occurrence leaves mislabelled _first/_last.
    taxon_map = build_taxon_map(mcc_file)

    def load(tree_file):
        num = renumber(tree_file, tree_file.replace(".tre", "_num.tre"), taxon_map)
        # read_nexus_trees drops [&orientation=...] etc.; put them back as node attributes.
        return annotate(read_nexus_trees(Path(num))[0], num)

    trees = {"ccd": load(ccd_file), "mcc": load(mcc_file)}

    true_status = "no true trees file"
    if true_file is not None:
        try:
            write_nexus(read_true_newick(true_file, true_line), TRUE_NEXUS)
            trees["true"] = load(TRUE_NEXUS)
        except Exception as e:
            traceback.print_exc()
            true_status = clean(f"true tree failed: {e}")

    rows = []
    for pair in PAIRS:
        first, second = pair.split("_")
        row = {"pair": pair, "status": "", "rf": "NA", "rf_norm": "NA",
               "ext_rf_clades": "NA", "ext_rf_splits": "NA"}
        if first not in trees or second not in trees:
            row["status"] = true_status
        else:
            try:
                row.update(compare(trees[first], trees[second], taxon_map))
                row["status"] = "ok"
            except Exception as e:
                traceback.print_exc()
                row["status"] = clean(f"failed: {e}")
        rows.append(row)
        if row["status"] == "ok":
            print(f"{pair}: RF: {row['rf']} normalised: {row['rf_norm']} | "
                  f"extended RF clades: {row['ext_rf_clades']} "
                  f"splits: {row['ext_rf_splits']}")
        else:
            print(f"{pair}: {row['status']}")

    with open("rf_distance.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(rows)

    for row in rows:
        print("RF_RESULT\t" + "\t".join(clean(row[f]) for f in FIELDS), flush=True)


if __name__ == "__main__":
    main()