from copy import deepcopy
from pathlib import Path
import csv
import os
import sys

from brokilon.core import read_nexus_trees
from brokilon.metrics import robinson_foulds_ccd, extended_robinson_foulds
from brokilon.ccd.domain.topology.ccd import get_clades

BASE = Path("/home/ket581/skyline/SUMMER/skyline-SRFBD-validation/CCD")
OUT_CSV = BASE / "penguin_crocs_rf_results.csv"

os.chdir(BASE / "crocodilians")
sys.path.insert(0, ".")  # so renumber_trees.py / annotate_tree.py (put them in this folder) is importable
from renumber_trees import build_taxon_map, renumber
from annotate_tree import annotate


def find_RF(tree_1, tree_2, taxon_map):
    # Work on copies: the clade/RF helpers can relabel nodes in place, which is
    # probably where names like 'Allodaposuchus_megadontos_firstinternal1' came from.
    c1, c2 = get_clades(deepcopy(tree_1)), get_clades(deepcopy(tree_2))
    rf = robinson_foulds_ccd(deepcopy(tree_1), deepcopy(tree_2))
    rf_norm = rf / (len(c1) + len(c2))
    print("RF:", rf, "normalised:", rf_norm)

    rf_clades = extended_robinson_foulds(deepcopy(tree_1), deepcopy(tree_2),
                                         kind="sranges", taxon_map=taxon_map, use_clades=True)
    print("Extended RF with clades (sranges):", rf_clades)
    rf_splits = extended_robinson_foulds(deepcopy(tree_1), deepcopy(tree_2),
                                         kind="sranges", taxon_map=taxon_map, use_clades=False)
    print("Extended RF with splits (sranges):", rf_splits)

    return {"rf": rf, "rf_normalised": rf_norm,
            "rf_clades": rf_clades, "rf_splits": rf_splits}


results = []


def run(dataset, label, ccd_file, mcc_file):
    print(label)
    # One shared numbering for BOTH trees (the two files order their taxa differently).
    taxon_map = build_taxon_map(mcc_file)
    ccd_num = renumber(ccd_file, ccd_file.replace(".tre", "_num.tre"), taxon_map)
    mcc_num = renumber(mcc_file, mcc_file.replace(".tre", "_num.tre"), taxon_map)
    # read_nexus_trees drops [&orientation=...] etc.; put them back as node attributes.
    ccd_tree = annotate(read_nexus_trees(Path(ccd_num))[0], ccd_num)
    mcc_tree = annotate(read_nexus_trees(Path(mcc_num))[0], mcc_num)
    res = find_RF(ccd_tree, mcc_tree, taxon_map)
    results.append({"dataset": dataset, "label": label, **res})


print("CROCODILIANS")
run("crocodilians", "chronostratigraphic", "chronostratigraphic_ccd_1.tre", "chronostratigraphic_mcc_formatted.tre")
run("crocodilians", "constant", "constant_ccd_1.tre", "constant_mcc_formatted.tre")

os.chdir(BASE / "penguins")
print("\n\nPENGUINS")
run("penguins", "chronostratigraphic", "chronostratigraphic_ccd_1.tre", "chronostratigraphic_mcc_formatted.tre")
run("penguins", "constant", "constant_ccd_1.tre", "constant_mcc_formatted.tre")
run("penguins", "equal", "equal_ccd_1.tre", "equal_mcc_formatted.tre")

with open(OUT_CSV, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=["dataset", "label", "rf", "rf_normalised",
                                           "rf_clades", "rf_splits"])
    writer.writeheader()
    writer.writerows(results)

print(f"\nWrote {len(results)} rows to {OUT_CSV}")