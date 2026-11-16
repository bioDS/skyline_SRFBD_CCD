from copy import deepcopy
from pathlib import Path
import os
import sys

from brokilon.core import read_nexus_trees
from brokilon.metrics import robinson_foulds_ccd, extended_robinson_foulds
from brokilon.ccd.domain.topology.ccd import get_clades

os.chdir("/home/ket581/skyline/SUMMER/skyline-SRFBD-validation/CCD/crocodilians")
sys.path.insert(0, ".")  # so renumber_trees.py / annotate_tree.py (put them in this folder) is importable
from renumber_trees import build_taxon_map, renumber
from annotate_tree import annotate


import re

def find_RF(tree_1, tree_2, taxon_map):
    # Work on copies: the clade/RF helpers can relabel nodes in place, which is
    # probably where names like 'Allodaposuchus_megadontos_firstinternal1' came from.
    c1, c2 = get_clades(deepcopy(tree_1)), get_clades(deepcopy(tree_2))
    rf = robinson_foulds_ccd(deepcopy(tree_1), deepcopy(tree_2))
    print("RF:", rf, "normalised:", rf / (len(c1) + len(c2)))

    rf_2 = extended_robinson_foulds(deepcopy(tree_1), deepcopy(tree_2),
                                    kind="sranges", taxon_map=taxon_map, use_clades=True)
    print("Extended RF (sranges):", rf_2)


def run(label, ccd_file, mcc_file):
    print(label)
    # One shared numbering for BOTH trees (the two files order their taxa differently).
    taxon_map = build_taxon_map(mcc_file)
    ccd_num = renumber(ccd_file, ccd_file.replace(".tre", "_num.tre"), taxon_map)
    mcc_num = renumber(mcc_file, mcc_file.replace(".tre", "_num.tre"), taxon_map)
    # read_nexus_trees drops [&orientation=...] etc.; put them back as node attributes.
    ccd_tree = annotate(read_nexus_trees(Path(ccd_num))[0], ccd_num)
    mcc_tree = annotate(read_nexus_trees(Path(mcc_num))[0], mcc_num)
    find_RF(ccd_tree, mcc_tree, taxon_map)


run("CHRONOSTRATIGRAPHIC", "chronostratigraphic_ccd_1.tre", "chronostratigraphic_mcc_formatted.tre")
run("CONSTANT", "constant_ccd_1.tre", "constant_mcc_formatted.tre")

os.chdir("../penguins")
print("\n\nPENGUINS")
#run("CHRONOSTRATIGRAPHIC", "chronostratigraphic_ccd_1.tre", "chronostratigraphic_mcc_formatted.tre")
run("CONSTANT", "constant_ccd_1.tre", "constant_mcc_formatted.tre")
run("EQUAL", "equal_ccd_1.tre", "equal_mcc_formatted.tre")
