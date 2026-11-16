from pathlib import Path

from brokilon.core.taxon_map_utils import get_mapping_dict, change_taxon_map
import nexus_labels as nl

BASE = Path("/home/ket581/skyline/SUMMER/skyline-SRFBD-validation/CCD")

# Runs to process in each folder; edit per dataset if they differ.
DATASETS = {
    "penguins":     ["chronostratigraphic", "constant"],
    "crocodilians": ["chronostratigraphic", "constant"],
}

for folder, runs in DATASETS.items():
    d = BASE / folder
    print(f"\n=== {folder.upper()} ===")
    for run in runs:
        mcc = d / f"{run}_mcc.tre"
        ccd = d / f"{run}_ccd_1.tre"
        missing = [p.name for p in (mcc, ccd) if not p.exists()]
        if missing:
            print(f"  {run}: skipping, missing {', '.join(missing)}")
            continue

        no_internal = d / f"{run}_mcc_no_internal_labels.tre"
        renumbered  = d / f"{run}_mcc_renumbered.tre"
        formatted   = d / f"{run}_mcc_formatted.tre"

        print(f"  {run}")
        nl.strip_internal_labels(str(mcc), str(no_internal))

        new_map = get_mapping_dict(ccd)
        change_taxon_map(no_internal, renumbered, new_map)

        nl.ids_to_labels(str(renumbered), str(formatted), translate_from=str(ccd))
        nl.check(str(formatted))   # expect: taxon labels, 0 internal labels