from pathlib import Path

from brokilon.core.taxon_map_utils import get_mapping_dict, change_taxon_map
import nexus_labels as nl

BASE = Path("/home/ket581/skyline/SUMMER/skyline-SRFBD-validation/CCD")

DATASETS = {
    "penguins":     ["chronostratigraphic", "constant", "equal"],
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