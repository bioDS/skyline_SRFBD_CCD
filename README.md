# CCD and MCC summary trees and RF distances

This folder builds CCD MAP summary trees for the empirical datasets (`penguins`,
`crocodilians`) and for the simulated datasets. It then calculates
Robinson–Foulds (RF) distances between each CCD MAP tree and the matching MCC
tree. The empirical MCC trees come from
[bioDS/skyline-SRFBD-analyses](https://github.com/bioDS/skyline-SRFBD-analyses),
and the simulated trees are made with the setup in
[bioDS/skyline-SRFBD-validation](https://github.com/bioDS/skyline-SRFBD-validation).

## Pipeline

Run these steps in order from this folder:

```bash
# 1. Empirical datasets: CCD MAP trees
bash find_ccd_trees_penguins_crocs.sh

# 2. Simulated datasets: ESS filtering and CCD MAP trees
Rscript find_ccd_trees_simulations.R four_interval_simulation both
Rscript find_ccd_trees_simulations.R one_interval_simulation both

# 3. Empirical datasets: reformat MCC trees to match the CCD trees
python format_mcc.py

# 4. Empirical datasets: RF distances
python calculate_RF_penguins_and_crocs.py

# 5. Simulated datasets: RF distances
Rscript find_ccd_trees_simulations.R four_interval_simulation rf
Rscript find_ccd_trees_simulations.R one_interval_simulation rf
```

### Summary trees

- `find_ccd_trees_penguins_crocs.sh` builds the CCD MAP trees for the empirical
  datasets with `fossirange` (`<folder>/<run>_ccd_1.tre`). The input
  `*_with_burn_in.trees` files have already had burn-in removed, so `fossirange`
  is run with `--burn-in 0`. For penguins, a few problem trees (listed by
  `STATE_` number in the script) are first removed with `grep`, giving
  `penguins/*_with_burn_in_no_bad.trees`.
  - Penguins runs: `chronostratigraphic`, `constant`, `equal`.
  - Crocodilians runs: `chronostratigraphic`, `constant`.
- `find_ccd_trees_simulations.R <simulation> both` applies a 10% burn-in to each
  simulated posterior and, if the posterior ESS after burn-in is at least 625,
  builds the CCD MAP tree with `fossirange`. The simulated trees and their
  posteriors are made with the setup in
  [bioDS/skyline-SRFBD-validation](https://github.com/bioDS/skyline-SRFBD-validation).
  - `four_interval_simulation`: trees simulated with four skyline intervals.
  - `one_interval_simulation`: trees simulated with one interval (constant
    diversification rates).
- The empirical MCC trees (`<run>_mcc.tre`) are not made here. They come from
  [bioDS/skyline-SRFBD-analyses](https://github.com/bioDS/skyline-SRFBD-analyses).

### Formatting the empirical MCC trees

The MCC and CCD trees list their taxa in different orders, and the MCC trees
carry internal node labels. For each run in `penguins` and `crocodilians`,
`format_mcc.py`:

1. strips the internal node labels (`<run>_mcc_no_internal_labels.tre`),
2. renumbers the taxa to match the CCD tree's Translate block
   (`<run>_mcc_renumbered.tre`),
3. converts the IDs back to taxon labels (`<run>_mcc_formatted.tre`),
4. checks that the result has taxon labels and no internal labels.

Runs with a missing MCC or CCD file are skipped with a message.

### RF distances

Both RF scripts compare each CCD MAP tree with its MCC tree:

- `renumber_trees.py` relabels both trees with one shared integer taxon map
  taken from the MCC tree (`*_num.tre`). Some single-occurrence taxa appear as
  `X_last` in the CCD tree when the taxon list only has `X_first` (for example,
  *Marplesornis novaezealandiae* in the penguin constant tree). Leaves like
  these are renamed automatically, and each rename is printed.
- `annotate_tree.py` puts back the `[&orientation=...]` and `ancestral_range`
  metadata that `read_nexus_trees` drops.
- The RF distance (`robinson_foulds_ccd`), the normalised RF (RF divided by the
  total number of clades in both trees) and the extended sampled-ancestor RF
  (`extended_robinson_foulds`, `kind="sranges"`) are calculated with `brokilon`.

`calculate_RF_penguins_and_crocs.py` prints these values for all five empirical
runs: penguins `chronostratigraphic`, `constant` and `equal`, and crocodilians
`chronostratigraphic` and `constant`. `find_ccd_trees_simulations.R <simulation>
rf` runs `compute_rf.py` in each simulation folder, which writes
`rf_distance.csv` (`rf`, `rf_norm`, `ext_rf`) to that folder.

## Folders

- `penguins/`, `crocodilians/`: posterior tree files, CCD MAP trees, MCC trees,
  formatted MCC trees and renumbered trees for the empirical datasets.
- `four_interval_simulation/`, `one_interval_simulation/`: per-simulation
  summary trees and `rf_distance.csv` files.

## Requirements

R (for `find_ccd_trees_simulations.R`), `fossirange`, and Python 3.8+ with
`brokilon`. `renumber_trees.py`, `annotate_tree.py` and `nexus_labels.py` must
be importable from this folder.
