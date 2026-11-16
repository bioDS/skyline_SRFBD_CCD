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
  simulated posterior and, if the posterior ESS after burn-in is above 625,
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
`chronostratigraphic` and `constant`.

For the simulated datasets, `find_ccd_trees_simulations.R <simulation> rf` runs
`compute_rf.py` in each simulation folder with posterior ESS above 625. It
compares three pairs of trees:

| Pair       | Trees compared                              |
|------------|---------------------------------------------|
| `ccd_mcc`  | CCD MAP tree vs translated MCC tree         |
| `true_ccd` | true (simulated) tree vs CCD MAP tree       |
| `true_mcc` | true (simulated) tree vs translated MCC tree |

The true tree of simulation *i* is line *i* of `<simulation>/out_oriented`. This
file holds the simulated trees with orientation and `ancestral_range` metadata
added, and is made by `simulate.R` or `annotate_true_trees.R` in
skyline-SRFBD-validation. `compute_rf.py` saves the true tree in the simulation
folder as `true_tree.tre`.

Each simulation folder gets an `rf_distance.csv` with one row per pair
(`pair`, `status`, `rf`, `rf_norm`, `ext_rf_clades`, `ext_rf_splits`).

#### `rf_summary.csv`

`<simulation>/inf/rf_summary.csv` collects the results of all simulations in
one table. It has one row per simulation where at least one pair was
calculated, sorted by simulation number. Its columns are:

| Column                 | Description |
|------------------------|-------------|
| `sim`                  | Simulation number (the folder `inf/<sim>`). |
| `posterior_ess`        | Posterior ESS after 10% burn-in. |
| `rf_status`            | `rf run` when all three pairs were calculated. Otherwise `rf partial (...)`, with the status of each pair in brackets. |
| `rf_<pair>`            | RF distance (`robinson_foulds_ccd`). |
| `rf_norm_<pair>`       | RF divided by the total number of clades in both trees. |
| `ext_rf_clades_<pair>` | Extended sampled-ancestor RF (`kind="sranges"`) using clades. |
| `ext_rf_splits_<pair>` | Extended sampled-ancestor RF (`kind="sranges"`) using splits. |

`<pair>` is `ccd_mcc`, `true_ccd` or `true_mcc`, so there are four columns per
pair, in that order. A pair that could not be calculated has `NA` in its
columns, and the reason appears in `rf_status`. Only `rf_status` is quoted in
the file.

The same RF columns are also added to `<simulation>/inf/ess_summary.csv`. That
file has a row for every simulation, including those below the ESS cut, and
also has the tree-building status columns.

When the script is run with `test` (for example
`Rscript find_ccd_trees_simulations.R four_interval_simulation rf test`), only
simulation 1 is processed, whatever its ESS. The results go to
`rf_summary_test.csv` and `ess_summary_test.csv`, so the full summaries are
left unchanged.

## Folders

- `penguins/`, `crocodilians/`: posterior tree files, CCD MAP trees, MCC trees,
  formatted MCC trees and renumbered trees for the empirical datasets.
- `four_interval_simulation/`, `one_interval_simulation/`: the true trees
  (`out_oriented`), per-simulation summary trees and `rf_distance.csv` files
  in `inf/<sim>/`, and the `ess_summary.csv` and `rf_summary.csv` tables in
  `inf/`.

## Requirements

R (for `find_ccd_trees_simulations.R`), `fossirange`, and Python 3.8+ with
`brokilon`. `renumber_trees.py`, `annotate_tree.py` and `nexus_labels.py` must
be importable from this folder.
