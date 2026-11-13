Use `Rscript find_ccd_trees_simulations.R four_interval_simulation` to find MCC and CCD MAP trees for simulated trees with four skyline intervals, or `Rscript find_ccd_trees_simulations.R one_interval_simulation` to do the same for simulated trees with one interval (so using constant diversification rates).

We apply a 10% burn-in to simulated trees, and only create the MAP trees if the posterior ESS after burn-in is at least 625.

The `crocodilians` and `penguins` folders contain the MCC and CCD MAP trees for their respective empirical datasets.
