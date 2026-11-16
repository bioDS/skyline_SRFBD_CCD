F=constant_with_burn_in.trees
OUT=constant_with_burn_in_no_bad.trees

grep -v -E '^tree STATE_(1091850000|1590470000|2600470000|2805510000|2972290000) ' "$F" > "$OUT"

F=equal_with_burn_in.trees
OUT=equal_with_burn_in_no_bad.trees

grep -v -E '^tree STATE_(113500000|1498600000|1965730000|2040130000|2799860000|2867760000|2882070000) ' "$F" > "$OUT"

fossirange --verbose --trees-file penguins/chronostratigraphic_with_burn_in.trees --burn-in 0 --output-file penguins/chronostratigraphic_ccd_1.tre
fossirange --verbose --trees-file penguins/constant_with_burn_in_no_bad.trees --burn-in 0 --output-file penguins/constant_ccd_1.tre
fossirange --verbose --trees-file penguins/equal_with_burn_in_no_bad.trees --burn-in 0 --output-file penguins/equal_ccd_1.tre

# fossirange --verbose --trees-file crocodilians/chronostratigraphic_with_burn_in.trees --burn-in 0 --output-file crocodilians/ccd_1_map_chronostratigraphic.tre
fossirange --verbose --trees-file crocodilians/constant_with_burn_in.trees --burn-in 0 --output-file crocodilians/constant_ccd_1.tre


# fossigraph --verbose --trees-file crocodilians/chronostratigraphic_with_burn_in.trees --burn-in 0 --output-file crocodilians/chronostratigraphic_ccd_1_graph