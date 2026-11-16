######################################################
# Posterior ESS check + fossirange CCD + MCC summary + RF distance
# for sRanges runs
#
# Usage:
#   Rscript find_ccd_trees_simulations.R PARENT [MODE] [test]
#
#   PARENT: four_interval_simulation, one_interval_simulation,
#           or a full path to either
#   MODE:   ess        - compute posterior ESS only, write ess_summary.csv
#           fossirange - read ESS from PARENT/inf/ess_summary.csv and build
#                        the CCD and MCC trees where ESS > ess_cut
#           both       - compute ESS and build the trees in one pass (default)
#           rf         - read ess_summary.csv and compute the RF distances
#                        between the CCD, translated MCC and true trees where
#                        ESS > ess_cut
#           all        - ESS, trees and RF distances in one pass
#   test:   optional debugging run on simulation test_sim only (default 1):
#           its ESS is always computed, the ESS cut is ignored, and the
#           summaries go to ess_summary_test.csv / rf_summary_test.csv so the
#           full-run summaries are left untouched
#
# Examples:
#   Rscript find_ccd_trees_simulations.R four_interval_simulation ess
#   Rscript find_ccd_trees_simulations.R four_interval_simulation fossirange
#   Rscript find_ccd_trees_simulations.R four_interval_simulation rf
#   Rscript find_ccd_trees_simulations.R four_interval_simulation all test
#
# ESS mode, for every PARENT/inf/<i>/sRanges.42.log:
#   - discards the first 10% of samples as burn-in
#   - computes the effective sample size of the posterior column
#   - writes it to PARENT/inf/<i>/posterior_ess.txt
# Tree building is run inside PARENT/inf/<i>/ as (existing outputs are
# overwritten):
#   fossirange --verbose --trees-file sRanges.42.trees --burn-in 0.1
#              --output-file ccd_1_map_sRanges_42.tre
#   ~/beast/bin/applauncher SRTreeAnnotator -trees sRanges.42.trees
#              -out mcc_sRanges_42.tre -burnin 10
#   python translate_tree.py mcc_sRanges_42.tre mcc_sRanges_42_translated.tre
#              (only when the MCC tree was built)
# RF mode, for each eligible PARENT/inf/<i>/:
#   python compute_rf.py <sim_dir> ccd_1_map_sRanges_42.tre
#              mcc_sRanges_42_translated.tre <helper_dir> <true_trees> <line>
#   - three pairs: ccd_mcc (CCD vs MCC), true_ccd (true vs CCD) and
#     true_mcc (true vs MCC)
#   - the true tree of simulation i is line (i - true_first_sim + 1) of the
#     true trees file (PARENT/out: one Newick string per line). The script
#     stops if that file can't be found.
#   - per pair: RF (robinson_foulds_ccd), RF normalised by the total number
#     of clades in both trees, and extended RF (kind = "sranges") with
#     clades and with splits
#   - writes PARENT/inf/<i>/rf_distance.csv (one row per pair)
#   - needs brokilon in the python on PATH, and renumber_trees.py and
#     annotate_tree.py in helper_dir
# The summary table PARENT/inf/ess_summary.csv
# (sim, posterior_ess, status, ccd_status, mcc_status, translate_status,
#  rf_status, then rf_<pair>, rf_norm_<pair>, ext_rf_clades_<pair>,
#  ext_rf_splits_<pair> for each pair)
# is written by every mode. The rf and all modes also write
# PARENT/inf/rf_summary.csv (sim, posterior_ess, rf_status and the RF
# columns), one row per simulation where at least one pair was computed.
######################################################
library(coda)

base_dir  <- "/home/ket581/skyline/SUMMER/skyline-SRFBD-validation/CCD/"
allowed   <- c("four_interval_simulation", "one_interval_simulation")
modes     <- c("ess", "fossirange", "both", "rf", "all")
burnin    <- 0.1
ess_cut   <- 625
log_name  <- "sRanges.42.log"
tree_name <- "sRanges.42.trees"
out_name  <- "ccd_1_map_sRanges_42.tre"
mcc_bin   <- path.expand("~/beast/bin/applauncher")   # runs SRTreeAnnotator from sRanges
mcc_name  <- "mcc_sRanges_42.tre"
# translate_tree.py: set this to where the script lives
translate_script <- file.path(base_dir, "translate_tree.py")
translated_name  <- "mcc_sRanges_42_translated.tre"
python_bin       <- "python"   # the python on PATH (e.g. your conda env)
# RF distance: compute_rf.py, and the folder with renumber_trees.py and
# annotate_tree.py (set both to where the scripts live)
rf_script   <- file.path(base_dir, "compute_rf.py")
helper_dir  <- base_dir
rf_out_name <- "rf_distance.csv"
# True trees: one Newick string per line, relative to PARENT. May be the
# file itself or a folder holding exactly one file.
true_trees_path <- file.path("out_oriented")
true_first_sim  <- 1   # simulation number whose true tree is on line 1
test_sim        <- 1   # simulation processed when run with the 'test' argument

# tree pairs compared by compute_rf.py, and the columns stored for each
rf_pairs     <- c("ccd_mcc", "true_ccd", "true_mcc")
rf_fields    <- c("rf", "rf_norm", "ext_rf_clades", "ext_rf_splits")
rf_cols      <- as.vector(outer(rf_fields, rf_pairs, paste, sep = "_"))
rf_num_cols  <- rf_cols[!startsWith(rf_cols, "ext_")]
rf_chr_cols  <- rf_cols[startsWith(rf_cols, "ext_")]
legacy_cols  <- c("rf", "rf_norm", "ext_rf")   # from the CCD-vs-MCC-only version

######################################################
# Arguments and paths
######################################################
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1)
  stop("Usage: Rscript find_ccd_trees_simulations.R PARENT [",
       paste(modes, collapse = "|"), "] [test]")

mode <- if (length(args) >= 2) tolower(args[2]) else "both"
if (!mode %in% modes)
  stop("MODE must be one of: ", paste(modes, collapse = ", "))
test_mode <- length(args) >= 3 && tolower(args[3]) == "test"
if (length(args) >= 3 && !test_mode)
  stop("Third argument must be 'test' (got '", args[3], "')")

do_ess   <- mode %in% c("ess", "both", "all")
do_trees <- mode %in% c("fossirange", "both", "all")
do_rf    <- mode %in% c("rf", "all")

parent <- sub("/+$", "", args[1])
if (!basename(parent) %in% allowed)
  stop("Parent folder must be one of: ", paste(allowed, collapse = ", "))
# accept either a path that exists as given or a name relative to base_dir
if (!dir.exists(parent)) parent <- paste0(base_dir, basename(parent))
parent  <- normalizePath(parent, mustWork = TRUE)   # absolute, so setwd can't break it
inf_dir <- file.path(parent, "inf")
if (!dir.exists(inf_dir)) stop("No inference folder at ", inf_dir)
suffix <- if (test_mode) "_test" else ""
summary_file    <- file.path(inf_dir, paste0("ess_summary", suffix, ".csv"))
rf_summary_file <- file.path(inf_dir, paste0("rf_summary", suffix, ".csv"))

if (do_trees || do_rf) {
  if (Sys.which(python_bin) == "")
    stop(python_bin, " not found on PATH")
}
if (do_trees) {
  if (Sys.which("fossirange") == "")
    stop("fossirange not found on PATH")
  if (!file.exists(mcc_bin) || file.access(mcc_bin, 1) != 0)
    stop("applauncher not found or not executable at ", mcc_bin)
  if (!file.exists(translate_script))
    stop("translate_tree.py not found at ", translate_script)
  translate_script <- normalizePath(translate_script)   # absolute, used after setwd
}

# the true trees file; stops if it can't be found
find_true_trees <- function() {
  p <- file.path(parent, true_trees_path)
  if (dir.exists(p)) {
    f <- list.files(p, full.names = TRUE)
    f <- f[!dir.exists(f)]
    if (length(f) == 1) return(normalizePath(f))
    stop("Expected one true trees file in ", p, ", found ", length(f),
         " - set true_trees_path to the file")
  }
  if (file.exists(p)) return(normalizePath(p))
  stop("No true trees file at ", p)
}

true_trees_file <- NA_character_
if (do_rf) {
  if (!file.exists(rf_script))
    stop("compute_rf.py not found at ", rf_script)
  for (h in c("renumber_trees.py", "annotate_tree.py"))
    if (!file.exists(file.path(helper_dir, h)))
      stop(h, " not found in ", helper_dir)
  rf_script  <- normalizePath(rf_script)
  helper_dir <- normalizePath(helper_dir)
  true_trees_file <- find_true_trees()
}

cat(basename(parent), "- mode:", mode,
    if (test_mode) paste0("(TEST: simulation ", test_sim, " only)"), "\n")

######################################################
# Helpers
######################################################
# posterior ESS after burn-in for one simulation folder
compute_ess <- function(sim_dir) {
  log_file <- file.path(sim_dir, log_name)
  if (!file.exists(log_file)) return(list(ess = NA_real_, status = "no log"))

  t <- tryCatch(read.table(log_file, header = TRUE, sep = "\t"),
                error = function(e) NULL)
  if (is.null(t) || !"posterior" %in% names(t))
    return(list(ess = NA_real_, status = "unreadable log"))

  t <- t[-seq_len(ceiling(nrow(t) * burnin)), , drop = FALSE]
  if (nrow(t) < 2) return(list(ess = NA_real_, status = "too few samples"))

  ess <- as.numeric(effectiveSize(as.mcmc(t$posterior)))
  writeLines(format(ess, digits = 7), file.path(sim_dir, "posterior_ess.txt"))
  list(ess = ess, status = if (ess > ess_cut) "ESS ok" else "ESS too low")
}

# run a command inside sim_dir; the working directory is always restored.
# Any existing output file is deleted first, so it is always overwritten.
# Returns "<label> run", "no <input>" or "<label> failed (exit n)".
run_in_dir <- function(sim_dir, input, output, cmd, cmd_args, label) {
  if (!file.exists(file.path(sim_dir, input))) return(paste("no", input))
  out_file <- file.path(sim_dir, output)
  if (file.exists(out_file) && !file.remove(out_file))
    return(paste0(label, " failed (could not remove old ", output, ")"))
  old_wd <- setwd(sim_dir)
  on.exit(setwd(old_wd), add = TRUE)
  exit <- tryCatch(system2(cmd, args = cmd_args),
                   error = function(e) { message(conditionMessage(e)); -1 })
  if (identical(as.integer(exit), 0L)) paste(label, "run")
  else paste0(label, " failed (exit ", exit, ")")
}

run_fossirange <- function(sim_dir)
  run_in_dir(sim_dir, tree_name, out_name, "fossirange",
             c("--verbose",
               "--trees-file", tree_name,
               "--burn-in", "0.1",
               "--output-file", out_name),
             "fossirange")

run_mcc <- function(sim_dir)
  run_in_dir(sim_dir, tree_name, mcc_name, mcc_bin,
             c("SRTreeAnnotator",
               "-trees", tree_name,
               "-out", mcc_name,
               "-burnin", "10"),
             "mcc")

# replace numeric taxon IDs in the MCC tree with names (Nexus, metadata kept)
run_translate <- function(sim_dir)
  run_in_dir(sim_dir, mcc_name, translated_name, python_bin,
             c(shQuote(translate_script), mcc_name, translated_name),
             "translate")

# empty RF values: NA for the numeric columns, "" for the extended RF ones
empty_rf_values <- function() {
  v <- vector("list", length(rf_cols))
  names(v) <- rf_cols
  for (col in rf_num_cols) v[[col]] <- NA_real_
  for (col in rf_chr_cols) v[[col]] <- ""
  v
}

# add any missing RF columns (and drop those of the older version)
init_rf_cols <- function(df) {
  df[intersect(legacy_cols, names(df))] <- NULL
  for (col in rf_num_cols) {
    if (!col %in% names(df)) df[[col]] <- NA_real_
    df[[col]] <- suppressWarnings(as.numeric(df[[col]]))
  }
  for (col in rf_chr_cols) {
    if (!col %in% names(df)) df[[col]] <- ""
    df[[col]] <- ifelse(is.na(df[[col]]), "", as.character(df[[col]]))
  }
  df
}

# write a summary CSV with quotes only around rf_status (and around any other
# text column that holds a comma, quote or newline, which would otherwise
# break the CSV). Extended RF values that are numbers are written as numbers.
write_summary <- function(df, file) {
  for (col in intersect(rf_chr_cols, names(df))) {
    x   <- df[[col]]
    num <- suppressWarnings(as.numeric(x))
    if (all(is.na(num) == (is.na(x) | x == ""))) df[[col]] <- num
  }
  needs_quotes <- vapply(df, function(x) is.character(x) && any(grepl('[,"\n]', x)),
                         logical(1))
  q <- which(names(df) == "rf_status" | needs_quotes)
  writeLines(paste(names(df), collapse = ","), file)   # header without quotes
  write.table(df, file, sep = ",", row.names = FALSE, col.names = FALSE,
              quote = if (length(q)) q else FALSE, qmethod = "double",
              append = TRUE)
}

# RF distances between the CCD, translated MCC and true trees of one folder.
# Returns list(values = named list over rf_cols, status).
run_rf <- function(sim_dir, sim) {
  res <- list(values = empty_rf_values(), status = "")
  for (f in c(out_name, translated_name))
    if (!file.exists(file.path(sim_dir, f))) {
      res$status <- paste("no", f)
      return(res)
    }
  out_file <- file.path(sim_dir, rf_out_name)
  if (file.exists(out_file) && !file.remove(out_file)) {
    res$status <- paste0("rf failed (could not remove old ", rf_out_name, ")")
    return(res)
  }

  py_args <- c(shQuote(rf_script), shQuote(sim_dir),
               out_name, translated_name, shQuote(helper_dir))
  true_line <- sim - true_first_sim + 1
  if (!is.na(true_trees_file) && true_line >= 1)
    py_args <- c(py_args, shQuote(true_trees_file), true_line)

  out <- tryCatch(
    suppressWarnings(system2(python_bin, py_args, stdout = TRUE, stderr = "")),
    error = function(e) { message(conditionMessage(e)); NULL })
  if (is.null(out)) { res$status <- "rf failed (could not run python)"; return(res) }

  # echo the script's own output, minus the machine-readable lines
  is_result <- startsWith(out, "RF_RESULT\t")
  if (any(!is_result)) cat(paste0("    ", out[!is_result]), sep = "\n")

  exit <- attr(out, "status")
  if ((!is.null(exit) && exit != 0) || !any(is_result)) {
    res$status <- paste0("rf failed (exit ", if (is.null(exit)) 0 else exit, ")")
    return(res)
  }

  # RF_RESULT <pair> <status> <rf> <rf_norm> <ext_rf_clades> <ext_rf_splits>
  pair_status <- setNames(rep("not reported", length(rf_pairs)), rf_pairs)
  for (line in out[is_result]) {
    f <- strsplit(line, "\t", fixed = TRUE)[[1]]
    if (length(f) < 7 || !f[2] %in% rf_pairs) next
    p <- f[2]
    pair_status[p] <- f[3]
    if (f[3] != "ok") next
    res$values[[paste0("rf_", p)]]            <- suppressWarnings(as.numeric(f[4]))
    res$values[[paste0("rf_norm_", p)]]       <- suppressWarnings(as.numeric(f[5]))
    res$values[[paste0("ext_rf_clades_", p)]] <- f[6]
    res$values[[paste0("ext_rf_splits_", p)]] <- f[7]
  }

  n_ok   <- sum(pair_status == "ok")
  detail <- paste(paste0(names(pair_status), ": ", pair_status), collapse = "; ")
  res$status <- if (n_ok == length(rf_pairs)) "rf run"
                else if (n_ok > 0) paste0("rf partial (", detail, ")")
                else paste0("rf failed (", detail, ")")
  res
}

######################################################
# Get ESS: compute it (ess / both / all, or any mode in test) or read it
# back (fossirange / rf)
######################################################
if (do_ess || test_mode) {
  sims <- list.dirs(inf_dir, full.names = FALSE, recursive = FALSE)
  sims <- sims[grepl("^[0-9]+$", sims)]
  sims <- sims[order(as.integer(sims))]
  if (test_mode) {
    sims <- sims[as.integer(sims) == test_sim]
    if (length(sims) == 0)
      stop("No simulation folder ", file.path(inf_dir, test_sim))
  }
  cat(length(sims), "simulation folder(s)\n")

  summary_df <- do.call(rbind, lapply(sims, function(i) {
    res <- compute_ess(file.path(inf_dir, i))
    cat("sim", i, ": ESS =", round(res$ess, 1), "-", res$status, "\n")
    data.frame(sim = as.integer(i), posterior_ess = res$ess, status = res$status)
  }))
  summary_df$ccd_status <- ""
  summary_df$mcc_status <- ""
  summary_df$translate_status <- ""
  summary_df$rf_status <- ""
  summary_df <- init_rf_cols(summary_df)
  write_summary(summary_df, summary_file)
  cat("ESS summary written to", summary_file, "\n")
} else {
  if (!file.exists(summary_file))
    stop("No ", summary_file, " - run with mode 'ess' first")
  summary_df <- read.csv(summary_file, stringsAsFactors = FALSE)
  if (!all(c("sim", "posterior_ess") %in% names(summary_df)))
    stop(summary_file, " needs columns sim and posterior_ess")
  for (col in c("status", "ccd_status", "mcc_status", "translate_status",
                "rf_status")) {
    if (!col %in% names(summary_df)) summary_df[[col]] <- ""
    summary_df[[col]] <- ifelse(is.na(summary_df[[col]]), "",
                                as.character(summary_df[[col]]))
  }
  summary_df <- init_rf_cols(summary_df)
  cat("Read ESS for", nrow(summary_df), "simulation(s) from", summary_file, "\n")
}

todo <- which(!is.na(summary_df$posterior_ess) & summary_df$posterior_ess > ess_cut)
if (test_mode) {
  todo <- seq_len(nrow(summary_df))   # the test simulation, whatever its ESS
  cat("TEST: ESS cut ignored for simulation", test_sim, "\n")
}

######################################################
# Build CCD (fossirange) and MCC (SRTreeAnnotator) trees where ESS > ess_cut
######################################################
if (do_trees) {
  cat("\nBuilding CCD and MCC trees for", length(todo), "of", nrow(summary_df),
      "simulation(s) with ESS >", ess_cut, "\n")

  for (r in todo) {
    i       <- summary_df$sim[r]
    sim_dir <- file.path(inf_dir, i)
    cat("sim", i, ": ESS =", round(summary_df$posterior_ess[r], 1), "\n")

    summary_df$ccd_status[r] <- run_fossirange(sim_dir)
    cat("  CCD:", summary_df$ccd_status[r], "\n")

    summary_df$mcc_status[r] <- run_mcc(sim_dir)
    cat("  MCC:", summary_df$mcc_status[r], "\n")

    summary_df$translate_status[r] <-
      if (summary_df$mcc_status[r] == "mcc run") run_translate(sim_dir)
      else "skipped (no MCC tree)"
    cat("  translate:", summary_df$translate_status[r], "\n")

    summary_df$status[r] <- paste(summary_df$ccd_status[r],
                                  summary_df$mcc_status[r],
                                  summary_df$translate_status[r], sep = "; ")
  }

  write_summary(summary_df, summary_file)
  cat("\nfossirange run in", sum(summary_df$ccd_status == "fossirange run"), "of",
      length(todo), "eligible folder(s)\n")
  cat("SRTreeAnnotator run in", sum(summary_df$mcc_status == "mcc run"), "of",
      length(todo), "eligible folder(s)\n")
  cat("translate_tree.py run in", sum(summary_df$translate_status == "translate run"),
      "of", length(todo), "eligible folder(s)\n")
  cat("Summary updated in", summary_file, "\n")
}

######################################################
# RF distances (CCD vs MCC, true vs CCD, true vs MCC) where ESS > ess_cut
######################################################
if (do_rf) {
  cat("\nComputing RF distances for", length(todo), "of", nrow(summary_df),
      "simulation(s) with ESS >", ess_cut, "\n")
  if (!is.na(true_trees_file)) cat("True trees from", true_trees_file, "\n")

  for (r in todo) {
    i       <- summary_df$sim[r]
    sim_dir <- file.path(inf_dir, i)
    cat("sim", i, ":\n")

    res <- run_rf(sim_dir, i)
    for (col in rf_cols) summary_df[[col]][r] <- res$values[[col]]
    summary_df$rf_status[r] <- res$status
    cat("  RF:", res$status, "\n")
  }

  write_summary(summary_df, summary_file)
  ran <- startsWith(summary_df$rf_status, "rf run") |
         startsWith(summary_df$rf_status, "rf partial")
  cat("\nRF computed for all pairs in", sum(summary_df$rf_status == "rf run"),
      "and some pairs in", sum(startsWith(summary_df$rf_status, "rf partial")),
      "of", length(todo), "eligible folder(s)\n")
  cat("Summary updated in", summary_file, "\n")

  # one table of every simulation with at least one pair computed
  done <- summary_df[ran, c("sim", "posterior_ess", "rf_status", rf_cols),
                     drop = FALSE]
  done <- done[order(done$sim), , drop = FALSE]
  write_summary(done, rf_summary_file)
  cat("RF results for", nrow(done), "simulation(s) written to", rf_summary_file, "\n")
  for (p in rf_pairs) {
    rf_p <- done[[paste0("rf_", p)]]
    if (sum(!is.na(rf_p)) == 0) next
    cat("  ", p, ": n =", sum(!is.na(rf_p)),
        "; mean RF =", signif(mean(rf_p, na.rm = TRUE), 4),
        "; mean normalised RF =",
        signif(mean(done[[paste0("rf_norm_", p)]], na.rm = TRUE), 4), "\n")
  }
}