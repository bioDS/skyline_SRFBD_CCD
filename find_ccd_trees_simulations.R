######################################################
# Posterior ESS check + fossirange CCD + MCC summary for sRanges runs
#
# Usage:
#   Rscript find_ccd_trees_simulations.R PARENT [MODE]
#
#   PARENT: four_interval_simulation, one_interval_simulation,
#           or a full path to either
#   MODE:   ess        - compute posterior ESS only, write ess_summary.csv
#           fossirange - read ESS from PARENT/inf/ess_summary.csv and build
#                        the CCD and MCC trees where ESS > ess_cut
#           both       - compute ESS and build the trees in one pass (default)
#
# Examples:
#   Rscript find_ccd_trees_simulations.R four_interval_simulation ess
#   Rscript find_ccd_trees_simulations.R four_interval_simulation fossirange
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
# The summary table PARENT/inf/ess_summary.csv
# (sim, posterior_ess, status, ccd_status, mcc_status, translate_status)
# is written by every mode.
######################################################
library(coda)

base_dir  <- "/home/ket581/skyline/SUMMER/skyline-SRFBD-validation/"
allowed   <- c("four_interval_simulation", "one_interval_simulation")
modes     <- c("ess", "fossirange", "both")
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

######################################################
# Arguments and paths
######################################################
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1)
  stop("Usage: Rscript find_ccd_trees_simulations.R PARENT [ess|fossirange|both]")

mode <- if (length(args) >= 2) tolower(args[2]) else "both"
if (!mode %in% modes)
  stop("MODE must be one of: ", paste(modes, collapse = ", "))

parent <- sub("/+$", "", args[1])
if (!basename(parent) %in% allowed)
  stop("Parent folder must be one of: ", paste(allowed, collapse = ", "))
# accept either a path that exists as given or a name relative to base_dir
if (!dir.exists(parent)) parent <- paste0(base_dir, basename(parent))
parent  <- normalizePath(parent, mustWork = TRUE)   # absolute, so setwd can't break it
inf_dir <- file.path(parent, "inf")
if (!dir.exists(inf_dir)) stop("No inference folder at ", inf_dir)
summary_file <- file.path(inf_dir, "ess_summary.csv")

if (mode %in% c("fossirange", "both")) {
  if (Sys.which("fossirange") == "")
    stop("fossirange not found on PATH")
  if (!file.exists(mcc_bin) || file.access(mcc_bin, 1) != 0)
    stop("applauncher not found or not executable at ", mcc_bin)
  if (Sys.which(python_bin) == "")
    stop(python_bin, " not found on PATH")
  if (!file.exists(translate_script))
    stop("translate_tree.py not found at ", translate_script)
  translate_script <- normalizePath(translate_script)   # absolute, used after setwd
}

cat(basename(parent), "- mode:", mode, "\n")

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

######################################################
# Get ESS: compute it (ess / both) or read it back (fossirange)
######################################################
if (mode %in% c("ess", "both")) {
  sims <- list.dirs(inf_dir, full.names = FALSE, recursive = FALSE)
  sims <- sims[grepl("^[0-9]+$", sims)]
  sims <- sims[order(as.integer(sims))]
  cat(length(sims), "simulation folder(s)\n")

  summary_df <- do.call(rbind, lapply(sims, function(i) {
    res <- compute_ess(file.path(inf_dir, i))
    cat("sim", i, ": ESS =", round(res$ess, 1), "-", res$status, "\n")
    data.frame(sim = as.integer(i), posterior_ess = res$ess, status = res$status)
  }))
  summary_df$ccd_status <- ""
  summary_df$mcc_status <- ""
  summary_df$translate_status <- ""
  write.csv(summary_df, summary_file, row.names = FALSE)
  cat("ESS summary written to", summary_file, "\n")
} else {
  if (!file.exists(summary_file))
    stop("No ", summary_file, " - run with mode 'ess' first")
  summary_df <- read.csv(summary_file, stringsAsFactors = FALSE)
  if (!all(c("sim", "posterior_ess") %in% names(summary_df)))
    stop(summary_file, " needs columns sim and posterior_ess")
  for (col in c("status", "ccd_status", "mcc_status", "translate_status")) {
    if (!col %in% names(summary_df)) summary_df[[col]] <- ""
    summary_df[[col]] <- ifelse(is.na(summary_df[[col]]), "",
                                as.character(summary_df[[col]]))
  }
  cat("Read ESS for", nrow(summary_df), "simulation(s) from", summary_file, "\n")
}

######################################################
# Build CCD (fossirange) and MCC (SRTreeAnnotator) trees where ESS > ess_cut
######################################################
if (mode %in% c("fossirange", "both")) {
  todo <- which(!is.na(summary_df$posterior_ess) & summary_df$posterior_ess > ess_cut)
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

  write.csv(summary_df, summary_file, row.names = FALSE)
  cat("\nfossirange run in", sum(summary_df$ccd_status == "fossirange run"), "of",
      length(todo), "eligible folder(s)\n")
  cat("SRTreeAnnotator run in", sum(summary_df$mcc_status == "mcc run"), "of",
      length(todo), "eligible folder(s)\n")
  cat("translate_tree.py run in", sum(summary_df$translate_status == "translate run"),
      "of", length(todo), "eligible folder(s)\n")
  cat("Summary updated in", summary_file, "\n")
}