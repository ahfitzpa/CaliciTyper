#!/usr/bin/env bash
#SBATCH --job-name=bmj_pal
#SBATCH --partition=shared
#SBATCH --cpus-per-task=2
#SBATCH --time=30
#SBATCH --output=/home/people/fitzpatria/scratch/calicityper_runs/logs/bmj_pal_%j.out
set -uo pipefail
export R_LIBS_USER=/home/people/fitzpatria/scratch/calicityper_runs/manuscript/rlib
mkdir -p "$R_LIBS_USER"
module load r/4.4.1-gcc-11.5.0-plddwx2 2>/dev/null || module load R 2>/dev/null
echo "### $(hostname) $(date); R: $(which Rscript)"
Rscript -e '
  lib <- Sys.getenv("R_LIBS_USER"); .libPaths(c(lib, .libPaths()))
  if (!requireNamespace("ggsci", quietly=TRUE))
    install.packages("ggsci", lib=lib, repos="https://cloud.r-project.org", quiet=TRUE)
  library(ggsci)
  cat("ggsci version:", as.character(packageVersion("ggsci")), "\n")
  db <- getFromNamespace("ggsci_db", "ggsci")
  for (nm in c("bmj","jama","nejm","lancet")) {
    if (!is.null(db[[nm]])) {
      for (v in names(db[[nm]])) cat(sprintf("%s/%s: %s\n", nm, v, paste(db[[nm]][[v]], collapse=" ")))
    } else cat(nm, ": not in this ggsci version\n")
  }
'
echo "### DONE"
