#!/usr/bin/env Rscript
# Script to help refactor library() calls to requireNamespace()
#
# This script provides guidance for manual refactoring

cat("Library Call Refactoring Guide\n")
cat("================================\n\n")

cat("Pattern to follow:\n\n")

cat("BEFORE:\n")
cat("  library(ggplot2)\n")
cat("  library(dplyr)\n")
cat("  plot <- ggplot(data) + geom_point()\n\n")

cat("AFTER:\n")
cat("  if (!requireNamespace('ggplot2', quietly = TRUE)) {\n")
cat("    stop('Package ggplot2 required. Install with: install.packages(\"ggplot2\")')\n")
cat("  }\n")
cat("  if (!requireNamespace('dplyr', quietly = TRUE)) {\n")
cat("    stop('Package dplyr required. Install with: install.packages(\"dplyr\")')\n")
cat("  }\n")
cat("  plot <- ggplot2::ggplot(data) + ggplot2::geom_point()\n\n")

cat("Note: For packages already in Imports in DESCRIPTION, you can use :: directly\n")
cat("without the check, but requireNamespace() is best practice for optional deps.\n\n")

cat("Files to refactor (in priority order):\n")
cat("1. R/sctype_visualize.R (4 library calls)\n")
cat("2. R/sctype_uncertainty.R (11 library calls)\n")
cat("3. R/sctype_pathway_enrichment.R (9 library calls)\n")
cat("4. R/sctype_visualize_sce.R (8 library calls)\n")
cat("5. R/sctype_uncertainty_sce.R (10 library calls)\n")
cat("6. R/sctype_hierarchical_sce.R (2 library calls)\n")
cat("7. R/sctype_wrapper_sce.R (1 library call)\n")
