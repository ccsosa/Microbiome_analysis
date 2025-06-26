# Load required libraries
pacman::p_load(phyloseq, ggpubr, gridExtra, parallel, microbiome, microViz)

#' Identify Differentially Abundant Taxa at a Given Taxonomic Level
#'
#' This function identifies differentially abundant taxa using either
#' the Wilcoxon rank-sum test (2 groups) or the Kruskal-Wallis test (>2 groups).
#' It works in parallel using multiple cores.
#'
#' @param physeq_obj A phyloseq object.
#' @param level Taxonomic level to aggregate to (e.g., "Phylum", "Genus", "OTU").
#' @param col_fill_var The column name in sample_data used as grouping variable.
#' @param ncores Number of CPU cores to use for parallel processing.
#'
#' @return A data.frame with taxonomy, raw p-values, FDR-corrected p-values (q-values) for significant taxa (FDR < 0.05). Returns NULL if no significant taxa are found.

Differential_abundant_taxa <- function(physeq_obj, level, col_fill_var, ncores) {
  
  # Validate level argument
  if (is.null(level) | !level %in% c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "OTU")) {
    stop("Please use a valid taxonomic level.")
  }
  
  #Fixing ? issue
  physeq_obj <- microViz::tax_fix(physeq_obj, unknowns = c("?"))
  
  
  # Aggregate data to taxonomic level if not OTU
  if (level == "OTU") {
    pG <- physeq_obj
  } else {
    message(paste0("Aggregating taxa to level: ", level))
    pG <- microViz::tax_agg(ps = physeq_obj, level)
    pG <- microbiome::transform(pG, "compositional")  # Convert to relative abundance
  }
  
  # Extract grouping variable from sample_data
  groups <- sample_data(pG)[, col_fill_var]
  
  # Transpose OTU table to samples x taxa
  message("Obtaining OTU table from the Phyloseq object")
  otu_tab <- as.data.frame(t(otu_table(pG)))
  for (i in seq_len(ncol(otu_tab))) {
    otu_tab[, i] <- as.numeric(otu_tab[, i])
  }
  rm(i)
  
  # Start parallel cluster
  cl <- parallel::makeCluster(ncores)
  parallel::clusterExport(cl, varlist = c("otu_tab", "groups"), envir = environment())
  
  # Choose test based on number of groups
  if (length(unique(groups[[1]])) > 2) {
    message(paste0("Using Kruskal-Wallis test...","Processing: ",ncol(otu_tab)," operations"))
    graph_db1 <- parallel::parLapplyLB(cl,
                                       X = seq_len(ncol(otu_tab)),
                                       fun = function(i) {
                                         x <- data.frame(
                                           taxon = as.numeric(otu_tab[[i]]),
                                           group = as.factor(groups[[1]])
                                         )
                                         tryCatch(
                                           kruskal.test(taxon ~ group, data = x)$p.value,
                                           error = function(e) NA
                                         )
                                       })
  } else {
    message(paste0("Using Wilcoxon rank-sum test...","Processing: ",ncol(otu_tab)," operations"))
    graph_db1 <- parallel::parLapplyLB(cl,
                                       X = seq_len(ncol(otu_tab)),
                                       fun = function(i) {
                                         x <- data.frame(
                                           taxon = as.numeric(otu_tab[[i]]),
                                           group = as.factor(groups[[1]])
                                         )
                                         tryCatch(
                                           wilcox.test(taxon ~ group, data = x, exact = FALSE)$p.value,
                                           error = function(e) NA
                                         )
                                       })
  }
  
  # Stop the parallel cluster
  parallel::stopCluster(cl)
  
  message("Adding results (p-values and FDR values)")
  
  # Combine with taxonomy
  results_df <- as.data.frame(tax_table(pG))
  results_df$p <- unlist(graph_db1)
  results_df$fdr <- p.adjust(results_df$p, method = "fdr")
  
  # Filter for significant taxa
  results_df <- results_df[which(results_df$fdr < 0.05), ]
  
  # Return results or NULL
  message("Saving results if they are available")
  if (nrow(results_df) > 0) {
    return(results_df)
  } else {
    message(paste0("No differentially abundant taxa for level: ", level))
    return(NULL)
  }
}
