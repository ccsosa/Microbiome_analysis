################################################################################
#Load pacman to load complementary libraries
library(pacman)

#Load libraries or install if it is the case from CRAN and github respetively
pacman::p_load(
  phyloseq,
  microbiome,
  ggplot2,
  ggpubr,
  xlsx,
  readr,
  dplyr,
  compositions,
  ggrepel,
  iNEXT,
  gridExtra,
  ggsci,
  FactoMineR,
  factoextra,
  hrbrthemes,
  gcookbook,
  Hmisc,
  corrplot
)
pacman::p_load_gh(
  "gauravsk/ranacapa",
  "gmteunisse/fantaxtic",
  "kasperskytte/ampvis2",
  "david-barnett/microViz@0.12.0"
)#,"twbattaglia/btools")

corrplot_function <- function(x,graph_dir,csv_dir){
  # Use the rcorr() function from the Hmisc package.
  # It returns both correlation coefficients (r) and p-values.
  correlation_results <- Hmisc::rcorr(as.matrix(x), type = "spearman")
  
  # Extract the correlation matrix (r) and the p-value matrix (P)
  r_values <- correlation_results$r
  p_values <- correlation_results$P
  #getting cors < 0.05 to save in a CSV
  r_values_fix <- r_values
  r_values_fix[which(p_values>0.05)] <- NA
  
  # Create the plot
  x_cor <- corrplot::corrplot(
    r_values,
    p.mat = p_values,      # Provide the p-value matrix
    method = "circle",     # Use circles to represent correlations (can also be "number", "color", "pie")
    type = "upper",        # Show only the upper triangle of the matrix
    sig.level = 0.05,      # Set the significance level (e.g., 0.05)
    insig = "blank",       # Leave non-significant correlations blank
    order = "hclust",      # Reorder variables based on hierarchical clustering
    tl.col = "black",      # Color of text labels
    tl.srt = 90,           # Rotation of text labels
    diag = FALSE,          # Do not show the diagonal
    title = "Correlation between Alpha Diversity and Soil Variables",
    mar = c(0,0,1,0)       # Adjust plot margins
  )
  
  png(
    filename = file.path(graph_dir, "/COR_PCA.png"),
    type = "cairo",
    width = 10,
    height = 10,
    units = "in",
    res = 600
  )
  corrplot::corrplot(
    r_values,
    p.mat = p_values,      # Provide the p-value matrix
    method = "circle",     # Use circles to represent correlations (can also be "number", "color", "pie")
    type = "upper",        # Show only the upper triangle of the matrix
    sig.level = 0.05,      # Set the significance level (e.g., 0.05)
    insig = "blank",       # Leave non-significant correlations blank
    order = "hclust",      # Reorder variables based on hierarchical clustering
    tl.col = "black",      # Color of text labels
    tl.srt = 90,           # Rotation of text labels
    diag = FALSE,          # Do not show the diagonal
    title = "Correlation between Alpha Diversity and Soil Variables",
    mar = c(0,0,1,0)  
  )
  dev.off()
  #saving CSV
  write.csv(data.frame(r_values_fix),paste0(csv_dir, "/", "000_COR_PCA.csv"))
  return("DONE!")
}
