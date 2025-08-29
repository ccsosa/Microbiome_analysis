#' Plot Regressions of Alpha Diversity Indices Against a Continuous Variable
#'
#' Generates scatter plots with linear regression lines and Pearson correlation coefficients
#' between alpha diversity indices (Observed species, Shannon, Inverse Simpson, Pielou) 
#' and a continuous response variable (e.g., height, pH). Each index is plotted separately, 
#' grouped by a categorical variable, and arranged into a grid. The result is saved as a high-resolution PDF.
#'
#' @param Alpha A `data.frame` containing alpha diversity indices and metadata. Must include columns: 
#' `"Observed"`, `"Shannon"`, `"InvSimpson"`, `"pielou"`, the grouping variable (`col_fill_var`), and the response variable (`y_var`).
#' @param col_fill_var Character string. The name of the column in `Alpha` used to color and group samples (e.g., `"Ecosystem"` or `"Fungi_system"`).
#' @param y_var Character string. The name of the continuous response variable to regress against alpha diversity indices.
#' @param saved_name Character string. Prefix used in the filename for the saved PDF plot.
#'
#' @return No object is returned. A combined plot of all regressions is saved in the `graphics` folder as a `.pdf` file.
#'
#' @details
#' This function uses `ggscatter()` from the `ggpubr` package to generate regression plots for:
#' - Observed species richness
#' - Shannon diversity
#' - Inverse Simpson index
#' - Pielou's evenness
#'
#' It overlays regression lines with confidence intervals and prints Pearson correlation coefficients.
#'
#' @examples
#' \dontrun{
#' data(Alpha)  # Your alpha diversity + metadata table
#' Regression_alpha_function(Alpha, col_fill_var = "Fungi_system", y_var = "Height", saved_name = "experiment1")
#' }
#'
#' @importFrom ggpubr ggscatter stat_cor
#' @importFrom gridExtra grid.arrange
#' @import ggpubr
#' @import gridExtra
#' @export

# Load required libraries
pacman::p_load(phyloseq, ggpubr, gridExtra)


Regression_alpha_function <- function(Alpha, col_fill_var, y_var, saved_name) {
  
  # Dictionary of alpha diversity metrics and their labels for x-axis
  index_names <- c(
    Observed = "Observed species",
    Shannon = "Shannon index",
    InvSimpson = "Effective Number of species (1/D)",
    pielou = "Pielou evenness"
  )
  
  # Safety check: ensure all necessary columns are in the input data
  required_cols <- c(names(index_names), col_fill_var, y_var)
  missing <- setdiff(required_cols, names(Alpha))
  if (length(missing) > 0) stop("Missing columns in 'Alpha': ", paste(missing, collapse = ", "))
  
  # Use lapply to iterate over each index and generate the corresponding ggscatter plot
  plot_list <- lapply(names(index_names), function(index) {
    
    # Create scatter plot with regression line and confidence interval
    ggscatter(
      data = Alpha,
      x = index,
      y = y_var,
      color = col_fill_var,   # Coloring by group
      fill = col_fill_var,    # Filling points by group
      palette = "jco",        # Color palette
      shape = 19,             # Point shape
      size = 3,               # Point size
      add = "reg.line",       # Add linear regression line
      conf.int = TRUE,        # Add 95% confidence interval
      cor.method = "pearson", # Method for correlation (used for line fitting)
      xlab = index_names[index], # Dynamic x-axis label
      ylab = y_var,           # y-axis label based on input
      ggtheme = theme_pubr(base_size = 20) +  # Set base font size for theme
        theme(
          axis.text = element_text(size = 18),         # Axis tick label size
          axis.title = element_text(size = 22, face = "bold"), # Axis title style
          legend.text = element_text(size = 16),       # Legend text size
          legend.title = element_text(size = 18)       # Legend title size
        )
    ) +
      # Add Pearson correlation coefficient with increased size
      stat_cor(aes(color = Alpha[[col_fill_var]]), size = 6)
  })
  
  # Combine all individual plots into a grid layout (2x2 by default)
  combined_plot <- gridExtra::grid.arrange(grobs = plot_list)
  
  # Display message to console
  message(paste("Alpha diversity regressions with", y_var, "saved to graphics folder"))
  
  # Save the combined plot as a high-resolution PDF
  ggsave(
    filename = paste0(graph_dir, "/", "003_regression_", y_var, "_", saved_name, ".pdf"),
    plot = combined_plot,
    width = 18, height = 18, units = "in", dpi = 600
  )
}
