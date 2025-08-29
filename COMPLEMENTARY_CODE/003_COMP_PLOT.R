library(phyloseq)
library(microViz)
library(dplyr)

comp_plot_function <- function(physeq,species_option,name){
  
  theme_custom <- theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 14),
    strip.text = element_text(size = 18, face = "bold"),
    # Facet labels (AM, ECM)
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold")
  )
  
  p1 <- phyloseq::phyloseq(otu_table(physeq),
                           tax_table(physeq),
                           sample_data(physeq))
  
  p1 <- microViz::tax_fix(p1)
 # tax_fix p1 <- tax_fix(p1, unknowns = c("?"))
  
  #Order barplot
  bp1 <- microViz::comp_barplot(
    tax_fix(p1),
    tax_level = "Order",
    bar_outline_colour = NA,
    sample_order = "bray",
    facet_by = "file",
    bar_width = 0.7,
    taxon_renamer = toupper
  ) +
    coord_flip() +
    # facet_wrap("Fungi_system", nrow = 1, scales = "free") +
    labs(x = NULL, y = NULL) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
    theme_custom
  
  #family barplot
  
  bp2 <- microViz::comp_barplot(
    tax_fix(p1),
    tax_level = "Family",
    n_taxa = 10,
    bar_outline_colour = NA,
    sample_order = "bray",
    facet_by = "file",
    bar_width = 0.7,
    taxon_renamer = toupper
  ) +
    coord_flip() +
    # facet_wrap("Fungi_system", nrow = 1, scales = "free") +
    labs(x = NULL, y = NULL) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
    theme_custom
  
  
  #genus barplot
  
  bp3 <- microViz::comp_barplot(
    tax_fix(p1),
    tax_level = "Genus",
    n_taxa = 10,
    bar_outline_colour = NA,
    facet_by = "file",
    sample_order = "bray",
    bar_width = 0.7,
    taxon_renamer = toupper
  ) +
    coord_flip() +
    # facet_wrap("Fungi_system", nrow = 1, scales = "free") +
    labs(x = NULL, y = NULL) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
    theme_custom
  
  #species barplot
  
  if(isTRUE(species_option)){
    bp4 <- microViz::comp_barplot(
      tax_fix(p1),
      tax_level = "Species",
      n_taxa = 10,
      bar_outline_colour = NA,
      facet_by = "file",
      sample_order = "bray",
      bar_width = 0.7,
      taxon_renamer = toupper
    ) +
      coord_flip() +
      # facet_wrap("Fungi_system", nrow = 1, scales = "free") +
      labs(x = NULL, y = NULL) +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
      theme_custom  
    
    bp_total <- gridExtra::grid.arrange(bp1,
                                        bp2, bp3,bp4)
    
    ggsave(
      paste0(graph_dir, "/", name,"_003_BARPLOT_top10_Species.pdf"),
      bp4,
      width = 20,
      height = 25,
      units = "in",
      dpi = 600
    )
    
  } else {
    bp_total <- gridExtra::grid.arrange(bp1,
                                        bp2, bp3)
  }
  
  ggsave(
    paste0(graph_dir, "/",name, "_003_BARPLOT_top10_levels.pdf"),
    bp_total,
    width = 20,
    height = 25,
    units = "in",
    dpi = 600
  )
  
  ggsave(
    paste0(graph_dir, "/", name,"_003_BARPLOT_top10_Genus.pdf"),
    bp3,
    width = 20,
    height = 25,
    units = "in",
    dpi = 600
  )
}

