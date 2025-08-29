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
  grDevices
)
pacman::p_load_gh(
  "gauravsk/ranacapa",
  "gmteunisse/fantaxtic",
  "kasperskytte/ampvis2",
  "david-barnett/microViz@0.12.0"
)#,"twbattaglia/btools")


compo_plot <- function(physeq,level,ntop,column){
  
  #ABUNDANCE BARPLOTS
  # https://microbiome.github.io/tutorials/Composition.html
  p <- microViz::tax_agg(physeq,rank = level)
  # Obtener top 20 taxones por abundancia total
  top_taxa <- names(sort(taxa_sums(p), decreasing = TRUE))[1:ntop]
  
  # Filtrar phyloseq solo con esos taxones
  physeq_top <- prune_taxa(top_taxa, p)
  physeq_top <- phyloseq::prune_taxa(taxa_sums(physeq_top) > 0, physeq_top)
  physeq_top <- microbiome::transform(physeq_top, "compositional")
  # Averaged by group
  p2 <- plot_composition(physeq_top,
                         average_by = column 
                         # transform = "compositional"
  ) +
    #scale_fill_brewer("Genera", palette = "Paired") +
    theme_ipsum(grid="Y") +
    xlab("")+
    theme(axis.text.x = element_text(angle=90, hjust=1),
          legend.text = element_text(face = "italic"),
          legend.title = element_text(size = 0))+
    scale_y_percent() +
    ggsci::scale_color_jco() +   # jco palette for lines/points
    ggsci::scale_fill_jco() +    # jco palette for boxplots
  theme(text = element_text(family = "Arial"))
  
  
  return(p2)
}