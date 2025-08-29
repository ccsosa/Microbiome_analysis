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
  UpSetR
)
pacman::p_load_gh(
  "gauravsk/ranacapa",
  "gmteunisse/fantaxtic",
  "kasperskytte/ampvis2",
  "david-barnett/microViz@0.12.0"
)#,"twbattaglia/btools")


Upsetplot_compo <- function(physeq4,graph_dir,column,levels,name){
  #Subsetting phyloseq in three tables
  for(i in 1:length(levels)){
    level <- levels[[i]]
    message(paste0("Processing...",levels[[i]]))
    if(level=="OTU"){
      physeq_to_common <- physeq4
      taxa_names(physeq_to_common) <- paste0(row.names(tax_table(physeq_to_common)),"-",
                                             tax_table(physeq_to_common)[,7])
    } else {
      physeq_to_common <- physeq4  
      physeq_to_common <- microViz::tax_agg(physeq_to_common,level)
    }
    
    #subseeting
    otu_mat <- data.frame(otu_table(physeq_to_common))
    taxa_names_vec <- data.frame(taxa_names(physeq_to_common))
    samp_mat <- data.frame(sample_data(physeq_to_common))
    
    places <- unique(samp_mat[,column])
    
    # changing otu to presence - absense
    otu_pa <- otu_mat
    otu_pa[otu_pa > 0] <- 1
    
    unit_taxa_list <- lapply(1:length(places),function(j){
      #group
      unit <- places[[j]]
      #getting ids for the group
      unit_ids <- samp_mat$SampleID[which((samp_mat[,column]==unit))]
      #OTU for the unit
      unit_otu <- as.data.frame(otu_pa[,colnames(otu_pa) %in% unit_ids])
      row.names(unit_otu) <- row.names(otu_pa)
      #returning otu_names
      otu_avail <- rowSums(unit_otu)
      x <- names(otu_avail[otu_avail>0])
      return(x)
    })
    
    #assigning names to the list
    names(unit_taxa_list) <- places
    
    #Jaccard UpsetPlot
    venn_patch <- ggVennDiagram(
      unit_taxa_list,
      set_size = 30,
      label_geom = "label",
      force_upset = TRUE
    )
    
    # Edit individual components
    venn_patch[[1]] <- venn_patch[[1]] +
      theme(
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_blank(), 
        axis.line.y = element_blank()
      )
    
    # Just adjust existing label sizes, do not add geom_text
    venn_patch[[2]] <- venn_patch[[2]] +
      theme(
        text = element_text(size = 28),        # General text size
        axis.text.x = element_text(size = 14), # X-axis
        axis.title.x = element_text(size = 16)
      ) 
    venn_patch[[2]]$layers[[2]]$aes_params$size <- 8
    
    venn_patch[[3]] <- venn_patch[[3]] +
      theme(
        axis.text.y = element_text(size = 20),
        axis.title.y = element_text(size = 22),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 22)
      )
    
    # Convert to ggplot-compatible object for export
    venn_plot_combined <- as.ggplot(venn_patch)
    venn_plot_combined <- venn_plot_combined +
      ggtitle(paste0(level))+
      theme(plot.title = element_text(size = 24, face = "bold"))
    
    # Save to high-resolution PDF
    ggsave(
      filename = paste0(graph_dir,"/",name,"", "_Upset_", level, ".pdf"),
      plot = venn_plot_combined,
      width = 15,
      height = 12,
      units = "in",
      dpi = 600,
      device = cairo_pdf
    )
  }
  return("DONE!")
}
