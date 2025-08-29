library(ggVennDiagram)
library(ggplot2)
library(ggpubr)
library(phyloseq)
library(microViz)



Jaccard_upsetplot <- function(levels, physeq,markers_name,name){
  message("Starting upsetplot/Jaccard Script")
   i <- 2
  x_jacc <- lapply(1:length(levels), function(i) {

    level <- levels[[i]]
    message(paste0("Processing: ",level," for ", name))
    
    #Filtering Illumina
    ILL_filtered <- phyloseq::subset_samples(physeq,file=="ITS_ILLUMINA")
    ILL_filtered <-   phyloseq::prune_taxa(taxa_sums(ILL_filtered) > 0, ILL_filtered)
    
    #Filtering PACBIO
    PAC_filtered <- phyloseq::subset_samples(physeq,file=="ITS_PACBIO")
    PAC_filtered <-   phyloseq::prune_taxa(taxa_sums(PAC_filtered) > 0, PAC_filtered)
    
    #Joining in a list
    physeq_list <- list(ILL_filtered,PAC_filtered)
    names(physeq_list) <- markers_name
    
    phy_aggregate <- lapply(1:length(physeq_list),function(j){
      #Calling object, filter, and pruning
      phy1 <- physeq_list[[j]]
      phy1 <- tax_fix(phy1,unknowns = c("?"))
      phy1 <- microViz::tax_agg(ps = phy1, level)
      phy1 <- phyloseq::prune_taxa(taxa_sums(phy1) > 0, phy1)
      return(phy1)
    })
    
    #get taxonomic names
    phy_names <- lapply(1:length(phy_aggregate), function(j) {
      # j <- 1
      tax1 <- tax_table(phy_aggregate[[j]])
      otu1 <- as.character(unlist(row.names(tax1)))
      return(otu1)
    })
    
    names(phy_names) <- markers_name
    message("Plotting taxonomic similarities")
    
    #Jaccard UpsetPlot
    venn_patch <- ggVennDiagram(
      phy_names,
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
      ggtitle(paste0("Overlap of ", level))+
      theme(plot.title = element_text(size = 24, face = "bold"))
    
    # Save to high-resolution PDF
    ggsave(
      filename = paste0(graph_dir,"/",name,"_VENN", "_", level, ".pdf"),
      plot = venn_plot_combined,
      width = 15,
      height = 10,
      units = "in",
      dpi = 600,
      device = cairo_pdf
    )
    ##############################################################################
    ##############################################################################
    ##############################################################################
    #Jaccard
    message("Calculating Jaccard index")
    
    jacc_mat <- matrix(ncol=length(phy_names),nrow=length(phy_names))
    jacc_mat <- as.data.frame(jacc_mat)
    colnames(jacc_mat) <- markers_name
    row.names(jacc_mat) <- markers_name
    for(k in 1:length(phy_names)){
      for(l in 2:length(phy_names)){
        if(l>k){
          jacc_mat[k,l] <- 
            jaccard(phy_names[[k]],phy_names[[l]])
        }
      }
    }
    diag(jacc_mat) <- 1
    write.csv(jacc_mat,paste0(csv_dir, "/",name,"_JACCARD", "_", level, ".csv"))
    
    #Getting phy
    phy_names_unique_filtered <- lapply(1:length(phy_names), function(j) {
      tax1 <- tax_table(phy_aggregate[[j]])
      otu1 <- unlist(row.names(tax1))
      otu1 <- as.character(otu1)
      # otu1 <- otu1[!otu1 %in% common]
      if(length(otu1)>0){
        otu1 <- data.frame(taxon=otu1,dataset=markers_name[[j]])        
      } else {
        otu1 <- data.frame(taxon=NA,dataset=markers_name[[j]])
      }
      return(otu1)
    })
    
    phy_names_unique_filtered <- do.call(rbind,phy_names_unique_filtered)
    write.csv(phy_names_unique_filtered,paste0(csv_dir, "/",name,"_UNIQUE_TAXON", "_","_", level, ".csv"),row.names=F)
  })
}
