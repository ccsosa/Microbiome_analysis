library(phyloseq)
library(ggpubr)
library(ggrepel)
library(ggsci)
library(microViz)
library(microbiome)


BETA_func <- function(physeq_clr2,name){
  
  # Perform PCoA ordination using Bray-Curtis distance
  psd5.mds.euc2 <-
    ordinate(physeq_clr2, method = "PCoA", distance = "bray")
  evals <- psd5.mds.euc2$values$Eigenvalues
  
  
  # Extract files
  ord_df <- as.data.frame(plot_ordination(physeq_clr2, psd5.mds.euc2, justDF = TRUE))
  
  #Extract sequencing technology
  ord_df$file <- sample_data(physeq_clr2)$file
  #Extracting labels
  ord_df$SampleID <- sample_data(physeq_clr2)$SampleID
  
  # Obtaining convex hull
  hulls <- ord_df %>%
    group_by(file) %>%
    slice(chull(Axis.1, Axis.2))  # Convex hull para cada grupo
  
  
  # Create the ordination plot and apply the 'jco' color palette
  pord2 <-
    plot_ordination(physeq_clr2, psd5.mds.euc2, color = "file", shape = "file") +
    geom_polygon(data = hulls, aes(x = Axis.1, y = Axis.2, fill = file), alpha = 0.2, color = NA) +
    
    # labs(col = "Fungi System") +
    geom_point(size = 8) +  # Increase point size
    
    ggsci::scale_color_jco() +         # Color de puntos
    ggsci::scale_fill_jco() +          # Color de polígonos (igual que los puntos)
    # coord_fixed(sqrt(evals[2] / evals[1])) +
    #theme_pubr(base_size = 20) +  # Increase all base text size
    theme(
      axis.title = element_text(size = 22),
      axis.text = element_text(size = 20),
      strip.text = element_text(size = 20),
      legend.title = element_text(size = 0),
      legend.text = element_text(size = 18)
    ) + 
    
    # ggrepel to add labels
    geom_text_repel(
      data = ord_df,
      aes(x = Axis.1, y = Axis.2, label = SampleID, color = file),
      size = 10,max.overlaps = 100000,
      show.legend = FALSE
    )
  
  # Show the plot
  message("Beta diversity...Saved in graphics subfolder")
  
  ggsave(
    paste0(graph_dir, "/", name,"_004_Beta_PCoA_EUC.pdf"),
    pord2,
    width = 20,
    height = 20,
    units = "in",
    dpi = 600
  )
  
  
  #perform PERMANOVA
  #USING NO RAREFIED OBJECT + CLR
  #Additionally, our betadisper results are not significant,
  #meaning we cannot reject the null hypothesis that our groups have the same dispersions
  dist = phyloseq::distance(physeq_clr2, method="bray")
  metadata <- data.frame(sample_data(physeq_clr2))
  #using file to perform PERMANOVA
  test.adonis <- adonis2(dist ~ file, data = metadata,permutations = 10000,parallel = T)
  test.adonis <- as.data.frame(test.adonis)
  test.adonis$interpretation <- NA
  if(test.adonis$`Pr(>F)`[1]<0.05){
    test.adonis$interpretation[1] <- "Differences of groups"
  }
  
  write.csv(test.adonis,paste0(csv_dir,"/",name,"_004_PERMANOVA_CLR.csv"))
  
  aov_beta <- as.data.frame(anova(betadisper(dist,ord_df$file)))
  aov_beta$interpretation <- NA
  if(aov_beta$`Pr(>F)`[1]>0.05){
    aov_beta$interpretation[1] <- "No reject null hypothesis (same dispersions)"
  }
  write.csv(aov_beta,paste0(csv_dir,"/",name,"_004_BETADISPERSER_CLR.csv"))
}