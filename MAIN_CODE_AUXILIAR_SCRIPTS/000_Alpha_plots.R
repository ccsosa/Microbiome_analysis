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
  gcookbook
)
pacman::p_load_gh(
  "gauravsk/ranacapa",
  "gmteunisse/fantaxtic",
  "kasperskytte/ampvis2",
  "david-barnett/microViz@0.12.0"
)#,"twbattaglia/btools")


alpha_plots <- function(physeq4,ps_rare,csv_dir){
  

################################################################################
#Alpha richness analysis
message("Getting CSVs of Alpha diversity")
# #obtaining richness for non rarefied Phyloseq object
Alpha <- phyloseq::estimate_richness(physeq4)
Alpha$pielou <- microbiome::evenness(physeq4, index = "pielou")[, 1]
calc <- iNEXT::ChaoRichness(x = otu_table(physeq4), datatype = "abundance", conf = 0.95)
#Adding sample data
Alpha <- cbind(sample_data(physeq4), Alpha)
message("Saving Alpha diversity indexes CSV files")
write.csv(
  Alpha,
  paste0(csv_dir, "/", "002_Alpha.csv"),
  na = "",
  row.names = F,
  quote = T
)
################################################################################
#Alpha richness analysis (Rarefied)
# #obtaining richness for non rarefied Phyloseq object
Alpha_rare <- phyloseq::estimate_richness(ps_rare3)
Alpha_rare$pielou <- microbiome::evenness(ps_rare3, index = "pielou")[, 1]
calc <- iNEXT::ChaoRichness(x = otu_table(ps_rare3), datatype = "abundance", conf = 0.95)
#Adding sample data
Alpha_rare <- cbind(sample_data(ps_rare3), Alpha_rare)
message("Saving Alpha diversity indexes CSV files")
write.csv(
  Alpha_rare,
  paste0(csv_dir, "/", "002_Alpha_rare.csv"),
  na = "",
  row.names = F,
  quote = T
)
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
message("Plot Alpha diversity")

#PLOTS (NO RAREFIED)
#M_TYPE
# message("Plotting boxplots of Alpha diversity for rarefied and filtered Phyloseq object")
#Plot boxplots with Alpha indexes and including Wilcoxon tests (Habit) Michorryzal type
p_r_p4 <-
  phyloseq::plot_richness(
    physeq4,
    x = "Habit",
    color = "Habit",
    nrow = 2,
    measures = c("Observed", "Shannon","InvSimpson","Chao1"),
    title = ""
  ) +
  ylab("")+
  xlab("")+
  geom_boxplot(aes(fill = Habit), alpha = 0.7,outlier.shape = NA,coef=0) +
  stat_compare_means(method = "wilcox.test",
                     size = 6,
                     aes(label = paste0("p = ", after_stat(p.format)))) +
  ggsci::scale_color_jco() +   # jco palette for lines/points
  ggsci::scale_fill_jco() +    # jco palette for boxplots
  theme(
    legend.position = "None",
    strip.text = element_text(size = 16),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 22),
    plot.title = element_text(size = 24, face = "bold", hjust = 0.5)
  )


message("Alpha indicators...Saved in graphics subfolder")
ggsave(
  paste0(graph_dir, "/", "002_Alpha_boxplots_M_TYPE.pdf"),
  p_r_p4,
  width = 20,
  height = 25,
  units = "in",
  dpi = 600
)

################################################################################
jco_palette <- ggsci::pal_jco()(7)
jco_more <- colorRampPalette(jco_palette)(length(unique(samp_data$Site_name.x))) # generar 20 colores
################################################################################
#SITES
p_r_p4 <-
  phyloseq::plot_richness(
    physeq4,
    x = "Site_name.x",
    color = "Site_name.x",
    nrow = 2,
    measures = c("Observed", "Shannon","InvSimpson","Chao1"),
    title = "",sortby = "Observed"
  ) +
  ylab("")+
  xlab("")+
  geom_boxplot(aes(fill = Site_name.x), alpha = 0.7,outlier.shape = NA,coef=0) +
  # stat_compare_means(method = "wilcox.test",
  #                    size = 6,
  #                    comparisons = my_comparisons,
  #                    
  #                    aes(label = paste0("p = ", after_stat(p.format)))) +
  # ggsci::scale_color_jco() +   # jco palette for lines/points
  # ggsci::scale_fill_jco() +    # jco palette for boxplots
  scale_color_manual(values = jco_more) +
  scale_fill_manual(values = jco_more) +
  
  theme(
    legend.position = "None",
    strip.text = element_text(size = 16),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 22),
    plot.title = element_text(size = 24, face = "bold", hjust = 0.5)
  )


message("Alpha indicators...Saved in graphics subfolder")
ggsave(
  paste0(graph_dir, "/", "002_Alpha_boxplots_Sites.pdf"),
  p_r_p4,
  width = 20,
  height = 25,
  units = "in",
  dpi = 600
)
################################################################################
################################################################################
################################################################################
################################################################################
message("Plot Alpha diversity (RAREFIED)")

# message("Plotting boxplots of Alpha diversity for rarefied and filtered Phyloseq object")
#Plot boxplots with Alpha indexes and including Wilcoxon tests (Habit) Michorryzal type
p_r_p4 <-
  plot_richness(
    ps_rare3,
    x = "Habit",
    color = "Habit",
    nrow = 2,
    measures = c("Observed", "Shannon","InvSimpson","Chao1"),
    title = ""
  ) +
  ylab("")+
  xlab("")+
  geom_boxplot(aes(fill = Habit), alpha = 0.7,outlier.shape = NA,coef=0) +
  stat_compare_means(method = "wilcox.test",
                     size = 6,
                     aes(label = paste0("p = ", after_stat(p.format)))) +
  ggsci::scale_color_jco() +   # jco palette for lines/points
  ggsci::scale_fill_jco() +    # jco palette for boxplots
  theme(
    legend.position = "None",
    strip.text = element_text(size = 16),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 22),
    plot.title = element_text(size = 24, face = "bold", hjust = 0.5)
  )


message("Alpha indicators...Saved in graphics subfolder")
ggsave(
  paste0(graph_dir, "/", "002_Alpha_boxplots_M_TYPE_RARE.pdf"),
  p_r_p4,
  width = 20,
  height = 25,
  units = "in",
  dpi = 600
)

################################################################################
jco_palette <- ggsci::pal_jco()(7)
jco_more <- colorRampPalette(jco_palette)(length(unique(samp_data$Site_name.x))) # generar 20 colores
################################################################################
# my_comparisons <- combn(
#   unique(sample_data(physeq4)$Site_name.x), 
#   2, 
#   simplify = FALSE
# )
p_r_p4 <-
  plot_richness(
    ps_rare3,
    x = "Site_name.x",
    color = "Site_name.x",
    nrow = 2,
    measures = c("Observed", "Shannon","InvSimpson","Chao1"),
    title = "",sortby = "Observed"
  ) +
  ylab("")+
  xlab("")+
  geom_boxplot(aes(fill = Site_name.x), alpha = 0.7,outlier.shape = NA,coef=0) +
  # stat_compare_means(method = "wilcox.test",
  #                    size = 6,
  #                    comparisons = my_comparisons,
  #                    
  #                    aes(label = paste0("p = ", after_stat(p.format)))) +
  # ggsci::scale_color_jco() +   # jco palette for lines/points
  # ggsci::scale_fill_jco() +    # jco palette for boxplots
  scale_color_manual(values = jco_more) +
  scale_fill_manual(values = jco_more) +
  
  theme(
    legend.position = "None",
    strip.text = element_text(size = 16),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 22),
    plot.title = element_text(size = 24, face = "bold", hjust = 0.5)
  )


message("Alpha indicators...Saved in graphics subfolder")
ggsave(
  paste0(graph_dir, "/", "002_Alpha_boxplots_Sites_RARE.pdf"),
  p_r_p4,
  width = 20,
  height = 25,
  units = "in",
  dpi = 600
)

message("returning Alpha tables")
x <- list(Alpha = Alpha,
          Alpha_rare = Alpha_rare)
return(x)
}