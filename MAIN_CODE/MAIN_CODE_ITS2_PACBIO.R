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
  ggsci
)
pacman::p_load_gh(
  "gauravsk/ranacapa",
  "gmteunisse/fantaxtic",
  "kasperskytte/ampvis2",
  "david-barnett/microViz@0.12.0"
)#,"twbattaglia/btools")
################################################################################
#Loading scripts
message("Loading complementary scripts")
source("~/SCRIPTS/000_Prepare_dirs.R")
source("~/SCRIPTS/001_Accumulation_Species.R")
source("~/SCRIPTS/002_regression_alpha_vars.R")
source("~/SCRIPTS/003_Differentially_abundant_taxon.R")


#Home folder
dir <- "~"
#how to name the results
marker="ITS_PACBIO"
#Add directories
dirs <- dir_function(dir,marker=marker)
################################################################################
#Defining folder to use in the pipeline
SCRIPTS_dir <- dirs$SCRIPTS_dir
out_dir <-dirs$out_dir
graph_dir <-dirs$graph_dir
csv_dir <-dirs$csv_dir
RDS_dir <-dirs$RDS_dir
################################################################################
#Fixing a seed. Just in case
set.seed(1000)
################################################################################
#Load Phyloseq object
message("Loading Phyloseq object")
physeq <- readRDS("~/data/pipeline/kazakhstan-experimental/pacbio/ITS9MUN/ITS2/ecm_physeq.Rdata")
#Checking again if there are strange Phylum, Class, Order or families
#sample_data(physeq)

if(!file.exists(paste0(RDS_dir, "/001_filtered_PSD.RDS"))){
  
  ################################################################################
  #0 Loading extra metadata
  message("Loading complementary metadata to include into the analysis")
  metadata <-
    as.data.frame(
      read.csv(
        "~/METADATA/metadata_sl-kazakhstan23.csv",header = T)
    )
  # metadata$Code <-
  #   sub(pattern = "[-]",
  #       replacement = "_",
  #       x = metadata$Code)
  # metadata$Code <-
  #   sub(pattern = "[-]",
  #       replacement = "_",
  #       x = metadata$Code)
  # metadata$...16 <- NULL
  # metadata$...17 <- NULL
  # metadata$...18 <- NULL
  # metadata$`AM/ECM`[which(metadata$`AM/ECM` == "oak")] <- "ECM"
  # metadata <-
  #   metadata[, c("Code",
  #                "Plot",
  #                "Individual",
  #                "Replica",
  #                "AM/ECM",
  #                "Height",
  #                "DBH",
  #                "Altitude")]
  # colnames(metadata)[5] <- "Fungi_system"
  # metadata$Height <- as.numeric(metadata$Height)
  # 
  ################################################################################
  #Fixing species as Genus + species
  
  message("Fixing species as Genus + species and removing sp.")
  
  tax <- as.data.frame(tax_table(physeq))
  #get unique species
  UT <- unique(as.character(tax[,7]))
  #getting taxa as Russula_sp or Russula Genus and replace by "?"
  target_positions <- grep("(_sp$| Genus$)", UT)
  UT[target_positions] <- "?"
  UT <- UT[which(UT!="?")]
  tax[,7][!tax[,7] %in% UT] <- "?"
  tax_table(physeq) <- as.matrix(tax)
  physeq <- microViz::tax_fix(physeq,unknowns = c("?"))
  # for(i in 1:nrow(tax)){
  #   if(tax$Genus[[i]] =="?" &
  #      tax$Species[[i]]=="?"
  #   ){
  #     tax$Species[[i]] <- "?"
  #   } else if(
  #     tax$Genus[[i]] !="?" &
  #     tax$Species[[i]]=="?"
  #   ){
  #     tax$Species[[i]] <- "?"
  #   } else if(
  #     tax$Genus[[i]] !="?" &
  #     tax$Species[[i]]!="?"
  #   ){
  #     tax$Species[[i]] <- paste0(tax$Genus[[i]]," ",tax$Species[[i]])
  #   }
  # };rm(i) 
  # 
  # physeq@tax_table@.Data[,7]  <- as.character(tax$Species)
  # rm(tax)
  # 1. Extract numeric ASV count table
  message("Starting analysis: Filtering Phyloseq object")
  
  asv_tab <- otu_table(physeq)
  if (taxa_are_rows(physeq)) {
    asv_tab <- t(asv_tab)
  }
  
  # 2. Convert to presence/absence (logical matrix)
  asv_bin <- asv_tab #asv_tab > 0
  
  # 3. Convert to data.frame if needed
  asv_bin_df <- as.data.frame(asv_bin)
  #converting numbers to 1 to get OTU in negative samples
  asv_bin_df[asv_bin_df > 0] <- 1
  #using metadata associated to samples table
  samp_data <- as.data.frame(sample_data(physeq))
  
  #get samples with zero counts and remove them
  physeq4 <-
    phyloseq::prune_samples(samples = names(phyloseq::sample_sums(physeq)[sample_sums(physeq) >
                                                                             0]),
                            x = physeq)
  #Remove samples with less than 1000 reads
   physeq4 <-
  phyloseq::prune_samples(samples = names(phyloseq::sample_sums(physeq4)[sample_sums(physeq4) >=
                                                                           1000]),
                          x = physeq4)

  #Prune taxa with zero counts
  physeq4 <- phyloseq::prune_taxa(taxa_sums(physeq4) > 0, physeq4)
  
  #Adding extra information to sample data
  #using metadata associated to samples table
  samp_data <- as.data.frame(as.matrix(sample_data(physeq4)))
  samp_data$id <- row.names(samp_data)
  samp_data <-
    dplyr::left_join(x = samp_data,
                     y = metadata,
                     by = c("id" = "SampleID"))
  row.names(samp_data) <- samp_data$id
  sample_data(physeq4) <- samp_data
  ################################################################################
  #Saving Phyloseq object
  message("Phyloseq object filtered obtained...Saving to RDS subfolder")
  saveRDS(physeq4, paste0(RDS_dir, "/001_filtered_PSD.RDS"))
} else {
  message("loading 001_filtered_PSD.RDS")
  physeq4 <- readRDS(paste0(RDS_dir, "/001_filtered_PSD.RDS"))
}

write.csv(data.frame(sample_data(physeq4)),paste0(csv_dir, "/", "000_Sample_data.csv"))
################################################################################
################################################################################
################################################################################
################################################################################
#9.Obtaining rarefied phyloseq object. Two options:"iNEXT","ggrare"

 if(!file.exists(paste0(RDS_dir, "/002_rarefied_PSD_4.RDS"))){
   
   # ggrare(physeq4)
  # ps_rare <- accumulation_curve_function(physeq4, option = "iNEXT")
  # ps_rare2 <- accumulation_curve_function(physeq4, option = "plateau")
  ps_rare3 <- accumulation_curve_function(physeq4, option = "SRS")
  ps_rare <- ps_rare3
  
 } else {
   ps_rare <- readRDS(paste0(RDS_dir, "/002_rarefied_PSD_4.RDS"))
 }
################################################################################
################################################################################
################################################################################
################################################################################

################################################################################
#10 Alpha richness analysis
message("Getting CSVs of Alpha diversity")
# #obtaining richness for non rarefied Phyloseq object
# Alpha <- phyloseq::estimate_richness(physeq4)
# Alpha$pielou <- evenness(physeq4, index = "pielou")[, 1]
# calc <- ChaoRichness(x = otu_table(phy_rare[[i]]), datatype = "abundance", conf = 0.95)
# 
# #Adding sample data
# Alpha <- cbind(sample_data(physeq4), Alpha)
# #obtaining richness for non rarefied Phyloseq object
Alpha_rare <- phyloseq::estimate_richness(ps_rare)
Alpha_rare$pielou <- evenness(ps_rare, index = "pielou")[, 1]
calc <- ChaoRichness(x = otu_table(ps_rare), datatype = "abundance", conf = 0.95)
Alpha_rare$Chao1 <- calc$Estimator
Alpha_rare$se.chao1 <- calc$Est_s.e.

#Adding sample data
Alpha_rare <- cbind(sample_data(ps_rare), Alpha_rare)

message("Saving Alpha diversity indexes CSV files")
# write.csv(
#   Alpha,
#   paste0(csv_dir, "/", "002_Alpha.csv"),
#   na = "",
#   row.names = F,
#   quote = T
# )
write.csv(
  Alpha_rare,
  paste0(csv_dir, "/", "002_Alpha_rarefaction.csv"),
  na = "",
  row.names = F,
  quote = T
)

################################################################################
#10 Alpha richness analysis
message("Getting CSVs of Alpha diversity")
# #obtaining richness for non rarefied Phyloseq object
# Alpha <- phyloseq::estimate_richness(physeq4)
# Alpha$pielou <- evenness(physeq4, index = "pielou")[, 1]
# calc <- ChaoRichness(x = otu_table(phy_rare[[i]]), datatype = "abundance", conf = 0.95)
# 
# #Adding sample data
# Alpha <- cbind(sample_data(physeq4), Alpha)
# #obtaining richness for non rarefied Phyloseq object
Alpha_N <- phyloseq::estimate_richness(physeq4)
Alpha_N$pielou <- evenness(physeq4, index = "pielou")[, 1]
calc <- ChaoRichness(x = otu_table(physeq4), datatype = "abundance", conf = 0.95)
Alpha_N$Chao1 <- calc$Estimator
Alpha_N$se.chao1 <- calc$Est_s.e.

#Adding sample data
Alpha_N <- cbind(sample_data(physeq4), Alpha_N)

message("Saving Alpha diversity indexes CSV files")
# write.csv(
#   Alpha,
#   paste0(csv_dir, "/", "002_Alpha.csv"),
#   na = "",
#   row.names = F,
#   quote = T
# )
write.csv(
  Alpha_N,
  paste0(csv_dir, "/", "002_Alpha.csv"),
  na = "",
  row.names = F,
  quote = T
)
################################################################################
# message("Plotting boxplots of Alpha diversity for rarefied and filtered Phyloseq object")
#Plot boxplots with Alpha indexes and including Wilcoxon tests
# p_r_p4 <-
#   plot_richness(
#     physeq4,
#     x = "Fungi_system",
#     color = "Fungi_system",
#     measures = c("Shannon", "Observed", "InvSimpson"),
#     title = "Non-rarefied"
#   ) +
#   ylab("") +
#   geom_boxplot(aes(fill = Fungi_system), alpha = 0.7) +
#   stat_compare_means(method = "wilcox.test",
#                      size = 8,
#                      aes(label = paste0("p = ", after_stat(p.format)))) +
#   ggsci::scale_color_jco() +   # jco palette for lines/points
#   ggsci::scale_fill_jco() +    # jco palette for boxplots
#   theme(
#     legend.position = "None",
#     strip.text = element_text(size = 16),
#     axis.title = element_text(size = 14),
#     axis.text = element_text(size = 14),
#     axis.title.x = element_text(size = 18),
#     axis.title.y = element_text(size = 18),
#     axis.text.x = element_text(size = 16),
#     axis.text.y = element_text(size = 22),
#     plot.title = element_text(size = 24, face = "bold", hjust = 0.5)
#   )

# my_comparisons <- list(
#   c("AM", "ECM"))
# 
# 
# p_r_p4_rare <-
#   plot_richness(
#     ps_rare,
#     x = "Fungi_system",
#     color = "Fungi_system",
#     measures = c("Shannon", "Observed", "InvSimpson"),
#     title = ""#"Rarefied to twice the minimum size"
#   ) +
#   ylab("") +
#   geom_boxplot(aes(fill = Fungi_system), alpha = 0.7) +
#   stat_compare_means(comparisons = my_comparisons,
#                      method = "wilcox.test",
#                      size = 8,
#                      aes(label = paste0("p = ", after_stat(p.format)))) +
#   ggsci::scale_color_jco() +
#   ggsci::scale_fill_jco() +
#   theme(
#     legend.position = "None",
#     strip.text = element_text(size = 16),
#     axis.title = element_text(size = 14),
#     axis.text = element_text(size = 14),
#     axis.title.x = element_text(size = 18),
#     axis.title.y = element_text(size = 18),
#     axis.text.x = element_text(size = 16),
#     axis.text.y = element_text(size = 22),
#     plot.title = element_text(size = 24, face = "bold", hjust = 0.5)
#   )
# 
# #p_Alpha <- gridExtra::grid.arrange(p_r_p4, p_r_p4_rare)
# 
# message("Alpha indicators...Saved in graphics subfolder")
# ggsave(
#   paste0(graph_dir, "/", "002_Alpha_boxplots.pdf"),
#   p_Alpha,
#   width = 20,
#   height = 25,
#   units = "in",
#   dpi = 600
# )

################################################################################
################################################################################

# ggsave(
#   paste0(graph_dir, "/", "002_Alpha_boxplots_rare.pdf"),
#   p_r_p4_rare,
#   width = 20,
#   height = 25,
#   units = "in",
#   dpi = 600
# )
################################################################################
#Using microbiome to transform in compositional (relative abundance to use in beta divesity)
physeq_clr <- microbiome::transform(ps_rare, "compositional")
physeq_clr2 <- microbiome::transform(physeq4, "clr")
################################################################################
#Trying reggression of Alpha and diameter at breast height (DBH)

#Alpha: It is  a data.frame with samples data and Alpha richness calculated
#col_fill_var: grouping variable for samples such as fungi system or ecosystems
#y_var: Variable to perform the regression
#saved_name: prefix to save files and recognize file in the graphics dir
# DBH <-
#   Regression_alpha_function(
#     Alpha = Alpha_rare,
#     col_fill_var = "Fungi_system",
#     y_var = "DBH",
#     saved_name = "rarefaction"
#   )
# 
# #Trying reggression of Alpha and diameter at tree height
# Height <-
#   Regression_alpha_function(
#     Alpha = Alpha_rare,
#     col_fill_var = "Fungi_system",
#     y_var = "Height",
#     saved_name = "rarefaction"
#   )

################################################################################
################################################################################
################################################################################
################################################################################
#BETA DIVERSITY



# Perform PCoA ordination using Bray-Curtis distance
# psd5.mds.bray <-
#   ordinate(physeq_clr, method = "PCoA", distance = "bray")
# evals <- psd5.mds.bray$values$Eigenvalues
# 
# # Create the ordination plot and apply the 'jco' color palette
# pord2 <-
#   plot_ordination(physeq_clr, psd5.mds.bray, color = "Fungi_system", shape = "Fungi_system") +
#   labs(col = "Fungi System") +
#   geom_point(size = 8) +  # Increase point size
#   
#   ggsci::scale_color_jco() +        # Apply jco palette for color
#   coord_fixed(sqrt(evals[2] / evals[1])) +
#   #theme_pubr(base_size = 20) +  # Increase all base text size
#   theme(
#     axis.title = element_text(size = 22),
#     axis.text = element_text(size = 20),
#     strip.text = element_text(size = 20),
#     legend.title = element_text(size = 20),
#     legend.text = element_text(size = 18)
#   )
# # Show the plot
# message("Beta diversity...Saved in graphics subfolder")
# ggsave(
#   paste0(graph_dir, "/", "004_Beta_PCoA_Bray.pdf"),
#   pord2,
#   width = 20,
#   height = 25,
#   units = "in",
#   dpi = 600
# )
################################################################################
# Perform PCoA ordination using Bray-Curtis distance
psd5.mds.euc2 <-
  ordinate(physeq_clr2, method = "PCoA", distance = "euclidean")
evals <- psd5.mds.euc2$values$Eigenvalues

# Create the ordination plot and apply the 'jco' color palette
pord2 <-
   plot_ordination(physeq_clr2, psd5.mds.euc2) + #, color = "Fungi_system", shape = "Fungi_system") +
  # labs(col = "Fungi System") +
  geom_point(size = 8) +  # Increase point size
  
  ggsci::scale_color_jco() +        # Apply jco palette for color
  coord_fixed(sqrt(evals[2] / evals[1])) +
  #theme_pubr(base_size = 20) +  # Increase all base text size
  theme(
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 20),
    strip.text = element_text(size = 20),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18)
  )
# Show the plot
message("Beta diversity...Saved in graphics subfolder")
ggsave(
  paste0(graph_dir, "/", "004_Beta_PCoA_EUC.pdf"),
  pord2,
  width = 20,
  height = 25,
  units = "in",
  dpi = 600
)

################################################################################
################################################################################
#perform PERMANOVA
##USING NO RAREFIED OBJECT + CLR
# Additionally, our betadisper results are not significant,
# meaning we cannot reject the null hypothesis that our groups have the same dispersions
# dist = phyloseq::distance(physeq_clr2, method="euclidean")
# metadata <- data.frame(sample_data(physeq_clr2))
# test.adonis <- adonis2(dist ~ Fungi_system, data = metadata,permutations = 10000,parallel = T)
# test.adonis <- as.data.frame(test.adonis)
# test.adonis$interpretation <- NA
# if(test.adonis$`Pr(>F)`[1]<0.05){
#   test.adonis$interpretation[1] <- "Differences of groups"  
# }
# write.csv(test.adonis,paste0(csv_dir,"/","004_PERMANOVA_CLR.csv"))
# aov_beta <- as.data.frame(anova(betadisper(dist, metadata$Fungi_system)))
# aov_beta$interpretation <- NA
# if(aov_beta$`Pr(>F)`[1]>0.05){
#   aov_beta$interpretation[1] <- "No reject null hypothesis (same dispersions)"  
# }
# write.csv(aov_beta,paste0(csv_dir,"/","004_BETADISPERSER_CLR.csv"))
# ################################################################################
# # #perform PERMANOVA
# # ##USING RAREFIED OBJECT
# # # Additionally, our betadisper results are not significant,
# # # meaning we cannot reject the null hypothesis that our groups have the same dispersions
# # dist = phyloseq::distance(physeq_clr, method="bray")
# # metadata <- data.frame(sample_data(physeq_clr))
# # test.adonis <- adonis2(dist ~ Fungi_system, data = metadata,permutations = 1000,parallel = T)
# # test.adonis <- as.data.frame(test.adonis)
# # test.adonis$interpretation <- NA
# # if(test.adonis$`Pr(>F)`[1]<0.05){
# #   test.adonis$interpretation[1] <- "Differences of groups"  
# # }
# # write.csv(test.adonis,paste0(csv_dir,"/","004_PERMANOVA_RA.csv"))
# # aov_beta <- as.data.frame(anova(betadisper(dist, metadata$Fungi_system)))
# # aov_beta$interpretation <- NA
# # if(aov_beta$`Pr(>F)`[1]>0.05){
# #   aov_beta$interpretation[1] <- "No reject null hypothesis (same dispersions)"  
# # }
# # write.csv(aov_beta,paste0(csv_dir,"/","004_BETADISPERSER_RA.csv"))
# ################################################################################
# ################################################################################
# ################################################################################
#ABUNDANCE BARPLOTS
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

p1 <-
  ps_rare  #phyloseq::merge_samples(ps_rare, group = "Fungi_system")
p1 <- tax_fix(p1, unknowns = c("?"))

# #phylum barplot
bp1 <- microViz::comp_barplot(
  tax_fix(p1),
  tax_level = "Order",
  bar_outline_colour = NA,
  sample_order = "bray",
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

bp4 <- microViz::comp_barplot(
  tax_fix(p1),
  tax_level = "Species",
  n_taxa = 10,
  bar_outline_colour = NA,
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
                                    bp2, bp3, bp4)


ggsave(
  paste0(graph_dir, "/", "005_BARPLOT_top10_Rare.pdf"),
  bp_total,
  width = 20,
  height = 25,
  units = "in",
  dpi = 600
)

################################################################################
#ABUNDANCE BARPLOTS
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

p1 <-
  physeq4  #phyloseq::merge_samples(ps_rare, group = "Fungi_system")
p1 <- tax_fix(p1, unknowns = c("?"))

# View(tax_table(ps_rare))

# #phylum barplot
bp1 <- microViz::comp_barplot(
  tax_fix(p1),
  tax_level = "Order",
  bar_outline_colour = NA,
  sample_order = "bray",
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

bp4 <- microViz::comp_barplot(
  tax_fix(p1),
  tax_level = "Species",
  n_taxa = 10,
  bar_outline_colour = NA,
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
                                    bp2, bp3, bp4)


ggsave(
  paste0(graph_dir, "/", "005_BARPLOT_top10.pdf"),
  bp_total,
  width = 20,
  height = 25,
  units = "in",
  dpi = 600
)

#Save RData for further steps
save.image(paste0(RDS_dir,"/","003_RData.RData"))


######
