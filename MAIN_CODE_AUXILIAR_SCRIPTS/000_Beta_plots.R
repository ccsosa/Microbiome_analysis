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


beta_plots <- function(physeq_clr2,csv_dir){
#BETA DIVERSITY
# Perform PCoA ordination using Bray-Curtis distance
psd5.mds.euc2 <-
  ordinate(physeq_clr2, method = "PCoA", distance = "bray")
evals <- psd5.mds.euc2$values$Eigenvalues

# Extract files
ord_df <- as.data.frame(plot_ordination(physeq_clr2, psd5.mds.euc2, justDF = TRUE))

#Extract sequencing technology
ord_df$Habit <- sample_data(physeq_clr2)$Habit
#Extracting labels
ord_df$SampleID <- sample_data(physeq_clr2)$SampleID

# Obtaining convex hull
hulls <- ord_df %>%
  group_by(Habit) %>%
  slice(chull(Axis.1, Axis.2))  # Convex hull para cada grupo



# Create the ordination plot and apply the 'jco' color palette
pord2 <-
  plot_ordination(physeq_clr2, psd5.mds.euc2, color = "Habit", shape = "Habit") +
  labs(col = "Habit") +
  geom_point(size = 8) +  # Increase point size
  geom_polygon(data = hulls, aes(x = Axis.1, y = Axis.2, fill = Habit), alpha = 0.2, color = NA) +
  
  # labs(col = "Fungi System") +
  # scale_color_manual(values = jco_more) +
  # scale_fill_manual(values = jco_more) +
  ggsci::scale_color_jco() +         # Color de puntos
  ggsci::scale_fill_jco() +          # Color de polígonos (igual que los puntos)
  
  coord_fixed(sqrt(evals[2] / evals[1])) +
  #theme_pubr(base_size = 20) +  # Increase all base text size
  theme(
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 20),
    strip.text = element_text(size = 20),
    legend.title = element_text(size = 0),
    legend.text = element_text(size = 18)
  )
pord2
# Show the plot
message("Beta diversity...Saved in graphics subfolder")
ggsave(
  paste0(graph_dir, "/", "004_Beta_PCoA_EUC_RA_MS.pdf"),
  pord2,
  width = 20,
  height = 25,
  units = "in",
  dpi = 600
)

################################################################################
#SITES

# BETA DIVERSITY
# Perform PCoA ordination using Bray-Curtis distance
psd5.mds.euc2 <-
  ordinate(physeq_clr2, method = "PCoA", distance = "bray")
evals <- psd5.mds.euc2$values$Eigenvalues

# Extract files
ord_df <- as.data.frame(plot_ordination(physeq_clr2, psd5.mds.euc2, justDF = TRUE))

#Extract sequencing technology
ord_df$Site_name.x <- sample_data(physeq_clr2)$Site_name.x
#Extracting labels
ord_df$SampleID <- sample_data(physeq_clr2)$SampleID

# Obtaining convex hull
hulls <- ord_df %>%
  group_by(Site_name.x) %>%
  slice(chull(Axis.1, Axis.2))  # Convex hull para cada grupo



# Create the ordination plot and apply the 'jco' color palette
pord2 <-
  plot_ordination(physeq_clr2, psd5.mds.euc2, color = "Site_name.x") +
  labs(col = "Site_name.x") +
  geom_point(size = 8) +  # Increase point size
  geom_polygon(data = hulls, aes(x = Axis.1, y = Axis.2, fill = Site_name.x), alpha = 0.2, color = NA) +
  
  # labs(col = "Fungi System") +
  scale_color_manual(values = jco_more) +
  scale_fill_manual(values = jco_more) +
  
  coord_fixed(sqrt(evals[2] / evals[1])) +
  #theme_pubr(base_size = 20) +  # Increase all base text size
  theme(
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 20),
    strip.text = element_text(size = 20),
    legend.title = element_text(size = 0),
    legend.text = element_text(size = 18)
  )
pord2
# Show the plot
message("Beta diversity...Saved in graphics subfolder")
ggsave(
  paste0(graph_dir, "/", "004_Beta_PCoA_EUC_RA_SITE.pdf"),
  pord2,
  width = 20,
  height = 25,
  units = "in",
  dpi = 600
)

################################################################################
#perform PERMANOVA
#USING NO RAREFIED OBJECT + CLR
# Additionally, our betadisper results are not significant,
# meaning we cannot reject the null hypothesis that our groups have the same dispersions
dist = phyloseq::distance(physeq_clr2, method="bray")
metadata <- data.frame(sample_data(physeq_clr2))
test.adonis <- adonis2(dist ~ Habit, data = metadata,permutations = 10000,parallel = T)
test.adonis <- as.data.frame(test.adonis)
test.adonis$interpretation <- NA
if(test.adonis$`Pr(>F)`[1]<0.05){
  test.adonis$interpretation[1] <- "Differences of groups"
}
write.csv(test.adonis,paste0(csv_dir,"/","004_PERMANOVA_RA_MS.csv"))
aov_beta <- as.data.frame(anova(betadisper(dist, metadata$Habit)))
aov_beta$interpretation <- NA
if(aov_beta$`Pr(>F)`[1]>0.05){
  aov_beta$interpretation[1] <- "No reject null hypothesis (same dispersions)"
}
write.csv(aov_beta,paste0(csv_dir,"/","004_BETADISPERSER_RA_MS.csv"))
# ################################################################################
# # #perform PERMANOVA
# # ##USING RAREFIED OBJECT
test.adonis <- adonis2(dist ~ Site_name.x, data = metadata,permutations = 10000,parallel = T)
test.adonis <- as.data.frame(test.adonis)
test.adonis$interpretation <- NA
if(test.adonis$`Pr(>F)`[1]<0.05){
  test.adonis$interpretation[1] <- "Differences of groups"
}
write.csv(test.adonis,paste0(csv_dir,"/","004_PERMANOVA_RA_SITE.csv"))
aov_beta <- as.data.frame(anova(betadisper(dist, metadata$Site_name.x)))
aov_beta$interpretation <- NA
if(aov_beta$`Pr(>F)`[1]>0.05){
  aov_beta$interpretation[1] <- "No reject null hypothesis (same dispersions)"
}
write.csv(aov_beta,paste0(csv_dir,"/","004_BETADISPERSER_RA_SITE.csv"))

return("BETA PLOTS DONE!")
}