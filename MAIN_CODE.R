################################################################################
#Load pacman to load complementary libraries
library(pacman);
#Load libraries or install if it is the case from CRAN and github respetively
pacman::p_load(phyloseq,microbiome,ggplot2,ggpubr,xlsx,readr,dplyr,compositions,ggrepel,iNEXT,gridExtra,ggsci)
pacman::p_load_gh("gauravsk/ranacapa", "gmteunisse/fantaxtic",
                  "kasperskytte/ampvis2","david-barnett/microViz@0.12.0"
                  )#,"twbattaglia/btools")
################################################################################
#Loading scripts
message("Loading complementary scripts")
source("~/SCRIPTS/000_Prepare_dirs.R")
source("~/SCRIPTS/001_Accumulation_Species.R")
source("~/SCRIPTS/002_regression_alpha_vars.R")
################################################################################
################################################################################
#Fixing a seed. Just in case
set.seed(1000)
################################################################################
#Load Phyloseq object
message("Loading Phyloseq object")
load("~/00101_20250411JC1K0Y/lotus2_report/ITS2/phyloseq.Rdata")
#Checking again if there are strange Phylum, Class, Order or families
#sample_data(physeq)
################################################################################
#0 Loading extra metadata
message("Loading complementary metadata to include into the analysis")
metadata <- as.data.frame(readxl::read_xlsx("~/METADATA/Database samples_JGI_20250612.xlsx",
                              sheet = "Soil",col_names = TRUE,skip = 1,na = "NA"))
metadata$Code <- sub(pattern = "[-]",replacement = "_",x = metadata$Code)
metadata$Code <- sub(pattern = "[-]",replacement = "_",x = metadata$Code)
metadata$...16 <- NULL
metadata$...17 <- NULL
metadata$...18 <- NULL
metadata$`AM/ECM`[which(metadata$`AM/ECM`=="oak")] <- "ECM"
metadata <- metadata[,c("Code","Plot","Individual","Replica","AM/ECM","Height",
                        "DBH","Altitude")]
colnames(metadata)[5] <- "Fungi_system"
metadata$Height <- as.numeric(metadata$Height)

################################################################################
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
asv_bin_df[asv_bin_df>0] <- 1 
#using metadata associated to samples table
samp_data <- as.data.frame(sample_data(physeq))

#4. getting negative controls to prune
negatives <- row.names(samp_data[which(samp_data$Sample_or_Control=="Control")])
asv_neg <- asv_bin_df[negatives,]
# asv_neg[3,] <- colSums(asv_neg,na.rm = T)
# #OTUs associated to negative samples
# colnames(asv_neg[,asv_neg[3,]>0])
# neg <- tax_table(physeq)[colnames(asv_neg[,asv_neg[3,]>0]),]
#remove sample data
#rm(samp_data)

#5. Prune negative controls
physeq2 <- phyloseq::prune_samples(samples=row.names(samp_data)[
  !row.names(samp_data) %in% negatives],
  x=physeq)

#6. Prune taxa with no counts
physeq3 <- phyloseq::prune_taxa(taxa_sums(physeq2) > 0, physeq2)
#get unique Phylum
#unique(tax_table(physeq3)[,"Phylum"])

filterPhyla = c("Annelida",
                "Apicomplexa",
                "Arthropoda",
                "Ascomycota",
                "Tardigrada",
                "Vertebrata",
                "Rotifera",
                "Mollusca",
                "Cnidaria",
                "Cercozoa",
                "Ciliophora",
                "Schizoplasmodiida",
                "Echinodermata",
                "Nematozoa",
                "Diatomea",
                "Platyhelminthes",
                "Ciliophora",
                "Porifera",
                "Kathablepharidae",
                "Phragmoplastophyta",
                "Protalveolata",
                "Nematoda",
                "Centrohelida",
                "Apusomonadidae",
                "Gracilipodida",
                "Rigifilida",
                "Katablepharidophyta",
                "Lycopodiophyta",
                "Eukaryota_phy_Incertae_sedis",
                "Fungi_phy_Incertae_sedis",
                "?"
                )
################################################################################
message("Filtering strange Phyla")

#7. Filter to remove strange Phylum
physeq3 <- phyloseq::subset_taxa(physeq3, !Phylum %in% filterPhyla)

#8. Filter to at least get samples for 10% of samples
physeq4 <- phyloseq::filter_taxa(physeq3, function(x) sum(x > 2) > (0.1*length(x)), TRUE)
#get samples with zero counts and remove them
physeq4 <- phyloseq::prune_samples(samples = names(phyloseq::sample_sums(physeq4)[sample_sums(physeq4)>0]),
                                   x =physeq4)
#Remove samples with less than 200 reads  
physeq4 <- phyloseq::prune_samples(samples = names(phyloseq::sample_sums(physeq4)[sample_sums(physeq4)>=200]),
                                   x =physeq4)
#Prune taxa with zero counts
physeq4 <- phyloseq::prune_taxa(taxa_sums(physeq4) > 0, physeq4)

#Adding extra information to sample data
#using metadata associated to samples table
samp_data <- as.data.frame(as.matrix(sample_data(physeq4)))
samp_data$id <- row.names(samp_data) 
samp_data <- dplyr::left_join(x = samp_data,y = metadata,by = c("id"="Code"))
row.names(samp_data) <- samp_data$id
sample_data(physeq4) <- samp_data
################################################################################
#Saving Phyloseq object
message("Phyloseq object filtered obtained...Saving to RDS subfolder")
saveRDS(physeq4,paste0(RDS_dir,"/001_filtered_PSD.RDS"))
################################################################################
#9.Obtaining rarefied phyloseq object. Two options:"iNEXT","ggrare"
ps_rare <- accumulation_curve_function(physeq4,option="iNEXT") 
################################################################################
#10 Alpha richness analysis
message("Getting CSVs of Alpha diversity")
#obtaining richness for non rarefied Phyloseq object
Alpha <- phyloseq::estimate_richness(physeq4)
Alpha$pielou <- evenness(physeq4,index = "pielou")[,1]
#Adding sample data
Alpha <-cbind(sample_data(physeq4),Alpha)
#obtaining richness for non rarefied Phyloseq object
Alpha_rare <- phyloseq::estimate_richness(ps_rare)
Alpha_rare $pielou <- evenness(ps_rare,index = "pielou")[,1]

#Adding sample data
Alpha_rare <-cbind(sample_data(ps_rare),Alpha_rare)

message("Saving Alpha diversity indexes CSV files")
write.csv(Alpha,paste0(csv_dir, "/", "002_Alpha.csv"),na = "",row.names = F,quote =T)
write.csv(Alpha_rare,paste0(csv_dir, "/", "002_Alpha_rarefaction.csv"),na = "",row.names = F,quote =T)

message("Plotting boxplots of Alpha diversity for rarefied and filtered Phyloseq object")
#Plot boxplots with Alpha indexes and including Wilcoxon tests
p_r_p4 <- plot_richness(physeq4, x = "Fungi_system", color = "Fungi_system",
                        measures = c("Shannon", "Observed", "InvSimpson"),
                        title = "Non-rarefied") +
  ylab("") +
  geom_boxplot(aes(fill = Fungi_system), alpha = 0.7) +
  stat_compare_means(
    method = "wilcox.test",
    size = 8,
    aes(label = paste0("p = ", after_stat(p.format)))
  ) +
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

p_r_p4_rare <- plot_richness(ps_rare, x = "Fungi_system", color = "Fungi_system",
                             measures = c("Shannon", "Observed", "InvSimpson"),
                             title = "Rarefied to twice the minimum size") +
  ylab("") +
  geom_boxplot(aes(fill = Fungi_system), alpha = 0.7) +
  stat_compare_means(
    method = "wilcox.test",
    size = 8,
    aes(label = paste0("p = ", after_stat(p.format)))
  ) +
  ggsci::scale_color_jco() +
  ggsci::scale_fill_jco() +
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

p_Alpha <- gridExtra::grid.arrange(p_r_p4, p_r_p4_rare)

message("Alpha indicators...Saved in graphics subfolder")
ggsave(
  paste0(graph_dir, "/", "002_Alpha_boxplots.pdf"), p_Alpha,
  width = 20, height = 25, units = "in",dpi = 600)

################################################################################
#Trying reggression of Alpha and diameter at breast height (DBH) 

#Alpha: It is  a data.frame with samples data and Alpha richness calculated
#col_fill_var: grouping variable for samples such as fungi system or ecosystems
#y_var: Variable to perform the regression
#saved_name: prefix to save files and recognize file in the graphics dir
DBH<- 
Regression_alpha_function(Alpha=Alpha_rare,
                          col_fill_var="Fungi_system",
                          y_var="DBH",
                          saved_name="rarefaction")
  
#Trying reggression of Alpha and diameter at tree height 
Height <- 
Regression_alpha_function(Alpha=Alpha_rare,
                          col_fill_var="Fungi_system",
                          y_var="Height",
                          saved_name="rarefaction")

################################################################################
#BETA DIVERSITY
#Using microbiome to transform in compositional (relative abundance to use in beta divesity)
physeq_clr <- microbiome::transform(physeq4, "compositional")


# Perform MDS ordination using Bray-Curtis distance
psd5.mds.bray <- ordinate(physeq_clr, method = "MDS", distance = "bray")
evals <- psd5.mds.bray$values$Eigenvalues

# Create the ordination plot and apply the 'jco' color palette
pord2 <- plot_ordination(physeq_clr, psd5.mds.bray, color = "Fungi_system", shape = "Fungi_system") +
  labs(col = "Fungi System") +
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




# Optional: Facet the plot by Fungi_system
#pord2 <- pord2 + facet_grid(. ~ Fungi_system)

# Show the plot

message("Beta diversity...Saved in graphics subfolder")
ggsave(
  paste0(graph_dir, "/", "004_Beta_NMDS_Bray.pdf"), pord2,
  width = 20, height = 25, units = "in",dpi = 600)

################################################################################
#BETA DIVERSITY


#Plotting richness per Fungi system




# correlation.table <- microbiome::associate(x=as.data.frame(t(physeq4@otu_table)), 
#                                            y=as.data.frame(meta(physeq4)[,c(19)]), 
#                                            method = "pearson",
#                                            mode = "table",
#                                            order=T,
#                                            p.adj.threshold = 0.05,
#                                            n.signif = 1,
#                                            p.adj.method = "fdr",
#                                            filter.self.correlations = T
#                                            )
# 
# correlation.table$X2 <- "DBH"
# hm_plot <- microbiome::heat(correlation.table, "X1", "X2", fill = "Correlation", 
#                             star = "p.adj", p.adj.threshold = 0.05,
#                             colours = c("darkred","red","white","blue","darkblue"),plot.values = F,star.size=40,legend.text = "Spearman"
#                             
# )+
#   theme(panel.background = element_rect(fill = "gray95"),
#         text=element_text(size=60),axis.text.x  = element_text(size=60,colour="black",angle = 90, hjust = 1),
#         axis.text.y  = element_text(size=60,colour="black")) +
#   #scale_fill_continuous(guide = guide_colorbar(direction = "horizontal")) +
#   theme(legend.position="right",legend.direction = "vertical",legend.key.size =  unit(1.5, "in"))


# ggsave(paste0(graph_dir,"/","003_heat_correlation_abund",".pdf"),hm_plot,
#        dpi=600,width =160,height=130,units = "cm",scale=1.2,limitsize = FALSE)














################################################################################






top15 <- top_taxa(physeqRarefied, n = 15, relative = T,
                      discard_other = T, other_label = "Other")

top15 <- name_na_taxa(top15, label = "", species = F, other_label = "Other")

fantaxtic_bar(top15, color_by = "Family", label_by = "Genus", facet_by = NULL, grid_by = NULL, other_color = "Grey") -> ptop15


plot_bar(physeqRarefied, x = "Sample", y = "Abundance", fill ="Phylum") +
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")


x<- merge_samples(x = physeqRarefied,group = "Fungi_system")

# Transformar a abundancia relativa
physeq_rel <- transform_sample_counts(physeqRarefied, function(x) x / sum(x))

# Graficar abundancia relativa por Phylum
# Graficar abundancia relativa agrupado por "Fungi_system"
plot_bar(physeq_rel, x = "Sample", y = "Abundance", fill = "Phylum") +
  geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position = "stack") +
  ylab("Relative Abundance") +
  facet_wrap(~Fungi_system, scales = "free_x") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))



####Plotting trajactories with InvSimpson
# 
# label_data <- Alph_list %>%
#   group_by(id) %>%
#   filter(simulated == "No rarefaction") %>%
#   ungroup()
# 
# 
# 
# p <- ggscatter(
#   data=Alph_list,
#   x="simulated",
#   y="InvSimpson",#"Height",
#   combine = FALSE,
#   merge = FALSE,
#   color = "Fungi_system",
#   fill = "Fungi_system",
#   palette = "jco",
#   shape = 19,
#   size = 2,
#   point = TRUE,
#   rug = FALSE,
#   title = NULL,
#   xlab = "Simulated sample size to rarify",
#   ylab = "Effective Number of Species",
#   facet.by = "Fungi_system",
#   panel.labs = NULL,
#   short.panel.labs = TRUE,
#   add = c("reg.line"),
#   conf.int = TRUE,
#   conf.int.level = 0.95,
#   fullrange = FALSE,
#   label = NULL,
#   font.label = c(12, "plain"),
#   font.family = "",
#   label.select = NULL,
#   repel = T,
#   label.rectangle = FALSE,
#   parse = FALSE,
#   cor.coef = FALSE,
#   cor.coeff.args = list(),
#   cor.method = "pearson",
#   cor.coef.coord = c(NULL, NULL),
#   cor.coef.size = 4,
#   ggp = NULL,
#   show.legend.text = NA,
#   ggtheme = theme_pubr()
# ) +
#   stat_cor(aes(color = Fungi_system)) +
#   geom_line(aes(group = id)) +
#   geom_text_repel(data = label_data, 
#                   aes(x = simulated, y = InvSimpson, label = id, color = Fungi_system),
#                   hjust = -0.1, 
#                   size = 3.2, 
#                   fontface = "italic",
#                   show.legend = FALSE,max.overlaps = 10000))
# 
# p <- ggpar(p, x.text.angle = 90)
# p
###############################################################################
################################################################################
# ####Plotting trajactories with library size
# p <- ggscatter(
#   data=Alph_list,
#   x="simulated",
#   y="library_size",#"Height",
#   combine = FALSE,
#   merge = FALSE,
#   color = "Fungi_system",
#   fill = "Fungi_system",
#   palette = "jco",
#   shape = 19,
#   size = 2,
#   point = TRUE,
#   rug = FALSE,
#   title = NULL,
#   xlab = "Simulated sample size to rarify",
#   ylab = "Library Size",
#   #facet.by = "Fungi_system",
#   panel.labs = NULL,
#   short.panel.labs = TRUE,
#   add = c("reg.line"),
#   conf.int = TRUE,
#   conf.int.level = 0.95,
#   fullrange = FALSE,
#   label = NULL,
#   font.label = c(12, "plain"),
#   font.family = "",
#   label.select = NULL,
#   repel = T,
#   label.rectangle = FALSE,
#   parse = FALSE,
#   cor.coef = FALSE,
#   cor.coeff.args = list(),
#   cor.method = "pearson",
#   cor.coef.coord = c(NULL, NULL),
#   cor.coef.size = 4,
#   ggp = NULL,
#   show.legend.text = NA,
#   ggtheme = theme_pubr()
# ) +
#   stat_cor(aes(color = Fungi_system)) +
#   geom_line(aes(group = id))
# p <- ggpar(p, x.text.angle = 90)
# p
# ###############################################################################
