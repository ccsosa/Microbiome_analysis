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
  Hmisc,
  corrplot,
  agricolae
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
source("~/SCRIPTS/INTERNAL_CODE/000_Prepare_dirs.R")
source("~/SCRIPTS/INTERNAL_CODE/001_Accumulation_Species_threshold.R")
source("~/SCRIPTS/INTERNAL_CODE/002_regression_alpha_vars.R")
source("~/SCRIPTS/INTERNAL_CODE/003_Differentially_abundant_taxon.R")
source("~/SCRIPTS/INTERNAL_CODE/000_Composition_plot.R")
source("~/SCRIPTS/INTERNAL_CODE/004_cor_taxa.R")
source("~/SCRIPTS/MAIN_CODE_AUXILIAR_SCRIPTS/000_Alpha_plots.R")
source("~/SCRIPTS/MAIN_CODE_AUXILIAR_SCRIPTS/000_Beta_plots.R")
source("~/SCRIPTS/MAIN_CODE_AUXILIAR_SCRIPTS/000_intersected_Species_per_column.R")
source("~/SCRIPTS/MAIN_CODE_AUXILIAR_SCRIPTS/000_summary_table.R")
source("~/SCRIPTS/MAIN_CODE_AUXILIAR_SCRIPTS/000_Corrplot.R")
source("~/SCRIPTS/MAIN_CODE_AUXILIAR_SCRIPTS/Alpha_soils_cor.R")


#Home folder
dir <- "~"
#how to name the results
marker="JGICOL_ITS2_EUKARYOME"
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
physeq <- readRDS("~/00101_20250411JC1K0Y/r_output/ITS2_eukaryome/ecm_physeq.Rdata")
#Saving sample data to double check soil data
# write.csv(data.frame(sample_data(physeq)),"PHYSEQ_EUKARYOME.csv",row.names = T)
if(!file.exists(paste0(RDS_dir, "/001_filtered_PSD.RDS"))){
  
  ################################################################################
  #0 Loading extra metadata
  message("Loading complementary metadata to include into the analysis")
  metadata <-
    as.data.frame(
      readxl::read_xlsx(
        "~/METADATA/JGI_PHYSEQ_METADATA.xlsx",
        sheet = "JGI_PHYSEQ",
        col_names = TRUE,
        # skip = 1,
        na = "NA"
      )
    )
  row.names(metadata) <- metadata$SampleID

  ################################################################################
  #Fixing species as Genus + species
  
  message("Fixing species as Genus + species and removing sp.")
  
  tax <- as.data.frame(tax_table(physeq))
  for(i in 1:nrow(tax)){
    if(tax$Genus[[i]] =="?" &
       tax$Species[[i]]=="?"
    ){
      tax$Species[[i]] <- "?"
    } else if(
      tax$Genus[[i]] !="?" &
      tax$Species[[i]]=="?"
    ){
      tax$Species[[i]] <- "?"
    } else if(
      tax$Genus[[i]] !="?" &
      tax$Species[[i]]!="?"
    ){
      tax$Species[[i]] <- paste0(tax$Genus[[i]]," ",tax$Species[[i]])
    }
  };rm(i)
  
  #get unique species
  UT <- unique(as.character(tax[,7]))
  #getting taxa as Russula_sp or Russula Genus and replace by "?"
  target_positions <- grep("(_sp$| Genus$)", UT)
  UT[target_positions] <- "?"
  UT <- UT[which(UT!="?")]
  tax[,7][!tax[,7] %in% UT] <- "?"
  tax_table(physeq) <- as.matrix(tax)
  physeq <- microViz::tax_fix(physeq,unknowns = c("?"))
  tax <- data.frame(tax_table(physeq))
  tax$Species <- sub(pattern = "Genus",replacement = "sp.",tax$Species)
  tax_table(physeq)[,7] <- tax$Species
  # 
  # physeq@tax_table@.Data[,7]  <- as.character(tax$Species)
  # rm(tax)
  # 1. Extract numeric ASV count table
  message("Starting analysis: Filtering Phyloseq object")
  
  samp_data <- as.data.frame(sample_data(physeq))

  physeq4 <-
    phyloseq::prune_samples(samples = row.names(samp_data[which(samp_data$Sample_or_Control!="Control")]),
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
  samp_data$SampleID <- row.names(samp_data)
  samp_data <-
    dplyr::left_join(x = samp_data,
                     y = metadata,
                     by = c("SampleID" = "SampleID"))
  row.names(samp_data) <- samp_data$SampleID
  sample_data(physeq4) <- samp_data
  ################################################################################
  #Saving Phyloseq object
  message("Phyloseq object filtered obtained...Saving to RDS subfolder")
  saveRDS(physeq4, paste0(RDS_dir, "/001_filtered_PSD.RDS"))
  write.csv(data.frame(sample_data(physeq4)),paste0(csv_dir, "/", "000_Sample_data.csv"))
} else {
  message("loading 001_filtered_PSD.RDS")
  physeq4 <- readRDS(paste0(RDS_dir, "/001_filtered_PSD.RDS"))
  samp_data <- sample_data(physeq4)
}

################################################################################
################################################################################
################################################################################
################################################################################
#Obtaining rarefied phyloseq object. Two options:"iNEXT","ggrare"

if(!file.exists(paste0(RDS_dir, "/002_rarefied_PSD_4.RDS"))){

  # ggrare(physeq4)
 # ps_rare <- accumulation_curve_function(physeq4, option = "iNEXT")
 # ps_rare2 <- accumulation_curve_function(physeq4, option = "plateau")
 ps_rare3 <- accumulation_curve_threshold_function(physeq4,
                                                   threshold_label = "20000",
                                                   option = "SRS",Cmin = 20000)
 saveRDS(ps_rare3, paste0(RDS_dir, "/002_rarefied_PSD_4.RDS"))
 
 # ps_rare <- ps_rare3

} else {
  ps_rare3 <- readRDS(paste0(RDS_dir, "/002_rarefied_PSD_4.RDS"))
}

################################################################################
################################################################################
################################################################################
################################################################################
#ALPHA
Alpha <- alpha_plots(physeq4,ps_rare3,csv_dir)
################################################################################
x_kr <- kruskal(Alpha$Alpha$Observed,Alpha$Alpha$Site_name.x,
                p.adj = "fdr",group = F,alpha=0.05)

comp <- x_kr$comparison
comp$Signif.


x_kr_rare <- kruskal(Alpha$Alpha_rare$Observed,Alpha$Alpha_rare$Site_name.x,
                p.adj = "fdr",group = F,alpha=0.05)

comp_rare <- x_kr_rare$comparison
comp_rare$Signif.


################################################################################
padj_threshold <- 0.05
indexes <- c("Observed","Chao1","ACE","Shannon","Simpson","InvSimpson","pielou")
vars <- c(34:45,47:73)
data <- Alpha$Alpha

x <- cor_alpha_function(data,indexes,vars,padj_threshold)
x$dataset <- "ALL"
x <- x[complete.cases(x),]

# data_AM <- Alpha$Alpha[which(Alpha$Alpha$Habit=="AM"),]
# x_AM <- cor_alpha_function(data_AM,indexes,vars,padj_threshold)
# x_AM$dataset <- "AM"
# x_AM <- x_AM[complete.cases(x_AM),]
# 
# data_ECM <- Alpha$Alpha[which(Alpha$Alpha$Habit=="ECM"),]
# x_ECM <- cor_alpha_function(data_ECM,indexes,vars,padj_threshold)
# x_ECM$dataset <- "ECM"
# x_ECM <- x_ECM[complete.cases(x_ECM),]

################################################################################
data <- Alpha$Alpha
ggscatter(data, x = "DBH", y = "Observed", color = "Site_name.x", #"Habit",#"Site_name.x",
          ylab = "Diameter at Breast Height (DBH)",
          
          xlab = "Observed richness")+
  # línea de regresión general
  geom_smooth(method = "lm", se = FALSE, color = "black", formula = y ~ x) +
  ggsci::scale_color_jco() +   # jco palette for lines/points
  ggsci::scale_fill_jco() +    # jco palette for boxplots
  
  # mostrar correlación en lugar de fórmula
  stat_cor(method = "spearman",
           # label.x = min(xy$var), label.y = max(xy$freq),
           # aes(x = freq, y = var),
           color = "black") +
  theme(legend.title = element_blank())

################################################################################
names(data) <- make.names(names(data))

ggscatter(data, x = "DBH", 
          y = "ACE", 
          color = "Habit", #"Habit",#"Site_name.x",
          ylab = "Iron (mg/kg)",
          
          xlab = "Observed richness")+
  # línea de regresión general
  geom_smooth(method = "lm", se = FALSE, color = "black", formula = y ~ x) +
  # mostrar correlación en lugar de fórmula
  stat_cor(method = "spearman",
           # label.x = min(xy$var), label.y = max(xy$freq),
           # aes(x = freq, y = var),
           color = "black") +
  theme(legend.title = element_blank())

# length(sample_sums(physeq4)[sample_sums(physeq4)<10000])
# length(sample_sums(physeq4)[sample_sums(physeq4)<40000])
# 
# x_samp <- sample_data(physeq4)
# x <-sample_sums(physeq4)[sample_sums(physeq4)<40000]
#testing loss of samples
# x_samp_2 <- x_samp[x_samp$SampleID %in% names(x),]
# tapply(x_samp_2$Habit,x_samp_2$Habit,length)
#if cor has p value < alpha value (Thus, rarefy!)
# cor.test(sample_sums(physeq4),Alpha$Observed)
# plot(sample_sums(physeq4),Alpha$Observed)

################################################################################
jco_palette <- ggsci::pal_jco()(7)
jco_more <- colorRampPalette(jco_palette)(length(unique(samp_data$Site_name.x))) # generar 20 colores
################################################################################
x <- Alpha$Alpha[,c(36:45,47:75,77,79:82)]
corrplot_function(x,graph_dir,csv_dir)
PCA <- FactoMineR::PCA(x,scale.unit = T,graph = F,ncp = 5)
PCA_viz <-fviz_pca_biplot(
  PCA,
  geom.ind = "point",     # Solo puntos
  col.ind = Alpha$Alpha$Site_name.x,        # Color por grupo
  # palette = c("#00AFBB", "#E7B800", "#FC4E07"), # Paleta personalizada
  addEllipses = TRUE,     # Agregar elipses
  ellipse.type = "convex", # "confidence" o "t"
  legend.title = "",repel = T,col.var = "black"
) +
  scale_color_manual(values = jco_more) +
  scale_fill_manual(values = jco_more) +
  
  theme_minimal() +
  theme(
    legend.position = "top",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

ggsave(
  paste0(graph_dir, "/", "003_PCA_sites.pdf"),
  PCA_viz,
  width = 20,
  height = 25,
  units = "in",
  dpi = 600
)
###############################################################################
#PCA (RARE)

x <- Alpha$Alpha_rare[,c(36:45,47:75,77,79:82)]
PCA <- FactoMineR::PCA(x,scale.unit = T,graph = F,ncp = 5)
PCA_viz <-fviz_pca_biplot(
  PCA,
  geom.ind = "point",     # Solo puntos
  col.ind = Alpha$Alpha_rare$Site_name.x,        # Color por grupo
  # palette = c("#00AFBB", "#E7B800", "#FC4E07"), # Paleta personalizada
  addEllipses = TRUE,     # Agregar elipses
  ellipse.type = "convex", # "confidence" o "t"
  legend.title = "",repel = T,col.var = "black"
) +
  scale_color_manual(values = jco_more) +
  scale_fill_manual(values = jco_more) +
  
  theme_minimal() +
  theme(
    legend.position = "top",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

ggsave(
  paste0(graph_dir, "/", "003_PCA_sites_RARE.pdf"),
  PCA_viz,
  width = 20,
  height = 25,
  units = "in",
  dpi = 600
)

################################################################################
################################################################################
################################################################################
#Using microbiome to transform in compositional (relative abundance to use in beta divesity)
physeq_clr2 <- microbiome::transform(physeq4, "compositional")
################################################################################
################################################################################
################################################################################
beta <- beta_plots(physeq_clr2,csv_dir)
# ################################################################################
# ################################################################################
# ################################################################################
#Habit
bp1 <- compo_plot(physeq=physeq4,
                       level="Order",
                       ntop=10,
                       column="Habit")
  
bp2 <- compo_plot(physeq=physeq4,
                  level="Family",
                  ntop=10,
                  column="Habit")

bp3 <- compo_plot(physeq=physeq4,
                  level="Genus",
                  ntop=10,
                  column="Habit")
bp4 <- compo_plot(physeq=physeq4,
                  level="Species",
                  ntop=10,
                  column="Habit")
bp_total <- gridExtra::arrangeGrob(bp1,
                                    bp2, bp3, bp4)


png(
  filename = file.path(graph_dir, "005_BARPLOT_top10_MT.png"),
  type = "cairo",
  width = 10,
  height = 10,
  units = "in",
  res = 600
)
grid::grid.draw(bp_total)
dev.off()
################################################################################
#Sites
bp1_S <- compo_plot(physeq=physeq4,
                  level="Order",
                  ntop=10,
                  column="Site_name.x")

bp2_S <- compo_plot(physeq=physeq4,
                  level="Family",
                  ntop=10,
                  column="Site_name.x")

bp3_S <- compo_plot(physeq=physeq4,
                  level="Genus",
                  ntop=10,
                  column="Site_name.x")
bp4_S <- compo_plot(physeq=physeq4,
                  level="Species",
                  ntop=10,
                  column="Site_name.x")
bp_total_S <- gridExtra::arrangeGrob(bp1_S,
                                   bp2_S, bp3_S, bp4_S)


png(
  filename = file.path(graph_dir, "005_BARPLOT_top10_SITES.png"),
  type = "cairo",
  width = 10,
  height = 10,
  units = "in",
  res = 600
)
grid::grid.draw(bp_total_S)
dev.off()
################################################################################
cor_table <- cor_var_function(ps = physeq4,padj_threshold=0.05)
write.csv(
  cor_table,
  paste0(csv_dir, "/", "006_Correlation_Table_general_F.csv"),
  na = "",
  row.names = F,
  quote = T
)

physeq_clr2_sp_sp <- microViz::tax_agg(physeq4,"Genus")
physeq_clr2_sp_sp <- transform(physeq_clr2_sp_sp,"compositional")
x <- data.frame(otu_table(physeq_clr2_sp_sp))
x <- x[which(row.names(x)=="Hortiboletus"),]
tol <- 1e-12
x[abs(x) < tol] <- NA
x <- data.frame(SampleID = names(x),
                freq = as.numeric(x))
y <- data.frame(sample_data(physeq_clr2))
y <- data.frame(SampleID = row.names(y),
                var =as.numeric(y$Zinc..mg.kg.),
                Site=y$Site_name.x,
                Habit = y$Habit
)

xy <-
  dplyr::left_join(x = x,
                   y = y,
                   by = c("SampleID" = "SampleID"))
xy <- xy[complete.cases(xy),]

xy <- xy[xy$var %in% na.omit(remove_outliers(xy$var)),]

ggscatter(xy, y = "freq", x = "var", color = "Site", #"Habit",#"Site_name.x",
          xlab = "Zinc (mg/Kg)",
          
          ylab = "Hortiboletus (relative abundance)")+
     ggsci::scale_color_jco() +   # jco palette for lines/points
     ggsci::scale_fill_jco() +    # jco palette for boxplots
  # línea de regresión general
  geom_smooth(method = "lm", se = FALSE, color = "black", formula = y ~ x) +
  # mostrar correlación en lugar de fórmula
  stat_cor(method = "spearman",
           # label.x = min(xy$var), label.y = max(xy$freq),
           # aes(x = freq, y = var),
           color = "black") +
  theme(legend.title = element_blank())

################################################################################
names(data) <- make.names(names(data))


################################################################################
AM_ps <- subset_samples(physeq4,Habit=="AM")
AM_ps <- phyloseq::prune_taxa(taxa_sums(AM_ps) > 0, AM_ps)

cor_table_AM <- cor_var_function(ps = AM_ps,padj_threshold=0.05)



# physeq_clr2_sp_sp <- microViz::tax_agg(physeq_clr2,"Genus")
#physeq_clr2_sp_sp <- physeq_clr2_sp
# AM_ps_RA <- transform(AM_ps,"compositional")
# x <- data.frame(otu_table(AM_ps_RA))
# x <- x[which(row.names(x)=="OTU937"),]
# tol <- 1e-12
# x[abs(x) < tol] <- NA  
# x <- data.frame(SampleID = names(x),
#                 freq = as.numeric(x))
# y <- data.frame(sample_data(AM_ps_RA))
# y <- data.frame(SampleID = row.names(y),
#                 var =as.numeric(y$Exchangeable.Calcium..mg.kg.
# 
# ))
# 
# xy <-
#   dplyr::left_join(x = x,
#                    y = y,
#                    by = c("SampleID" = "SampleID"))
# xy <- xy[complete.cases(xy),]
# xy <- xy[xy$var %in% remove_outliers(xy$var),]
# xy <- xy[xy$var<56.20,]
# cor.test(xy$freq,xy$var)
# plot(xy$freq,xy$var)

write.csv(
  cor_table_AM,
  paste0(csv_dir, "/", "006_Correlation_Table_AM.csv"),
  na = "",
  row.names = F,
  quote = T
)

####################################################################
ECM_ps <- subset_samples(physeq4,Habit=="ECM")
ECM_ps <- phyloseq::prune_taxa(taxa_sums(ECM_ps) > 0, ECM_ps)

cor_table_ECM <- cor_var_function(ps = ECM_ps,padj_threshold=0.05)
write.csv(
  cor_table_ECM,
  paste0(csv_dir, "/", "006_Correlation_Table_ECM.csv"),
  na = "",
  row.names = F,
  quote = T
)
####################################################################
column <- "Site_name.x"
levels <- c("Genus","Species","OTU")
name <- "Sites"
Upsetplot_compo(physeq4=physeq4,graph_dir = graph_dir,column=column,levels=levels,name=name)

column <- "Habit"
levels <- c("Genus","Species","OTU")
name <- "M_T"
Upsetplot_compo(physeq4=physeq4,graph_dir = graph_dir,column=column,levels=levels,name=name)

####################################################################
summary <- summary_table_function(physeq4)
  
#Save RData for further steps
save.image(paste0(RDS_dir,"/","Analysis_RData.RData"))
####################################################################
#load(paste0(RDS_dir,"/","Analysis_RData.RData"))

