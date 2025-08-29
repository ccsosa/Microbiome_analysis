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
  UpSetR,
  ggVennDiagram,
  ggplotify,
  aplot,
  qqman
)
pacman::p_load_gh(
  "gauravsk/ranacapa",
  "gmteunisse/fantaxtic",
  "kasperskytte/ampvis2",
  "david-barnett/microViz@0.12.0",
  "thomasp85/patchwork"
)#,"twbattaglia/btools")

# library(MicrobiotaProcess)
################################################################################
jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}
################################################################################
################################################################################
#Loading scripts
message("Loading complementary scripts")
#loading auxiliary code
source("~/SCRIPTS/INTERNAL_CODE/000_Prepare_dirs.R")
source("~/SCRIPTS/INTERNAL_CODE/000_wilcoxon_test_phyloseq_paired.R")
source("~/SCRIPTS/INTERNAL_CODE/001_Accumulation_Species_threshold.R")

#loading code to perform each of the comparison steps!
# source("~/SCRIPTS/002_regression_alpha_vars.R")
# source("~/SCRIPTS/003_Differentially_abundant_taxon.R")
source("~/SCRIPTS/COMPLEMENTARY_CODE/000_CREATING_FOLDERS.R")
source("~/SCRIPTS/COMPLEMENTARY_CODE/001_FILTERING.R")
source("~/SCRIPTS/COMPLEMENTARY_CODE/001_ggrare.R")
source("~/SCRIPTS/COMPLEMENTARY_CODE/002_ALPHA.R")
source("~/SCRIPTS/COMPLEMENTARY_CODE/002_BETA_PLOT.R")
source("~/SCRIPTS/COMPLEMENTARY_CODE/003_COMP_PLOT.R")
source("~/SCRIPTS/COMPLEMENTARY_CODE/004_summary_table.R")
source("~/SCRIPTS/COMPLEMENTARY_CODE/005_Wilcox_level.R")
source("~/SCRIPTS/COMPLEMENTARY_CODE/006_Jaccard_index.R")

################################################################################
#fix seed (Just in case)
set.seed(1000)
################################################################################
# Home folder
dir <- "~"
#Folder name
comparison_name <- "comparisons_ITS2_PACBIO_ILL_2"
#Create folder structure
dirs <- creat_dir_function(dir,comparison_name)
#Defining folder to use in the pipeline
out_dir <-dirs$out_dir
graph_dir <-dirs$graph_dir
dens_dir <- dirs$dens_dir
csv_dir <-dirs$csv_dir
RDS_dir <-dirs$RDS_dir

#how to name the results
markers_name <- c(
  "ITS_ILLUMINA",
  "ITS_PACBIO"
)

###########load files###########################################################
metadata <-
  as.data.frame(
    read.csv(
      "~/METADATA/metadata_sl-kazakhstan23.csv",header = T)
  )

#ILLUMINA
load("~/data/pipeline/kazakhstan-experimental/illumina/ITS2/phyloseq.Rdata")
ILL <- physeq;rm(physeq)
ILL_orig <- ILL
samp2 <- sample_names(ILL)
samp2 <- sub(pattern = "_ITS",replacement = "",x = samp2)
samp2 <- sub(pattern = "_",replacement = "",x = samp2)
sample_data(ILL)$SampleID <- samp2
sample_names(ILL) <- paste0(samp2,"_",markers_name[[1]])

#removing control
ILL <- prune_samples(sample_names(ILL) !="C1_KZ", ILL)
ILL <- phyloseq::prune_taxa(taxa_sums(ILL) > 0, ILL)
#filtering by fungi and read threshold
ILL <- filtering_reads_function(physeq = ILL,metadata,read_threshold=1000,name=markers_name[[1]])
ILL <- microViz::tax_agg(ps = ILL,"Species")

#PACBIO
load("~/data/pipeline/kazakhstan-experimental/pacbio/ITS9MUN/ITS2/phyloseq.Rdata")
PAC <- physeq;rm(physeq)
PAC_orig <- PAC
samp2 <- sample_names(PAC)
samp2 <- sub(pattern = "_ITS",replacement = "",x = samp2)
samp2 <- sub(pattern = "_",replacement = "",x = samp2)
sample_data(PAC)$SampleID <- samp2
sample_names(PAC) <- paste0(samp2,"_",markers_name[[2]])

#removing control
PAC <- prune_samples(sample_names(PAC) !="C1_KZ", PAC)
PAC <- phyloseq::prune_taxa(taxa_sums(PAC) > 0, PAC)
#filtering by fungi and read threshold
PAC <- filtering_reads_function(physeq = PAC,metadata,read_threshold=1000,name=markers_name[[2]])
PAC <- microViz::tax_agg(ps = PAC,"Species")

################################################################################
###############################################################################
#List original phyloseq objects
physeq_original <- list(ILL_orig,PAC_orig)
#List filtered phyloseq objects
physeq_list <- list(ILL,PAC)
################################################################################
###############################################################################
#obtaining shared samples
sample_list <- lapply(1:length(physeq_list),function(i){
  x <- sample_data(physeq_list[[i]])$SampleID
  x <- as.character(x)
  return(x) 
})

#getting common samples to compare
sam_to_test <- Reduce(intersect, sample_list)
################################################################################
###############################################################################
#subsetting to common samples
physeq_list <- lapply(1:length(physeq_list), function(i) {
  x <- prune_samples(sample_data(physeq_list[[i]])$SampleID %in% sam_to_test, physeq_list[[i]])
  x <- prune_taxa(taxa_sums(x) > 0, x)
  sample_data(x)$file <- markers_name[[i]]
  return(x)
})

################################################################################
################################################################################
if(!file.exists(paste0(RDS_dir,"/","002_joined_PSD.RDS"))){
  phy_joined <- merge_phyloseq(physeq_list[[1]],physeq_list[[2]])
  phy_joined <- prune_taxa(taxa_sums(phy_joined) > 0, phy_joined)
  
  saveRDS(phy_joined,paste0(RDS_dir,"/","002_joined_PSD.RDS"))
  
} else {
  phy_joined <- readRDS(paste0(RDS_dir, "/002_joined_PSD.RDS"))
  
}

################################################################################
#Creating first scenario
phy_joined_genus <-  microViz::tax_agg(ps = phy_joined,"Genus")
saveRDS(phy_joined_genus,paste0(RDS_dir,"/","GENUS_002_joined.RDS"))
################################################################################
#Creating second scenario
phy_joined_2 <- phy_joined

tax <- as.data.frame(tax_table(phy_joined_2))
#get unique species
UT <- unique(as.character(tax[,7]))
#getting taxa without taxonomic clarity and replacing by "?"
target_positions <- grep("(_sp$| Genus$| Family$| Order$| Class$| Phylum$| Domain$)", UT)
UT[target_positions] <- "?"
UT <- UT[which(UT!="?")]
tax[,7][!tax[,7] %in% UT] <- "?"
tax_table(phy_joined_2) <- as.matrix(tax)

phy_joined_sp <- phyloseq::subset_taxa(phy_joined_2, Species != "?")
saveRDS(phy_joined_sp,paste0(RDS_dir,"/","sp_002_joined.RDS"))
###############################################################################
################################################################################
###############################################################################
################################################################################
###############################################################################
################################################################################
###############################################################################
#rarefaction curves
#First scenario
ggrare_function(phy_joined_genus,name="Genus")
#Second scenario
ggrare_function(phy_joined_sp,name="Genus_sp_filtered")
###############################################################################
################################################################################
#ALPHA DIVESITY 
#First scenario
Alpha_function(physeq = phy_joined_genus,name="Genus")
#Second scenario
Alpha_function(physeq = phy_joined_sp,name="Genus_sp_filtered")
###############################################################################
################################################################################
#BETA DIVESITY 
#First scenario
physeq_clr2 <- microbiome::transform(phy_joined_genus, "compositional")
BETA_func(physeq_clr2,name="Genus")
#Second scenario
phy_joined_genus_sp <-  microViz::tax_agg(ps = phy_joined_sp,"Genus")
physeq_clr2_genus_sp <- microbiome::transform(phy_joined_genus_sp, "compositional")
BETA_func(physeq_clr2_genus_sp,name="Genus_sp_filtered")
###############################################################################
################################################################################
#Composition plots
#First scenario
comp_plot_function(physeq=phy_joined_genus,species_option=F,name="Genus")
#Second scenario
comp_plot_function(physeq=phy_joined_sp,species_option=T,name="Genus_sp_filtered")
################################################################################
################################################################################
#Summary tables
#First scenario
summary_table(markers = markers_name,
              physeq=phy_joined_genus,
              physeq_original=physeq_original,
              physeq_list=physeq_list,
              name="Genus")
#Second scenario
summary_table(markers = markers_name,
              physeq=phy_joined_sp,
              physeq_original=physeq_original,
              physeq_list=physeq_list,
              name="Genus_sp_filtered")
###############################################################################
################################################################################
#Differential abundant taxa
#Original
wilcox_per_level_func(levels=c("Family", "Genus", "Species"),
                      physeq=physeq_original,
                      markers_name=markers_name,
                      name="original",
                      dataset_name="original",
                      orig_status=T)

# x <- do.call(rbind,x)
# x2 <- data.frame(taxon = x$taxon,
#                  level = x$level,
#                  BP=NA,
#                  P = x$FDR)

#filtered
wilcox_per_level_func(levels=c("Family", "Genus", "Species"),
                      physeq=physeq_list,
                      markers_name=markers_name,
                      name="filtered",
                      dataset_name="Filtered_1000RC",orig_status=F)


AM_ps <- subset_samples(phy_joined_sp,file=="ITS_ILLUMINA")
AM_ps <- phyloseq::prune_taxa(taxa_sums(AM_ps) > 0, AM_ps)

ECM_ps <- subset_samples(phy_joined_sp,file=="ITS_PACBIO")
ECM_ps <- phyloseq::prune_taxa(taxa_sums(ECM_ps) > 0, ECM_ps)

ps_list <- list(AM_ps,ECM_ps)


x <-wilcox_per_level_func(levels=c("Family", "Genus", "Species"),
                      physeq=ps_list,
                      markers_name=markers_name,
                      name="test",
                      dataset_name="test_seq",orig_status=F)
for(i in 1:length(x)){
  x[[i]] <- x[[i]][order(x[[i]]$taxon,decreasing = F), ]
}
x2 <- do.call(rbind,x)
x2 <- data.frame(taxon = x2$taxon,
                 level = x2$level,
                 BP=NA,
                 P = x2$FDR)
################################################################################
###############################################################################
#Upset plots
#First scenario
Jaccard_upsetplot(levels=c("Family", "Genus"),
                  physeq=phy_joined_genus,
                  markers_name=markers_name,
                  name="Genus")
#Second scenario
Jaccard_upsetplot(levels=c("Family", "Genus", "Species"),
                  physeq=phy_joined_sp,
                  markers_name=markers_name,
                  name="Genus_sp_filtered")
################################################################################
###############################################################################