#' @title Microbiome Comparative Analysis
#' @description Script for comparing multiple microbial marker datasets across taxonomic levels using phyloseq objects. Includes diversity analysis, Jaccard similarity, and Venn diagram generation.
#' @author Chrystian Sosa
#' @keywords microbiome, phyloseq, diversity, jaccard, venn, visualization
#' @import pacman phyloseq microbiome ggplot2 ggpubr xlsx readr dplyr compositions ggrepel iNEXT gridExtra ggsci UpSetR ggVennDiagram ggplotify aplot
#' @importFrom microViz tax_fix tax_agg
#' @importFrom microbiome transform
#' @importFrom ggVennDiagram ggVennDiagram
#' @importFrom ggplotify as.ggplot
#' @importFrom ggplot2 ggsave ggtitle theme element_text element_blank
#' @importFrom stats sd
#' @importFrom utils read.csv write.csv
#' @export

#' @examples
#' # Not intended to be run as a standalone example.
#' # Source this file in a package or project where marker folders and paths exist.
#' source("microbiome_comparative_analysis.R")

#' @seealso \code{phyloseq}, \code{microViz}, \code{ggVennDiagram}

# Load necessary packages using pacman
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
library(pacman)

# Load or install packages from CRAN and GitHub
pacman::p_load(
  phyloseq, microbiome, ggplot2, ggpubr, xlsx, readr, dplyr,
  compositions, ggrepel, iNEXT, gridExtra, ggsci, UpSetR, 
  ggVennDiagram, ggplotify, aplot
)
pacman::p_load_gh(
  "gauravsk/ranacapa",
  "gmteunisse/fantaxtic",
  "kasperskytte/ampvis2",
  "david-barnett/microViz@0.12.0",
  "thomasp85/patchwork"
)

#' @title Jaccard index calculation
#' @description Computes the Jaccard index between two sets.
#' @param a First set (character vector)
#' @param b Second set (character vector)
#' @return Numeric Jaccard similarity
#' @export
jaccard <- function(a, b) {
  intersection <- length(intersect(a, b))
  union <- length(a) + length(b) - intersection
  return(intersection / union)
}

################################################################################
#Loading scripts
message("Loading complementary scripts")
source("~/SCRIPTS/000_Prepare_dirs.R")
source("~/SCRIPTS/001_Accumulation_Species.R")
source("~/SCRIPTS/002_regression_alpha_vars.R")
source("~/SCRIPTS/003_Differentially_abundant_taxon.R")


################################################################################
#Home folder
dir <- "~"
#how to name the results
markers <- c(
  "SSU_dada2_results",
  "SSU_eukaryome",
  "SSU_VSEARCH"
)

#original_paths
paths <- c(
  "~/00101_20250411JC1K0Y/r_output/SSU_dada2/amf_physeq.Rdata",
  "~/00101_20250411JC1K0Y/r_output/SSU_eukaryome/amf_physeq.Rdata",
  "~/00101_20250411JC1K0Y/r_output/SSU_vsearch//amf_physeq.Rdata"
)
#Loading phyloseq objects
physeq_list <- lapply(1:length(markers), function(i) {
  phy_i <- 
  readRDS(paste0(dir, "/", markers[[i]], "/outcomes/rds_dir/001_filtered_PSD.RDS"))
  phy_i <- microViz::tax_fix(phy_i, unknowns = c("?"))
  phy_i@sam_data$file <- markers[[i]]
  return(phy_i)
})

#Create a Summary file with the original Phyloseq
summary <- lapply(1:length(markers), function(i) {
#Loading rarefied Phyloseq, Alpha, Beta, Threshold, counts
phy_o <- readRDS(paths[[i]])
phy_i <-   readRDS(paste0(dir, "/", markers[[i]], "/outcomes/rds_dir/002_rarefied_PSD_4.RDS"))
Alpha_f <-read.csv(paste0(dir, "/", markers[[i]], "/outcomes/csv_dir/002_Alpha.csv"))
Beta_f <-read.csv(paste0(dir, "/", markers[[i]], "/outcomes/csv_dir/004_PERMANOVA_RA.csv"))

Thr_i <-read.csv(paste0(dir, "/", markers[[i]], "/outcomes/csv_dir/001_threshold_4.csv"))
Alpha_i <-read.csv(paste0(dir, "/", markers[[i]], "/outcomes/csv_dir/002_Alpha_rarefaction.csv"))
Beta_i <-read.csv(paste0(dir, "/", markers[[i]], "/outcomes/csv_dir/004_PERMANOVA_RARE_CLR.csv"))

counts <- sample_sums(physeq_list[[i]])
x <- data.frame(
  marker = markers[[i]],
  nsamples = length(sample_sums(phy_o)),
  read_counts = sum(sample_sums(phy_o)),
  min = min(sample_sums(phy_o),na.rm=T),
  max =  max(sample_sums(phy_o),na.rm=T),
  min_sample = names(sample_sums(phy_o)[which.min(sample_sums(phy_o))]),
  max_sample = names(sample_sums(phy_o)[which.max(sample_sums(phy_o))]),
  
  nsamples_filtered = length(counts),
  read_counts_filtered = sum(sample_sums(physeq_list[[i]])),
  min_filtered = min(counts,na.rm=T),
  max_filtered=  max(counts,na.rm=T),
  min_sample_filtered = names(counts[which.min(counts)]),
  max_sample_filtered = names(counts[which.max(counts)]),
  
  rare_threshold=Thr_i$threshold,
  nsamples_rare = length(sample_sums(phy_i)),
  obs_sp_av = mean(Alpha_f$Observed,na.rm=T),
  obs_sp_sd = sd(Alpha_f$Observed,na.rm=T),
  Shannon_av = mean(Alpha_f$Shannon,na.rm=T),
  Shannon_sd = sd(Alpha_f$Shannon,na.rm=T),
  InvSimpson_av = mean(Alpha_f$InvSimpson,na.rm=T),
  InvSimpson_sd = sd(Alpha_f$InvSimpson,na.rm=T),
  Pielou_av = mean(Alpha_f$pielou,na.rm=T),
  Pielou_sd = sd(Alpha_f$pielou,na.rm=T),
  
  obs_sp_rare_av = mean(Alpha_i$Observed,na.rm=T),
  obs_sp_sd_rare = sd(Alpha_i$Observed,na.rm=T),
  Shannon_av_rare = mean(Alpha_i$Shannon,na.rm=T),
  Shannon_sd_rare = sd(Alpha_i$Shannon,na.rm=T),
  InvSimpson_av_rare = mean(Alpha_i$InvSimpson,na.rm=T),
  InvSimpson_sd_rare = sd(Alpha_i$InvSimpson,na.rm=T),
  Pielou_av_rare = mean(Alpha_i$pielou,na.rm=T),
  Pielou_sd_rare = sd(Alpha_i$pielou,na.rm=T),
  PERMANOVA_comp = Beta_f$Pr..F.[1],
  PERMANOVA_rare = Beta_i$Pr..F.[1]
)
return(x)
})

summary <- do.call(rbind,summary)
write.csv(summary,paste0(dir, "/comparisons/SSU_summary.csv"),row.names = F)

################################################################################
#Define levels
levels <-
  c("Order","Family", "Genus", "Species", "OTU")

#Create Jaccard
x_jacc <- lapply(1:length(levels), function(i) {
  
  level <- levels[[i]]
  
  physeq_list2 <- lapply(1:length(physeq_list), function(j) {
    #j <- 3
    # print(j)
    phy1 <- physeq_list[[j]]
    if (level == "OTU") {
      phy1 <- phy1
    } else {
      phy1 <- microViz::tax_agg(ps = phy1, level)
    }
    if(nrow(tax_table(phy1))>1){
      phy1 <-
        microbiome::transform(phy1, "compositional")  # Convert to relative abundance
    } else {
      message("ONLY ONE TAXON")
      return(NULL)
    }
    return(phy1)
  })
    
  phy_names <- lapply(1:length(physeq_list2), function(j) {
    # j <- 1
    tax1 <- tax_table(physeq_list2[[j]])
    otu1 <- unlist(row.names(tax1))
    return(otu1)
      })
  
  names(phy_names) <- markers
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
      axis.title.y = element_blank(),  # ya no lo necesitas si es redundante
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
    ggtitle(paste0("Overlap of ", level, " taxonomic level"))+
    theme(plot.title = element_text(size = 24, face = "bold"))
  
  # Save to high-resolution PDF
  ggsave(
    filename = paste0(dir, "/comparisons/SSU_VENN", "_", level, ".pdf"),
    plot = venn_plot_combined,
    width = 15,
    height = 10,
    units = "in",
    dpi = 600,
    device = cairo_pdf
  )
  
  jacc_mat <- matrix(ncol=length(phy_names),nrow=length(phy_names))
  jacc_mat <- as.data.frame(jacc_mat)
  colnames(jacc_mat) <- markers
  row.names(jacc_mat) <- markers
  for(k in 1:length(phy_names)){
    for(l in 2:length(phy_names)){
      if(l>k){
        jacc_mat[k,l] <- 
          jaccard(phy_names[[k]],phy_names[[l]])
      }
    }
  }
  diag(jacc_mat) <- 1
  write.csv(jacc_mat,paste0(dir, "/comparisons/SSU_JACCARD", "_", level, ".csv"))
  return(physeq_list2)
})
