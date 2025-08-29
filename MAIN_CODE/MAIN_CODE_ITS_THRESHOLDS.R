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
  FactoMineR,
  factoextra,
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
source("~/SCRIPTS/001_Accumulation_Species_threshold.R")
source("~/SCRIPTS/002_regression_alpha_vars.R")
source("~/SCRIPTS/003_Differentially_abundant_taxon.R")


#Home folder
dir <- "~"
#how to name the results
marker="ITS2_results_threshold_test"
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
physeq <- readRDS("~/00101_20250411JC1K0Y/r_output/ITS2/ecm_physeq.Rdata")
ss <- data.frame(sample_data(physeq))
write.csv(ss,paste0(csv_dir,"/","ORIGINAL_SAMPLE_DATA.csv"),row.names = T)
#Checking again if there are strange Phylum, Class, Order or families
#sample_data(physeq)
################################################################################
thresholds <- c(500,1000,2000,5000,10000)
################################################################################
  #0 Loading extra metadata
  message("Loading complementary metadata to include into the analysis")
  metadata <-
    as.data.frame(
      readxl::read_xlsx(
        "~/METADATA/Database samples_JGI_20250612.xlsx",
        sheet = "Soil",
        col_names = TRUE,
        skip = 1,
        na = "NA"
      )
    )
  metadata$Code <-
    sub(pattern = "[-]",
        replacement = "_",
        x = metadata$Code)
  metadata$Code <-
    sub(pattern = "[-]",
        replacement = "_",
        x = metadata$Code)
  metadata$...16 <- NULL
  metadata$...17 <- NULL
  metadata$...18 <- NULL
  metadata$`AM/ECM`[which(metadata$`AM/ECM` == "oak")] <- "ECM"
  metadata <-
    metadata[, c("Code",
                 "Plot",
                 "Individual",
                 "Replica",
                 "AM/ECM",
                 "Height",
                 "DBH",
                 "Altitude")]
  colnames(metadata)[5] <- "Fungi_system"
  metadata$Height <- as.numeric(metadata$Height)
  
  ################################################################################
  #Fixing species as Genus + species
  
  message("Fixing species as Genus + species for eukaryome")
  
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
  
  physeq@tax_table@.Data[,7]  <- as.character(tax$Species)
  rm(tax)
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
  
  #4. getting negative controls to prune
  negatives <-
    row.names(samp_data[which(samp_data$Sample_or_Control == "Control")])
  asv_neg <- asv_bin_df[negatives, ]
  # asv_neg[3,] <- colSums(asv_neg,na.rm = T)
  # #OTUs associated to negative samples
  # colnames(asv_neg[,asv_neg[3,]>0])
  # neg <- tax_table(physeq)[colnames(asv_neg[,asv_neg[3,]>0]),]
  #remove sample data
  #rm(samp_data)
  
  #5. Prune negative controls
  physeq2 <- phyloseq::prune_samples(samples = row.names(samp_data)[!row.names(samp_data) %in% negatives],
                                     x = physeq)
  #6. Prune taxa with no counts
  physeq3 <- phyloseq::prune_taxa(taxa_sums(physeq2) > 0, physeq2)
  #get unique Phylum
  #unique(tax_table(physeq3)[,"Phylum"])
  
  # filterPhyla = c(
  #   "Annelida",
  #   "Apicomplexa",
  #   "Arthropoda",
  #   "Ascomycota",
  #   "Tardigrada",
  #   "Vertebrata",
  #   "Rotifera",
  #   "Mollusca",
  #   "Cnidaria",
  #   "Cercozoa",
  #   "Ciliophora",
  #   "Schizoplasmodiida",
  #   "Echinodermata",
  #   "Nematozoa",
  #   "Diatomea",
  #   "Platyhelminthes",
  #   "Ciliophora",
  #   "Porifera",
  #   "Kathablepharidae",
  #   "Phragmoplastophyta",
  #   "Protalveolata",
  #   "Nematoda",
  #   "Centrohelida",
  #   "Apusomonadidae",
  #   "Gracilipodida",
  #   "Rigifilida",
  #   "Katablepharidophyta",
  #   "Lycopodiophyta",
  #   "Eukaryota_phy_Incertae_sedis",
  #   "Fungi_phy_Incertae_sedis",
  #   "?"
  # )
  ################################################################################
  # message("Filtering strange Phyla")
  # 
  # #7. Filter to remove strange Phylum
  # physeq3 <- phyloseq::subset_taxa(physeq3,!Phylum %in% filterPhyla)
  
  #8. Filter to at least get samples for 10% of samples
  #physeq4 <- phyloseq::filter_taxa(physeq3, function(x) sum(x > 2) > (0.1*length(x)), TRUE)
  #get samples with zero counts and remove them
  #Remove samples with less than 200 reads
  
  #Adding extra information to sample data
  #using metadata associated to samples table
  samp_data <- as.data.frame(as.matrix(sample_data(physeq3)))
  samp_data$id <- row.names(samp_data)
  samp_data <-
    dplyr::left_join(x = samp_data,
                     y = metadata,
                     by = c("id" = "Code"))
  row.names(samp_data) <- samp_data$id
  sample_data(physeq3) <- samp_data
  
  ##############################################################################
  ##############################################################################
  ##############################################################################
  #Creating different phyloseq objects per read count thresholds
  phy_i <- list()
  for(i in 1:length(thresholds)){
    physeq5 <-
      phyloseq::prune_samples(samples = names(phyloseq::sample_sums(physeq3)[sample_sums(physeq3) >=
                                                    thresholds[[i]]]),x = physeq3)    
    ################################################################################
    #Saving Phyloseq object
    message("Phyloseq object filtered obtained...Saving to RDS subfolder")
    saveRDS(physeq5, paste0(RDS_dir, "/","001_filtered_PSD_",thresholds[[i]],".RDS"))
    phy_i[[i]] <- physeq5
  }
  
  ##############################################################################
  ##############################################################################
  ##############################################################################
  ##############################################################################
  #Rarefying using Cmin after thresholds
  
  phy_rare <-lapply(1:length(thresholds),function(i){
    # i <- 1
  phy_i_to_rare <- phy_i[[i]]
  #Cmin is the sample with minimum read counts
  Cmin=min(sample_sums(phy_i_to_rare))
  #rarefying using SRS and plotting
  ps_rare <- accumulation_curve_threshold_function(physeq4 = phy_i_to_rare,
                                                   option="SRS",
                                                   threshold_label=thresholds[[i]],
                                                   Cmin=Cmin)
    return(ps_rare)
  })
  
  ##############################################################################
  #Rarefying using Cmin after thresholds
  phy_rare2 <-lapply(1:length(thresholds),function(i){
    phy_i_to_rare <- phy_i[[i]]
    #Cmin is twice the number of mmininum read counts
    Cmin=min(sample_sums(phy_i_to_rare))*2
    ps_rare <- accumulation_curve_threshold_function(physeq4 = phy_i_to_rare,
                                                     option="SRS",
                                                     threshold_label=paste0(thresholds[[i]],"_","2CMin"),
                                                     Cmin=Cmin)
    return(ps_rare)
  })
  
  phy_rare2
  
################################################################################
################################################################################
################################################################################  
################################################################################  
#Calculating Alpha values  
  Alpha_value <- lapply(1:length(phy_rare),function(i){
    #richness
    Alpha_rare <- phyloseq::estimate_richness(phy_rare[[i]])
    #Pielou
    Alpha_rare$pielou <- evenness(phy_rare[[i]], index = "pielou")[, 1]
    #Chao1
    calc <- ChaoRichness(x = otu_table(phy_rare[[i]]), datatype = "abundance", conf = 0.95)
    Alpha_rare <- Alpha_rare[,-c(1:5)]
    Alpha_rare <- cbind(Alpha_rare,calc)
    colnames(Alpha_rare) [7:8] <- c("Chao1","Chao1_s.e.")
    Alpha_rare$id <- row.names(Alpha_rare)
    x <- data.frame(id = as.character(row.names(sample_data(phy_rare[[i]]))),
               Fungi_system = as.character(sample_data(phy_rare[[i]])$Fungi_system)
               )  
    #Joining in one file
    Alpha_rare <- left_join(x = Alpha_rare,y = x,by = c("id"="id"))
    row.names(Alpha_rare) <- Alpha_rare$id
    #Adding read counts
    Alpha_rare$min <- unique(sample_sums(phy_rare[[i]]))
    #Adding thresholds
    Alpha_rare$threshold <- thresholds[[i]]
    #Adding CMin to the dataset
    Alpha_rare$dataset <- paste0(thresholds[[i]],"_","CMin")
    return(Alpha_rare)
  })
  
  Alpha_value <- do.call(rbind,Alpha_value)
################################################################################  
  Alpha_value2 <- lapply(1:length(phy_rare2),function(i){
    #Richness
    Alpha_rare <- phyloseq::estimate_richness(phy_rare2[[i]])
    #Pielou
    Alpha_rare$pielou <- evenness(phy_rare2[[i]], index = "pielou")[, 1]
    #Chao1
    calc <- ChaoRichness(x = otu_table(phy_rare2[[i]]), datatype = "abundance", conf = 0.95)
    Alpha_rare <- Alpha_rare[,-c(1:5)]
    Alpha_rare <- cbind(Alpha_rare,calc)
    colnames(Alpha_rare) [7:8] <- c("Chao1","Chao1_s.e.")
    Alpha_rare$id <- row.names(Alpha_rare)
    x <- data.frame(id = as.character(row.names(sample_data(phy_rare2[[i]]))),
                    Fungi_system = as.character(sample_data(phy_rare2[[i]])$Fungi_system)
    )  
     #Joining in one file
    Alpha_rare <- left_join(x = Alpha_rare,y = x,by = c("id"="id"))
    row.names(Alpha_rare) <- Alpha_rare$id
   #adding Cmin
    Alpha_rare$min <- unique(sample_sums(phy_rare2[[i]]))
    #Adding thresold
    Alpha_rare$threshold <- thresholds[[i]]
    Alpha_rare$dataset <- paste0(thresholds[[i]],"_","2xCMin")
    return(Alpha_rare)
  })
  
  Alpha_value2 <- do.call(rbind,Alpha_value2)
################################################################################  
  #Saving Alpha in a unique file
  Alpha_v <- rbind(Alpha_value,Alpha_value2)
  
  message("saving Alpha values")
  write.csv(Alpha_value,paste0(csv_dir,"/","000_Alpha_table_thr.csv"))
################################################################################    
################################################################################
################################################################################
#Summary file for Alpha
n = tapply(Alpha_v$Observed,Alpha_v$dataset,length)
#Summary tables  
summ_table <- data.frame(
  n =n,
  #Observed species mean and sd
  Observed_mean = tapply(Alpha_v$Observed,Alpha_v$dataset,mean),
  Observed_sd = tapply(Alpha_v$Observed,Alpha_v$dataset,sd),
  #Shannon mean and sd
  Shannon_mean = tapply(Alpha_v$Shannon,Alpha_v$dataset,mean),
  Shannon_sd = tapply(Alpha_v$Shannon,Alpha_v$dataset,sd),
  #ENS mean and sd
  InvSimpson_mean = tapply(Alpha_v$InvSimpson,Alpha_v$dataset,mean),
  InvSimpson_sd = tapply(Alpha_v$InvSimpson,Alpha_v$dataset,sd),
  #Chao1 mean
  Chao1_mean = tapply(Alpha_v$Chao1,Alpha_v$dataset,mean),
  #Mininum read counts
  min= tapply(Alpha_v$min,Alpha_v$dataset,min),
  #dataset names
  dataset = names(n),
  #Read count thresholds
  threshold = tapply(Alpha_v$threshold,Alpha_v$dataset,min)
  )
  
#Adding factor level for plots
  Alpha_v$dataset <- factor(Alpha_v$dataset,levels =  
              c("500_CMin","500_2xCMin",
                "1000_CMin","1000_2xCMin",
                "2000_CMin","2000_2xCMin",
                "5000_CMin","5000_2xCMin",
                "10000_CMin","10000_2xCMin")
  )
  #     c("500_CMin", "1000_CMin","2000_CMin","5000_CMin","10000_CMin",
  # "500_2xCMin","1000_2xCMin","2000_2xCMin","5000_2xCMin","10000_2xCMin"))
################################################################################
################################################################################
################################################################################
#PLOTS

#Defining comparisons for Wilcoxon rank sum tests
  my_comparisons <- list(
    c("500_CMin", "500_2xCMin"),
    c("1000_CMin", "1000_2xCMin"),
    c("2000_CMin", "2000_2xCMin"),
    c("5000_CMin", "5000_2xCMin"),
    c("10000_CMin", "10000_2xCMin")
  )

###Observed species
p1 <-   ggpubr::ggboxplot(Alpha_v,x = "dataset",y = "Observed",fill="Fungi_system",xlab = "") +
    stat_compare_means(comparisons = my_comparisons,
                       method = "wilcox.test",
                       label = "p.format",
                       size = 5) +  # Tamaño del texto de p-valor    ggsci::scale_color_jco() +   # jco palette for lines/points
    ggsci::scale_fill_jco() +    # jco palette for boxplots
    theme(
      # legend.position = "None",
      strip.text = element_text(size = 16),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 14),
      axis.title.x = element_text(size = 18),
      axis.title.y = element_text(size = 18),
      axis.text.x = element_text(size = 16),
      axis.text.y = element_text(size = 22),
      plot.title = element_text(size = 24, face = "bold", hjust = 0.5)
    )+ 
  rotate_x_text()

###Chao1
p2 <-   ggpubr::ggboxplot(Alpha_v,x = "dataset",y = "Chao1",fill="Fungi_system",xlab = "") +
  stat_compare_means(comparisons = my_comparisons,
                     method = "wilcox.test",
                     label = "p.format",
                     size = 5) +  # Tamaño del texto de p-valor    ggsci::scale_color_jco() +   # jco palette for lines/points
  ggsci::scale_fill_jco() +    # jco palette for boxplots
  theme(
    # legend.position = "None",
    strip.text = element_text(size = 16),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 22),
    plot.title = element_text(size = 24, face = "bold", hjust = 0.5)
  )+ 
  rotate_x_text()

###Pielou
p3 <-   ggpubr::ggboxplot(Alpha_v,x = "dataset",y = "pielou",fill="Fungi_system",xlab = "") +
  stat_compare_means(comparisons = my_comparisons,
                     method = "wilcox.test",
                     label = "p.format",
                     size = 5) +  # Tamaño del texto de p-valor    ggsci::scale_color_jco() +   # jco palette for lines/points
  ggsci::scale_fill_jco() +    # jco palette for boxplots
  theme(
    # legend.position = "None",
    strip.text = element_text(size = 16),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 22),
    plot.title = element_text(size = 24, face = "bold", hjust = 0.5)
  ) + rotate_x_text()

###ENS
p4 <-   ggpubr::ggboxplot(Alpha_v,x = "dataset",y = "InvSimpson",fill="Fungi_system",xlab = "") +
  stat_compare_means(comparisons = my_comparisons,
                     method = "wilcox.test",
                     label = "p.format",
                     size = 5) +  # Tamaño del texto de p-valor    ggsci::scale_color_jco() +   # jco palette for lines/points
  ggsci::scale_fill_jco() +    # jco palette for boxplots
  theme(
    # legend.position = "None",
    strip.text = element_text(size = 16),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 22),
    plot.title = element_text(size = 24, face = "bold", hjust = 0.5)
  )+ 
  rotate_x_text()
#Join in one grid arranges
p_Alpha <- gridExtra::grid.arrange(p1,p2,p3,p4)

message("Alpha indicators...Saved in graphics subfolder")
ggsave(
  paste0(graph_dir, "/", "002_Alpha_boxplots_datasets.pdf"),
  p_Alpha,
  width = 20,
  height = 25,
  units = "in",
  dpi = 600
)
################################################################################
p1 <-   ggpubr::ggboxplot(Alpha_v,x = "dataset",y = "Observed",fill="dataset",xlab = "") +
  stat_compare_means(comparisons = my_comparisons,
                     method = "wilcox.test",
                     label = "p.format",
                     size = 5) +  # Tamaño del texto de p-valor    ggsci::scale_color_jco() +   # jco palette for lines/points
  ggsci::scale_fill_jco() +    # jco palette for boxplots
  theme(
    # legend.position = "None",
    strip.text = element_text(size = 16),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 22),
    plot.title = element_text(size = 24, face = "bold", hjust = 0.5)
  )+ 
  rotate_x_text()


p2 <-   ggpubr::ggboxplot(Alpha_v,x = "dataset",y = "Chao1",fill="dataset",xlab = "") +
  stat_compare_means(comparisons = my_comparisons,
                     method = "wilcox.test",
                     label = "p.format",
                     size = 5) +  # Tamaño del texto de p-valor    ggsci::scale_color_jco() +   # jco palette for lines/points
  ggsci::scale_fill_jco() +    # jco palette for boxplots
  theme(
    # legend.position = "None",
    strip.text = element_text(size = 16),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 22),
    plot.title = element_text(size = 24, face = "bold", hjust = 0.5)
  )+ 
  rotate_x_text()

p3 <-   ggpubr::ggboxplot(Alpha_v,x = "dataset",y = "pielou",fill="dataset",xlab = "") +
  stat_compare_means(comparisons = my_comparisons,
                     method = "wilcox.test",
                     label = "p.format",
                     size = 5) +  # Tamaño del texto de p-valor    ggsci::scale_color_jco() +   # jco palette for lines/points
  ggsci::scale_fill_jco() +    # jco palette for boxplots
  theme(
    # legend.position = "None",
    strip.text = element_text(size = 16),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 22),
    plot.title = element_text(size = 24, face = "bold", hjust = 0.5)
  ) + rotate_x_text()

p4 <-   ggpubr::ggboxplot(Alpha_v,x = "dataset",y = "InvSimpson",fill="dataset",xlab = "") +
  stat_compare_means(comparisons = my_comparisons,
                     method = "wilcox.test",
                     label = "p.format",
                     size = 5) +  # Tamaño del texto de p-valor    ggsci::scale_color_jco() +   # jco palette for lines/points
  ggsci::scale_fill_jco() +    # jco palette for boxplots
  theme(
    # legend.position = "None",
    strip.text = element_text(size = 16),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 22),
    plot.title = element_text(size = 24, face = "bold", hjust = 0.5)
  )+ 
  rotate_x_text()

p_Alpha <- gridExtra::grid.arrange(p1,p2,p3,p4)

message("Alpha indicators...Saved in graphics subfolder")
ggsave(
  paste0(graph_dir, "/", "002_Alpha_boxplots_datasets_NOF.pdf"),
  p_Alpha,
  width = 20,
  height = 25,
  units = "in",
  dpi = 600
)

# p5 <-   ggpubr::ggboxplot(Alpha_v,x = "dataset",y = "min",fill="Fungi_system") +
#   stat_compare_means(comparisons = my_comparisons,
#                      method = "wilcox.test",
#                      label = "p.format",
#                      size = 5) +  # Tamaño del texto de p-valor    ggsci::scale_color_jco() +   # jco palette for lines/points
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

  ggsave(
    paste0(graph_dir, "/", "002_Alpha_boxplots_.pdf"),
    p_Alpha,
    width = 20,
    height = 25,
    units = "in",
    dpi = 600
  )

  
  ##############################################################################
  
    
  summ_table$dataset <- factor(summ_table$dataset,levels =  
                                c("500_CMin","500_2xCMin",
                                  "1000_CMin","1000_2xCMin",
                                  "2000_CMin","2000_2xCMin",
                                  "5000_CMin","5000_2xCMin",
                                  "10000_CMin","10000_2xCMin")
)

write.csv(summ_table,paste0(csv_dir,"/","000_summary_table.csv"))

Alpha_v$rare_threshold <- NA
Alpha_v$rare_threshold[Alpha_v$dataset %in% c("500_CMin",
                                               "1000_CMin",
                                               "2000_CMin",
                                               "5000_CMin",
                                               "10000_CMin")] <- "CMin"
  
Alpha_v$rare_threshold[Alpha_v$dataset %in% c("500_2xCMin",
                                              "1000_2xCMin",
                                              "2000_2xCMin",
                                              "5000_2xCMin",
                                              "10000_2xCMin")] <- "2xCMin"
  
Alpha_PCA <- PCA(Alpha_v[,c(1,3,5,6,7)],scale.unit = T,graph = F)

# fviz_pca_biplot(Alpha_PCA,
#              habillage = Alpha_v$dataset,
#              repel = TRUE,
#              label = "var", # color by groups
#              # palette = c("#00AFBB", "#E7B800", "#FC4E07"),
#              addEllipses = TRUE, ellipse.type = "convex",col.var = "black") 
# 
# 
# Extraer coordenadas de individuos
ind_df <- as.data.frame(Alpha_PCA$ind$coord)
ind_df$dataset <- factor(as.character(Alpha_v$threshold),c("500","1000","2000","5000","10000"))
ind_df$rare_threshold <- Alpha_v$rare_threshold

# 3. Crear convex hulls por grupo
hull_df <- ind_df %>%
  group_by(dataset, rare_threshold) %>%
  slice(chull(Dim.1, Dim.2))
# 2. Extraer coordenadas de variables (para flechas)
var_df <- as.data.frame(Alpha_PCA$var$coord)
var_df$varname <- rownames(var_df)

# Escalado opcional de vectores de variables
scale_factor <- 5
var_df_scaled <- var_df %>%
  mutate(Dim.1 = Dim.1 * scale_factor,
         Dim.2 = Dim.2 * scale_factor)

# Graficar con ggrepel
ggplot(ind_df, aes(x = Dim.1, y = Dim.2, color = dataset)) +
  # convex hulls
  geom_polygon(data = hull_df, aes(fill = dataset, group = dataset), 
               alpha = 0.2, color = NA) +
  # puntos
  geom_point(size = 2, alpha = 0.8) +
  # flechas de variables
  geom_segment(data = var_df_scaled,
               aes(x = 0, y = 0, xend = Dim.1, yend = Dim.2),
               arrow = arrow(length = unit(0.3, "cm")),
               inherit.aes = FALSE, color = "black") +
  # etiquetas sin solapamiento
  geom_text_repel(data = var_df_scaled,
                  aes(x = Dim.1, y = Dim.2, label = varname),
                  color = "black",
                  size = 4,
                  max.overlaps = Inf,
                  box.padding = 0.3,
                  point.padding = 0.2,
                  segment.color = "gray30",
                  inherit.aes = FALSE) +
  facet_grid(. ~ rare_threshold) +
  theme_minimal() +
  labs(title = "PCA Biplot with Convex Hulls and Non-overlapping Variable Arrows",
       x = "PC1", y = "PC2") +
  theme(legend.position = "bottom")

save.image("~/ITS2_results_threshold_test/RDAta_test.RData")
#load("~/ITS2_results_threshold_test/RDAta_test.RData")
