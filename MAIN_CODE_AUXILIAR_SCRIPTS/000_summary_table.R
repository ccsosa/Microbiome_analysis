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
  grDevices
)
pacman::p_load_gh(
  "gauravsk/ranacapa",
  "gmteunisse/fantaxtic",
  "kasperskytte/ampvis2",
  "david-barnett/microViz@0.12.0"
)#,"twbattaglia/btools")


summary_table_function <- function(physeq4){
  #reads
  reads <- sample_sums(physeq4)
  ################################################################################
  #Getting taxonomical info and its levels
  tax_table_data <- data.frame(tax_table(physeq4))
  levels <- colnames(tax_table_data)
  #unique taxa per level 
  unique_tax <-  lapply(1:length(levels),function(i){
    #Getting info for levels
    x <- tax_table_data[,levels[[i]]]
    # Count the number of unique units in levels
    x <- length(unique(x))
    return(x)
  })
  ################################################################################
  #Subsetting by AM
  AM_ps <- subset_samples(physeq4,Habit=="AM")
  AM_ps <- phyloseq::prune_taxa(taxa_sums(AM_ps) > 0, AM_ps)
  reads_AM <- sample_sums(AM_ps)
  #Getting taxonomical info and its levels
  tax_table_data <- data.frame(tax_table(AM_ps))
  levels <- colnames(tax_table_data)
  #unique taxa per level 
  unique_tax_AM <-  lapply(1:length(levels),function(i){
    #Getting info for levels
    x <- tax_table_data[,levels[[i]]]
    # Count the number of unique units in levels
    x <- length(unique(x))
    return(x)
  })
  ################################################################################
  #Subsetting by ECM
  ECM_ps <- subset_samples(physeq4,Habit=="ECM")
  ECM_ps <- phyloseq::prune_taxa(taxa_sums(ECM_ps) > 0, ECM_ps)
  reads_ECM <- sample_sums(ECM_ps)
  
  tax_table_data <- data.frame(tax_table(ECM_ps))
  levels <- colnames(tax_table_data)
  #unique taxa per level 
  unique_tax_ECM <-  lapply(1:length(levels),function(i){
    #Getting info for levels
    x <- tax_table_data[,levels[[i]]]
    # Count the number of unique units in levels
    x <- length(unique(x))
    return(x)
  })
  
  ################################################################################
  #GENERAL
    df <- data.frame(
    dataset = "total",
    total = sum(reads),
    min = min(reads),
    max = max(reads),
    mean = mean(reads),
    sd = sd(reads),
    median = median(reads),
    mad = mad(reads),
    Domain = unique_tax[[1]],
    Phylum = unique_tax[[2]],
    Class = unique_tax[[3]],
    Order = unique_tax[[4]],
    Family = unique_tax[[5]],
    Genus = unique_tax[[6]],
    Species = unique_tax[[7]],
    OTU = length(unique(row.names(tax_table(physeq4))))
    )
    #AM
  df_AM <- data.frame(
    dataset = "AM",
    total = sum(reads_AM),
    min = min(reads_AM),
    max = max(reads_AM),
    mean = mean(reads_AM),
    sd = sd(reads_AM),
    median = median(reads_AM),
    mad = mad(reads_AM),
    Domain = unique_tax_AM[[1]],
    Phylum = unique_tax_AM[[2]],
    Class = unique_tax_AM[[3]],
    Order = unique_tax_AM[[4]],
    Family = unique_tax_AM[[5]],
    Genus = unique_tax_AM[[6]],
    Species = unique_tax_AM[[7]],
    OTU = length(unique(row.names(tax_table(AM_ps))))
  )
    #ECM
    df_ECM <- data.frame(
    dataset = "ECM",
    total = sum(reads_ECM),
    min = min(reads_ECM),
    max = max(reads_ECM),
    mean = mean(reads_ECM),
    sd = sd(reads_ECM),
    median = median(reads_ECM),
    mad = mad(reads_ECM),
    Domain = unique_tax_ECM[[1]],
    Phylum = unique_tax_ECM[[2]],
    Class = unique_tax_ECM[[3]],
    Order = unique_tax_ECM[[4]],
    Family = unique_tax_ECM[[5]],
    Genus = unique_tax_ECM[[6]],
    Species = unique_tax_ECM[[7]],
    OTU = length(unique(row.names(tax_table(ECM_ps))))
  )
  
    df_final <- rbind(df,df_AM)
    df_final <- rbind(df_final,df_ECM)
  write.csv(
    df_final,
    paste0(csv_dir, "/", "Summary_table.csv"),
    na = "",
    row.names = F,
    quote = T
  )
  return(df_final)
}