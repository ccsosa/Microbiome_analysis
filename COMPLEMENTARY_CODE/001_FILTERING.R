library(phyloseq)
library(microViz)
library(dplyr)

filtering_reads_function <- function(physeq,metadata,read_threshold,name){
  message("Filtering by using fungi domain")
  #only filtering at Fungi domain
  physeq <- phyloseq::subset_taxa(physeq, Domain == "Fungi")
  #removing possible otu without counts
  physeq <- phyloseq::prune_taxa(taxa_sums(physeq) > 0, physeq)
  #creating a column sampleID to join to
  # sample_data(physeq)$SampleID <- sample_names(physeq)
  message("Loading complementary metadata to include into the analysis")
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

 
  #get samples with zero counts and remove them
  physeq4 <-
    phyloseq::prune_samples(samples = names(phyloseq::sample_sums(physeq)[sample_sums(physeq) >
                                                                          0]),
                            x = physeq)
  
  #Remove samples with less than 1000 reads
  physeq4 <-
    phyloseq::prune_samples(samples = names(phyloseq::sample_sums(physeq4)[sample_sums(physeq4) >=
                                                                           read_threshold]),
                          x = physeq4)

  #Prune taxa with zero counts
  physeq4 <- phyloseq::prune_taxa(taxa_sums(physeq4) > 0, physeq4)
  sample_data(physeq4)$orig_id <- sample_names(physeq4)
  samp_data <- as.data.frame(as.matrix(sample_data(physeq4)))
  
  row.names(samp_data) <- sample_data(physeq4)$orig_id
  
  samp_data <-
    dplyr::left_join(x = samp_data,
                   y = metadata,
                   by = c("SampleID" = "SampleID"))
  row.names(samp_data) <- samp_data$orig_id
  
  sample_data(physeq4) <- samp_data
################################################################################
  #Saving Phyloseq object
  message("Phyloseq object filtered obtained...Saving to RDS subfolder")
  saveRDS(physeq4, paste0(RDS_dir, "/",name,"_001_filtered_PSD.RDS"))
  return(physeq4)
}
