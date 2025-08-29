library(phyloseq)
library(microViz)
library(dplyr)

wilcox_per_level_func <- function(levels,physeq,name,markers_name,dataset_name,orig_status){
  x_jacc <- lapply(1:length(levels), function(i) {
    # i <- 1
    level <- levels[[i]]
    #aggregating taxonomic level using original
    phy_aggregate <- lapply(1:length(physeq),function(j){
      phy1 <- physeq[[j]]
      phy1 <- tax_fix(phy1,unknowns = c("?"))
      phy1 <- microViz::tax_agg(ps = phy1, level)
      
      return(phy1)
    })
    
    #get taxonomic names
    phy_names <- lapply(1:length(phy_aggregate), function(j) {
      # j <- 1
      tax1 <- tax_table(phy_aggregate[[j]])
      otu1 <- as.character(unlist(row.names(tax1)))
      return(otu1)
    })
    #shared taxonomic names
    tax_to_test <- Reduce(intersect, phy_names)    
    ##############################################################################
    ##############################################################################
    message(paste0("Wilcoxon test for: ",level))
    if(length(phy_aggregate)>2){
      stop("NO IMPLEMENTED FOR KRUSKAL WALLIS YET")
    } else {
      #get otu table
      x_list <- lapply(1:length(phy_aggregate),function(l){
        x <- data.frame(otu_table(phy_aggregate[[l]]))
        # x$dataset <- markers[[l]]
        return(x)
      })
      #get samples
      sample_list <- lapply(1:length(x_list),function(l){
        x <- colnames(x_list[[l]])
        x <- as.character(x)
        return(x) 
      })
      
      if(isFALSE(orig_status)){
        sample_list[[1]] <- sub(pattern = "_ITS_ILLUMINA",replacement = "",x = sample_list[[1]])
        sample_list[[2]] <- sub(pattern = "_ITS_PACBIO",replacement = "",x = sample_list[[2]])
        colnames(x_list[[1]]) <- sample_list[[1]]
        colnames(x_list[[2]]) <- sample_list[[2]]
      } else if(isTRUE(orig_status)) {
        message("original phyloseq object used!")
        sample_list[[2]] <- sub(pattern = "_ITS",replacement = "",x = colnames(x_list[[2]]))
        sample_list[[2]] <- sub(pattern = "_",replacement = "",x = sample_list[[2]])
        colnames(x_list[[2]]) <- sample_list[[2]]
      }
      #fixing pacbio sample names
      

      # 
      #get sample intersected
      sam_to_test <- Reduce(intersect, sample_list)    
      sam_to_test <- sam_to_test[sam_to_test!="dataset"]
      ##############################################################################
      #get wilcoxon paired tests
      wilc_test <- wilcoxon_function_microb(tax_to_test_original = tax_to_test,
                                            x_list = x_list,
                                            sam_to_test = sam_to_test)
      wilc_test <- do.call(rbind,wilc_test)
      #Adjusting FDR
      wilc_test$FDR <- p.adjust(wilc_test$p.value,method = "fdr")
      #Ordering ascending
      # wilc_test <- wilc_test[order(wilc_test$FDR,decreasing = F),]
      wilc_test$level <- level
      wilc_test_total <- wilc_test
      wilc_test <- wilc_test[which(wilc_test$FDR<0.05),]
      
    }
    if(nrow(wilc_test)>0){
      wilc_test$method <- "COUNT"
      wilc_test$dataset <- paste0(dataset_name,wilc_test$method,"-",wilc_test$level)
      write.csv(wilc_test,paste0(csv_dir, "/",name,"_",dataset_name,"_WILCOXON", "_","ORIGINAL","COUNT","_", level, ".csv"),row.names=F)
    }
    ##############################################################################
    #Obtaining data from two datasets to create density plots
    if(nrow(wilc_test)>0){
      plot_density_function(dens_dir = dens_dir,
                            wilc_test = wilc_test,
                            markers = markers_name,
                            x_list = x_list,
                            sam_to_test=sam_to_test,
                            method="counts",
                            level=level,
                            dataset=dataset_name)      
    } else {
      message("no Diff taxa for original Phyloseq object")
    }
    return(wilc_test_total)
 })
  message("DONE!")
  return(x_jacc)
}