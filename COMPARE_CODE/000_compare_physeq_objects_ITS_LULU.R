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
  aplot
)
pacman::p_load_gh(
  "gauravsk/ranacapa",
  "gmteunisse/fantaxtic",
  "kasperskytte/ampvis2",
  "david-barnett/microViz@0.12.0",
  "thomasp85/patchwork"
)#,"twbattaglia/btools")
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
source("~/SCRIPTS/000_Prepare_dirs.R")
source("~/SCRIPTS/000_wilcoxon_test_phyloseq_paired.R")
# source("~/SCRIPTS/001_Accumulation_Species.R")
# source("~/SCRIPTS/002_regression_alpha_vars.R")
# source("~/SCRIPTS/003_Differentially_abundant_taxon.R")

################################################################################
set.seed(1000)
################################################################################
#Home folder
dir <- "~"
#MAKE DIR
out_comparison_dir <- paste0(dir,"/","comparisons_ITS2_LULU")
if(!dir.exists(out_comparison_dir)){
  dir.create(out_comparison_dir)
}

#how to name the results
markers <- c(
  "ITS2_1st_batch",
  "ITS2lulu"
)

#original_paths
paths <- c(
  "~/1stBatch/r_output/ITS2/ecm_physeq.Rdata",
  "~/data/pipeline/its2/101_sl-jgi_colombia25/ITS2lulu/ecm_physeq.Rdata"
)

compare_phyloseq <- function (dir, markers,out_comparison_dir,method){
  # method="counts"
  #MAKE graphic dir
  graph_dir <- paste0(out_comparison_dir,"/","graphics")
  if(!dir.exists(graph_dir)){
    dir.create(graph_dir)
  }
  
  #MAKE densities dir
  dens_dir <- paste0(graph_dir,"/","density")
  if(!dir.exists(dens_dir)){
    dir.create(dens_dir)
  }
  
  #MAKE csv dir
  csv_dir <- paste0(out_comparison_dir,"/","csv")
  if(!dir.exists(csv_dir)){
    dir.create(csv_dir)
  }
  #MAKE RDS dir
  RDS_dir <- paste0(out_comparison_dir,"/","RDS")
  if(!dir.exists(RDS_dir)){
    dir.create(RDS_dir)
  }
  
  # Validate transformation method argument
  if (!method %in% c("clr","compositional","counts")) {
    stop("Please use a valid transformation")
  }
  
  #Loading phyloseq objects
  physeq_list <- lapply(1:length(markers), function(i) {
    # i <- 1
    phy_i <- 
      readRDS(paste0(dir, "/", markers[[i]], "/outcomes/rds_dir/001_filtered_PSD.RDS"))
    phy_i <- microViz::tax_fix(phy_i, unknowns = c("?"))
    phy_i@sam_data$file <- markers[[i]]
    # #get taxonomic table to be modified
    tax_t <- tax_table(phy_i)
    #get unique species
    UT <- unique(as.character(tax_t[,7]))
    #getting taxa as Russula_sp or Russula Genus and replace by "?"
    target_positions <- grep("(_sp$| Genus$)", UT)
    UT[target_positions] <- "?"
    UT <- UT[which(UT!="?")]
    tax_t[,7][!tax_t[,7] %in% UT] <- "?"
    tax_table(phy_i) <- tax_t
    phy_i <- microViz::tax_fix(phy_i,unknowns = c("?"))
    # 
    return(phy_i)
  })

  #Aggregating rarefied object (physeq_rare)
  physeq_rare <- lapply(1:length(markers), function(i) {
    phy_i <- 
      readRDS(paste0(dir, "/", markers[[i]], "/outcomes/rds_dir/002_rarefied_PSD_4.RDS"))
    phy_i <- microViz::tax_fix(phy_i, unknowns = c("?"))
    phy_i@sam_data$file <- markers[[i]]
    #get taxonomic table to be modified
    tax_t <- tax_table(phy_i)
    #get unique species
    UT <- unique(as.character(tax_t[,7]))
    #getting taxa as Russula_sp or Russula Genus and replace by "?"
    target_positions <- grep("(_sp$| Genus$)", UT)
    UT[target_positions] <- "?"
    UT <- UT[which(UT!="?")]
    tax_t[,7][!tax_t[,7] %in% UT] <- "?"
    tax_table(phy_i) <- tax_t
    phy_i <- microViz::tax_fix(phy_i,unknowns = c("?"))
    return(phy_i)
  })
  
  #Loading original phyloseqs without rarefy
  physeq_original <- lapply(1:length(markers), function(i) {
    # i <- 1
    phy_o <- readRDS(paths[[i]])
    phy_o <- microViz::tax_fix(phy_o, unknowns = c("?"))
    #get taxonomic table to be modified
    tax_t <- tax_table(phy_o)
    #get unique species
    UT <- unique(as.character(tax_t[,7]))
    #getting taxa as Russula_sp or Russula Genus and replace by "?"
    target_positions <- grep("(_sp$| Genus$)", UT)
    UT[target_positions] <- "?"
    UT <- UT[which(UT!="?")]
    tax_t[,7][!tax_t[,7] %in% UT] <- "?"
    tax_table(phy_o) <- tax_t
    phy_o <- microViz::tax_fix(phy_o,unknowns = c("?"))
    return(phy_o)
  })
  ##############################################################################
  #SUMMARY TABLE
  message("creating summary, be patient")
  summary <- lapply(1:length(markers), function(i) {
    # i <- 1
    phy_o <- physeq_original[[i]]#readRDS(paths[[i]])
    #rarifed
    phy_i <-   physeq_rare[[i]] #readRDS(paste0(dir, "/", markers[[i]], "/outcomes/rds_dir/002_rarefied_PSD_4.RDS"))
    Alpha_f <-read.csv(paste0(dir, "/", markers[[i]], "/outcomes/csv_dir/002_Alpha_rarefaction.csv"))
    Beta_f <-read.csv(paste0(dir, "/", markers[[i]], "/outcomes/csv_dir/004_PERMANOVA_CLR.csv"))
    Thr_f <-read.csv(paste0(dir, "/", markers[[i]], "/outcomes/csv_dir/001_threshold_4.csv"))
    # Alpha_i <-read.csv(paste0(dir, "/", markers[[i]], "/outcomes/csv_dir/002_Alpha_rarefaction.csv"))
    # Beta_i <-read.csv(paste0(dir, "/", markers[[i]], "/outcomes/csv_dir/004_PERMANOVA_RARE_CLR.csv"))
    counts <- sample_sums(physeq_list[[i]])
    
    x <- data.frame(
      marker = markers[[i]],
      orig_nsamples = length(sample_sums(phy_o)),
      orig_read_counts = sum(sample_sums(phy_o)),
      orig_min = min(sample_sums(phy_o),na.rm=T),
      orig_max =  max(sample_sums(phy_o),na.rm=T),
      orig_min_sample = names(sample_sums(phy_o)[which.min(sample_sums(phy_o))]),
      orig_max_sample = names(sample_sums(phy_o)[which.max(sample_sums(phy_o))]),
      
      orig_families = length(unique(unlist(lapply(1:nrow(tax_table(physeq_original[[i]])),function(j){
        x <- paste(tax_table(physeq_original[[i]])[j,1:5],collapse = "-")
        return(x)
      })))),
      orig_genus = length(unique(unlist(lapply(1:nrow(tax_table(physeq_original[[i]])),function(j){
        x <- paste(tax_table(physeq_original[[i]])[j,1:6],collapse = "-")
        return(x)
      }))))
      ,
      orig_sp = length(unique(unlist(lapply(1:nrow(tax_table(physeq_original[[i]])),function(j){
        x <- paste(tax_table(physeq_original[[i]])[j,1:7],collapse = "-")
        return(x)
      })))),
      
      filtered_nsamples = length(counts),
      filtered_read_counts = sum(sample_sums(physeq_list[[i]])),
      filtered_min = min(counts,na.rm=T),
      filtered_max =  max(counts,na.rm=T),
      filtered_min_sample = names(counts[which.min(counts)]),
      filtered_max_sample = names(counts[which.max(counts)]),
      
      filtered_families = 
        length(unique(unlist(lapply(1:nrow(tax_table(physeq_list[[i]])),function(j){
          x <- paste(tax_table(physeq_list[[i]])[j,1:5],collapse = "-")
          return(x)
        })))),
      
      filtered_genus = 
        length(unique(unlist(lapply(1:nrow(tax_table(physeq_list[[i]])),function(j){
          x <- paste(tax_table(physeq_list[[i]])[j,1:6],collapse = "-")
          return(x)
        })))),
      
      filtered_sp = 
        length(unique(unlist(lapply(1:nrow(tax_table(physeq_list[[i]])),function(j){
          x <- paste(tax_table(physeq_list[[i]])[j,1:7],collapse = "-")
          return(x)
        })))),
      rare_threshold=Thr_f$threshold,
      rare_nsamples = length(sample_sums(phy_i)),
      rare_AM_nsample=sum(sample_data(phy_i)$Fungi_system=="AM"),
      rare_ECM_nsample=sum(sample_data(phy_i)$Fungi_system=="ECM"),
      
      rare_families = 
        length(unique(unlist(lapply(1:nrow(tax_table(phy_i)),function(j){
          x <- paste(tax_table(phy_i)[j,1:5],collapse = "-")
          return(x)
        })))),
      rare_genus = 
        length(unique(unlist(lapply(1:nrow(tax_table(phy_i)),function(j){
          x <- paste(tax_table(phy_i)[j,1:6],collapse = "-")
          return(x)
        })))),
      rare_sp = 
        length(unique(unlist(lapply(1:nrow(tax_table(phy_i)),function(j){
          x <- paste(tax_table(phy_i)[j,1:7],collapse = "-")
          return(x)
        })))),
      rare_obs_sp_av = mean(Alpha_f$Observed,na.rm=T),
      rare_obs_sp_sd = sd(Alpha_f$Observed,na.rm=T),
      rare_Shannon_av = mean(Alpha_f$Shannon,na.rm=T),
      rare_Shannon_sd = sd(Alpha_f$Shannon,na.rm=T),
      rare_InvSimpson_av = mean(Alpha_f$InvSimpson,na.rm=T),
      rare_InvSimpson_sd = sd(Alpha_f$InvSimpson,na.rm=T),
      rare_Pielou_av = mean(Alpha_f$pielou,na.rm=T),
      rare_Pielou_sd = sd(Alpha_f$pielou,na.rm=T),
      rare_Chao1_av = mean(Alpha_f$Chao1,na.rm = T),
      rare_Chao1_sd = sd(Alpha_f$Chao1,na.rm = T),
      # obs_sp_rare_av = mean(Alpha_i$Observed,na.rm=T),
      # obs_sp_sd_rare = sd(Alpha_i$Observed,na.rm=T),
      # Shannon_av_rare = mean(Alpha_i$Shannon,na.rm=T),
      # Shannon_sd_rare = sd(Alpha_i$Shannon,na.rm=T),
      # InvSimpson_av_rare = mean(Alpha_i$InvSimpson,na.rm=T),
      # InvSimpson_sd_rare = sd(Alpha_i$InvSimpson,na.rm=T),
      # Pielou_av_rare = mean(Alpha_i$pielou,na.rm=T),
      # Pielou_sd_rare = sd(Alpha_i$pielou,na.rm=T),
      rare_PERMANOVA_comp = Beta_f$Pr..F.[1]
      #PERMANOVA_rare = Beta_i$Pr..F.[1]
    )
    return(x)
  })
  
  summary <- do.call(rbind,summary)
  write.csv(summary,paste0(out_comparison_dir, "/comparison_summary.csv"),row.names = F)
  
  ################################################################################
  levels <-
    c("Class","Order","Family", "Genus", "Species")
  
  x_jacc <- lapply(1:length(levels), function(i) {
    
     i <- 5
    level <- levels[[i]]
    message(paste0("Level: ",level))
    ##############################################################################
    ##############################################################################
    ##############################################################################
    #aggregating taxonomic level using original
    phy_original <- lapply(1:length(physeq_original),function(j){
      phy1 <- physeq_original[[j]]
      if (level == "OTU") {
        phy1 <- phy1
      } else {
        phy1 <- microViz::tax_agg(ps = phy1, level)
        if(method!="counts"){
          phy1 <- microbiome::transform(phy1, method)          
        } else {
          phy1 <- phy1
        }
        
      }
      return(phy1)
    })
    
    #get taxonomic names
    phy_names_original <- lapply(1:length(phy_original), function(j) {
      # j <- 1
      tax1 <- tax_table(phy_original[[j]])
      otu1 <- as.character(unlist(row.names(tax1)))
      return(otu1)
    })
    #shared taxonomic names
    tax_to_test_original <- Reduce(intersect, phy_names_original)    
    ##############################################################################
    ##############################################################################
    message("Wilcoxon test using original Phyloseqs")
    if(length(phy_original)>2){
      stop("NO IMPLEMENTED FOR KRUSKAL WALLIS YET")
    } else {
      #get otu table
      x_list <- lapply(1:length(phy_original),function(l){
        x <- data.frame(otu_table(phy_original[[l]]))
        # x$dataset <- markers[[l]]
        return(x)
      })
      #get samples
      sample_list <- lapply(1:length(x_list),function(l){
        x <- colnames(x_list[[l]])
        x <- as.character(x)
        return(x) 
      })
      #get sample intersected
      sam_to_test <- Reduce(intersect, sample_list)    
      sam_to_test <- sam_to_test[sam_to_test!="dataset"]
      sam_to_test <- sam_to_test[!sam_to_test %in% c("NC_1","NC_2")]
      ##############################################################################
      #get wilcoxon paired tests
      wilc_test <- wilcoxon_function_microb(tax_to_test_original,x_list,sam_to_test)
      wilc_test <- do.call(rbind,wilc_test)
      #Adjusting FDR
      wilc_test$FDR <- p.adjust(wilc_test$p.value,method = "fdr")
      #Ordering ascending
      wilc_test <- wilc_test[order(wilc_test$FDR,decreasing = F),]
      wilc_test$level <- level
      wilc_test <- wilc_test[which(wilc_test$FDR<0.05),]
      
    }
    if(method=="clr"){
      wilc_test$method <- "CLR"
      wilc_test$dataset <- paste0("ORIGINAL-",wilc_test$method,"-",wilc_test$level)
      write.csv(wilc_test,paste0(csv_dir, "/WILCOXON", "_","ORIGINAL","CLR","_", level, ".csv"),row.names=F)
    } else if(method=="compositional"){
      wilc_test$method <- "RA"
      wilc_test$dataset <- paste0("ORIGINAL-",wilc_test$method,"-",wilc_test$level)
      write.csv(wilc_test,paste0(csv_dir, "/WILCOXON", "_","ORIGINAL","RA","_", level, ".csv"),row.names=F)
    } else if(method=="counts"){
      wilc_test$method <- "COUNT"
      wilc_test$dataset <- paste0("ORIGINAL-",wilc_test$method,"-",wilc_test$level)
      write.csv(wilc_test,paste0(csv_dir, "/WILCOXON", "_","ORIGINAL","COUNT","_", level, ".csv"),row.names=F)
    }
    ##############################################################################
    #Obtaining data from two datasets to create density plots
    plot_density_function(dens_dir,wilc_test,markers,x_list,sam_to_test,method,level,dataset="original")
    ##############################################################################
    #Aggregating taxa using rarefied objects
    message("Aggregating rarefied data to selected level")
    physeq_list2 <- lapply(1:length(physeq_list), function(j) {
      # print(j)
      phy1 <- physeq_list[[j]]
      if (level == "OTU") {
        phy1 <- phy1
      } else {
        phy1 <- microViz::tax_agg(ps = phy1, level)
      }
      if(nrow(tax_table(phy1))>1){
        
        #microbiome::transform(phy1, "compositional")  # Convert to relative abundance
        if(method=="clr"){
          message(paste0("using clr method for Wilcoxon test for:",level))
          phy1 <-         microbiome::transform(phy1, "clr")  # Convert to relative abundance
        } else if(method=="compositional"){
          message(paste0("using relative abundances method for Wilcoxon test for:",level))
          phy1 <-         microbiome::transform(phy1, "compositional")
        } else if(is.null(method)){
          message(paste0("using countsfor Wilcoxon test for:",level))
          phy1 <-         phy1
        }
        
      } else {
        message("ONLY ONE TAXON")
        return(NULL)
      }
      return(phy1)
    })
    
    #obtaining names for Jaccard approach
    phy_names <- lapply(1:length(physeq_list2), function(j) {
      # j <- 1
      tax1 <- tax_table(physeq_list2[[j]])
      otu1 <- unlist(row.names(tax1))
      otu1 <- as.character(otu1)
      return(otu1)
    })
    
    names(phy_names) <- markers
    
    
    
    #getting common to compare
    common <- Reduce(intersect, phy_names)
    #unique taxa
    phy_names_unique <- lapply(1:length(physeq_list2), function(j) {
      # j <- 1
      tax1 <- tax_table(physeq_list2[[j]])
      otu1 <- unlist(row.names(tax1))
      otu1 <- as.character(otu1)
      otu1 <- otu1[!otu1 %in% common]
      if(length(otu1)>0){
        otu1 <- data.frame(taxon=otu1,dataset=markers[[j]])        
      } else {
        otu1 <- data.frame(taxon=NA,dataset=markers[[j]])
      }
      
      return(otu1)
    })
    
    ##############################################################################
    ##############################################################################
    ##############################################################################
    message("Plotting taxonomic similarities")
    
    #Jaccard UpsetPlot
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
      filename = paste0(graph_dir,"/ITS2_VENN", "_", level, ".pdf"),
      plot = venn_plot_combined,
      width = 15,
      height = 10,
      units = "in",
      dpi = 600,
      device = cairo_pdf
    )
    ##############################################################################
    ##############################################################################
    ##############################################################################
    #Jaccard
    message("Calculating Jaccard index")
    
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
    write.csv(jacc_mat,paste0(csv_dir, "/JACCARD", "_", level, ".csv"))
    
    
    
    ##############################################################################
    ##############################################################################
    ##############################################################################
    #WILCOXON FILTERED!
    
    #aggregating taxonomic level using original
    phy_filtered <- lapply(1:length(physeq_list),function(j){
      phy1 <- physeq_list[[j]]
      if (level == "OTU") {
        phy1 <- phy1
      } else {
        phy1 <- microViz::tax_agg(ps = phy1, level)
        if(method!="counts"){
          phy1 <- microbiome::transform(phy1, method)          
        } else {
          phy1 <- phy1
        }
        
      }
      return(phy1)
    })
    
    
    
    #obtaining names for Jaccard approach
    phy_names <- lapply(1:length(phy_filtered), function(j) {
      # j <- 1
      tax1 <- tax_table(phy_filtered[[j]])
      otu1 <- unlist(row.names(tax1))
      otu1 <- as.character(otu1)
      return(otu1)
    })
    
    names(phy_names) <- markers
    
    
    #getting common to compare
    common <- Reduce(intersect, phy_names)
    #unique taxa
    phy_names_unique_filtered <- lapply(1:length(phy_filtered), function(j) {
      tax1 <- tax_table(phy_filtered[[j]])
      otu1 <- unlist(row.names(tax1))
      otu1 <- as.character(otu1)
      otu1 <- otu1[!otu1 %in% common]
      if(length(otu1)>0){
        otu1 <- data.frame(taxon=otu1,dataset=markers[[j]])        
      } else {
        otu1 <- data.frame(taxon=NA,dataset=markers[[j]])
      }
      return(otu1)
    })
    
    
    
    message("Wilcoxon test for filtered phyloseqs")
    
    #test differences in composition
    tax_to_test <- Reduce(intersect, phy_names)
    
    if(length(physeq_list)>2){
      stop("NO IMPLEMENTED FOR KRUSKAL WALLIS YET")
    } else {
      #get otu table
      x_list <- lapply(1:length(phy_filtered),function(l){
        x <- data.frame(otu_table(phy_filtered[[l]]))
        # x$dataset <- markers[[l]]
        return(x)
      })
      #get samples
      sample_list <- lapply(1:length(x_list),function(l){
        x <- colnames(x_list[[l]])
        x <- as.character(x)
        return(x) 
      })
      #get taxa to do
      tax_to_test <- Reduce(intersect, phy_names)
      #get sample intersected
      sam_to_test <- Reduce(intersect, sample_list)
      sam_to_test <- sam_to_test[sam_to_test!="dataset"]
      
      
      #get wilcoxon paired tests
      wilc_test <- wilcoxon_function_microb(tax_to_test,x_list,sam_to_test)
      wilc_test <- do.call(rbind,wilc_test)
      #Adjusting FDR
      wilc_test$FDR <- p.adjust(wilc_test$p.value,method = "fdr")
      #Ordering ascending
      wilc_test <- wilc_test[order(wilc_test$FDR,decreasing = F),]
      wilc_test$level <- level
      wilc_test <- wilc_test[which(wilc_test$FDR<0.05),]
    }
    if(method=="clr"){
      wilc_test$method <- "CLR"
      wilc_test$dataset <- paste0("FILTERED-",wilc_test$method,"-",wilc_test$level)
      write.csv(wilc_test,paste0(csv_dir, "/WILCOXON", "_","FILTERED_","CLR","_", level, ".csv"),row.names=F)
    } else if(method=="compositional"){
      wilc_test$method <- "RA"
      wilc_test$dataset <- paste0("FILTERED-",wilc_test$method,"-",wilc_test$level)
      write.csv(wilc_test,paste0(csv_dir, "/WILCOXON", "_","FILTERED_","RA","_", level, ".csv"),row.names=F)
    } else if(method=="counts"){
      wilc_test$method <- "COUNT"
      wilc_test$dataset <- paste0("FILTERED-",wilc_test$method,"-",wilc_test$level)
      write.csv(wilc_test,paste0(csv_dir, "/WILCOXON", "_","FILTERED_","COUNT","_", level, ".csv"),row.names=F)
    }
    
    #Obtaining data from two datasets to create density plots
    plot_density_function(dens_dir,wilc_test,markers,x_list,sam_to_test,method,level,dataset="filtered")
    ##############################################################################
    ##############################################################################
    ##############################################################################
    #WILCOXON RAREFIED!
    message("Wilcoxon test for rarefied phyloseqs")
    
    
    physeq_rare_2 <- lapply(1:length(physeq_rare),function(j){
      phy1 <- physeq_rare[[j]]
      if (level == "OTU") {
        phy1 <- phy1
      } else {
        phy1 <- microViz::tax_agg(ps = phy1, level)
        if(method!="counts"){
          phy1 <- microbiome::transform(phy1, method)          
        } else {
          phy1 <- phy1
        }
        
      }
      return(phy1)
    })
    
    #obtaining names for Jaccard approach
    phy_names <- lapply(1:length(physeq_rare_2), function(j) {
      # j <- 1
      tax1 <- tax_table(physeq_rare_2[[j]])
      otu1 <- unlist(row.names(tax1))
      otu1 <- as.character(otu1)
      return(otu1)
    })
    names(phy_names) <- markers
    
    
    #getting common to compare
    common <- Reduce(intersect, phy_names)
    #unique taxa
    phy_names_unique_rarified <- lapply(1:length(physeq_rare_2), function(j) {
      tax1 <- tax_table(physeq_rare_2[[j]])
      otu1 <- unlist(row.names(tax1))
      otu1 <- as.character(otu1)
      otu1 <- otu1[!otu1 %in% common]
      if(length(otu1)>0){
        otu1 <- data.frame(taxon=otu1,dataset=markers[[j]])        
      } else {
        otu1 <- data.frame(taxon=NA,dataset=markers[[j]])
      }
      return(otu1)
    })
    
    #test differences in composition
    tax_to_test <- Reduce(intersect, phy_names)
    
    
    if(length(physeq_rare_2)>2){
      stop("NO IMPLEMENTED FOR KRUSKAL WALLIS YET")
    } else {
      #get otu table
      x_list <- lapply(1:length(physeq_rare_2),function(l){
        x <- data.frame(otu_table(physeq_rare_2[[l]]))
        # x$dataset <- markers[[l]]
        return(x)
      })
      #get samples
      sample_list <- lapply(1:length(x_list),function(l){
        x <- colnames(x_list[[l]])
        x <- as.character(x)
        return(x) 
      })
      #get taxa to do
      tax_to_test <- Reduce(intersect, phy_names)
      #get sample intersected
      sam_to_test <- Reduce(intersect, sample_list)
      sam_to_test <- sam_to_test[sam_to_test!="dataset"]
      
      
      #get wilcoxon paired tests
      wilc_test <- wilcoxon_function_microb(tax_to_test_original = tax_to_test,x_list = x_list,
                                            sam_to_test =sam_to_test)
      wilc_test <- do.call(rbind,wilc_test)
      #Adjusting FDR
      wilc_test$FDR <- p.adjust(wilc_test$p.value,method = "fdr")
      #Ordering ascending
      wilc_test <- wilc_test[order(wilc_test$FDR,decreasing = F),]
      wilc_test$level <- level
      wilc_test <- wilc_test[which(wilc_test$FDR<0.05),]
    }
    
    if(method=="clr"){
      wilc_test$method <- "CLR"
      wilc_test$dataset <- paste0("RAREFIED-",wilc_test$method,"-",wilc_test$level)
      write.csv(wilc_test,paste0(csv_dir, "/WILCOXON", "_","RARIFIED_","CLR","_", level, ".csv"),row.names=F)
    } else if(method=="compositional"){
      wilc_test$method <- "RA"
      wilc_test$dataset <- paste0("RAREFIED-",wilc_test$method,"-",wilc_test$level)
      write.csv(wilc_test,paste0(csv_dir, "/WILCOXON", "_","RARIFIED_","RA","_", level, ".csv"),row.names=F)
    } else if(method=="counts"){
      wilc_test$method <- "COUNT"
      wilc_test$dataset <- paste0("RAREFIED-",wilc_test$method,"-",wilc_test$level)
      write.csv(wilc_test,paste0(csv_dir, "/WILCOXON", "_","RARIFIED_","COUNT","_", level, ".csv"),row.names=F)
    }
    
    ############################################################################
    message("detecting taxa to show density differences")
    #taxa to use for density plots
    taxa_to_test_density <- wilc_test$taxon
    
    #Obtaining data from two datasets to create density plots
    plot_density_function(dens_dir,wilc_test,markers,x_list,sam_to_test,method,level,dataset="rarified")
    ############################################################################
    
    
    #obtain unique taxa per dataset and object and saving a CSV file
    if(level %in% levels[4:5]){
      message("Saving unique taxa in a CSV file")
      phy_names_unique <- do.call(rbind,phy_names_unique)
      phy_names_unique_filtered <- do.call(rbind,phy_names_unique_filtered)
      phy_names_unique_rarified <- do.call(rbind,phy_names_unique_rarified)
      phy_names_unique$object <- "No filtered"
      phy_names_unique_filtered$object <- "Filtered"
      phy_names_unique_rarified$object <- "Rarefied"
      
      list_unique <- rbind(phy_names_unique,phy_names_unique_filtered)
      list_unique <- rbind(list_unique,phy_names_unique_rarified)
      list_unique <- list_unique[complete.cases(list_unique),]
      write.csv(list_unique,paste0(csv_dir, "/UNIQUE_TAXON", "_","_", level, ".csv"),row.names=F)
      ###########################################################################
    }
    return(physeq_list2)
  })
  
  message("Saving image")
  
  save.image(paste0(RDS_dir,"/","Comparisons_lotus_",method,".RData"))
  
}

# compare_phyloseq(dir, markers,out_comparison_dir,method="clr")
compare_phyloseq(dir, markers,out_comparison_dir,method="counts")
# compare_phyloseq(dir, markers,out_comparison_dir,method="counts")
