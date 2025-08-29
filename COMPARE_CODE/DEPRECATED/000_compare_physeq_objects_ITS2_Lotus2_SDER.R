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
# source("~/SCRIPTS/001_Accumulation_Species.R")
# source("~/SCRIPTS/002_regression_alpha_vars.R")
# source("~/SCRIPTS/003_Differentially_abundant_taxon.R")

################################################################################
set.seed(1000)
################################################################################
#Home folder
dir <- "~"

#MAKE DIR
out_comparison_dir <- paste0(dir,"/","comp_lotus_SEP_DEP")
if(!dir.exists(out_comparison_dir)){
  dir.create(out_comparison_dir)
}

#how to name the results
markers <- c(
  #"ITS2_chimechecktogether.lotus2tax",
  "ITS2_chimecheckseparate.lotus2tax",
  "ITS2_chimechecktogether.derepmin.lotus2tax"
)

#original_paths
paths <- c(
  "/data/pipeline/its2/101_sl-jgi_colombia25/chimecheckseparate.lotus2tax/ecm_physeq.Rdata",
  "/data/pipeline/its2/101_sl-jgi_colombia25/chimechecktogether.derepmin.lotus2tax/ecm_physeq.Rdata"
)


compare_phyloseq <- function (dir, markers,out_comparison_dir,method){
  # method="compositional"
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
    phy_i <- 
      readRDS(paste0(dir, "/", markers[[i]], "/outcomes/rds_dir/001_filtered_PSD.RDS"))
    phy_i <- microViz::tax_fix(phy_i, unknowns = c("?"))
    phy_i@sam_data$file <- markers[[i]]
    return(phy_i)
  })
  
  #Loading original phyloseqs without rarefy
  physeq_original <- lapply(1:length(markers), function(i) {
    phy_o <- readRDS(paths[[i]])
    phy_o <- microViz::tax_fix(phy_o, unknowns = c("?"))
    return(phy_o)
  })
  
  message("creating summary, be patient")
  summary <- lapply(1:length(markers), function(i) {
     # i <- 1
    phy_o <- readRDS(paths[[i]])
    phy_i <-   readRDS(paste0(dir, "/", markers[[i]], "/outcomes/rds_dir/002_rarefied_PSD_4.RDS"))
    Alpha_f <-read.csv(paste0(dir, "/", markers[[i]], "/outcomes/csv_dir/002_Alpha_rarefaction.csv"))
    Beta_f <-read.csv(paste0(dir, "/", markers[[i]], "/outcomes/csv_dir/004_PERMANOVA_RARE_CLR.csv"))
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
  
  # i <- 3
  x_jacc <- lapply(1:length(levels), function(i) {
    
    # i <- 5
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
        phy1 <- microbiome::transform(phy1, "compositional")
      }
      return(phy1)
    })
    
    #get taxonomic names
    phy_names_original <- lapply(1:length(phy_original), function(j) {
      # j <- 1
      tax1 <- tax_table(phy_original[[j]])
      otu1 <- unlist(row.names(tax1))
      return(otu1)
    })
    #shared taonomic names
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
      sample_list <- lapply(length(x_list),function(l){
        x <- colnames(x_list[[l]])
        return(x) 
      })
      #get taxa to do
      
      #get sample intersected
      sam_to_test <- Reduce(intersect, sample_list)
      sam_to_test <- sam_to_test[sam_to_test!="dataset"]
      sam_to_test <- sam_to_test[!sam_to_test %in% c("NC_1","NC_2")]
      
      #get wilcoxon paired tests
      wilc_test <- lapply(1:length(tax_to_test_original),function(m){
        # print(m)
        # m <- 2
        x <- lapply(1:length(x_list),function(n){
          x <- data.frame(value = as.numeric(x_list[[n]][tax_to_test_original[[m]],sam_to_test]),
                          marker = markers[[n]])
          return(x)
        })
        #x <- do.call(rbind, x)
        #Calculating wilcoxon paired test
        x2 <- data.frame(taxon=tax_to_test_original[[m]],
                         p.value=wilcox.test(x[[1]][,1],x[[2]][,1], exact = T,paired=T)$p.value,
                         ks=ks.test(as.numeric(x[[1]][,1]),as.numeric(x[[2]][,1]))$p.value,
                         median_1 = median(x[[1]][,1],na.rm = T),
                         mad_1 = mad(x[[1]][,1],na.rm = T),
                         mean_1 = mean(x[[1]][,1],na.rm = T),
                         sd_1 = sd(x[[1]][,1],na.rm = T),
                         median_2 = median(x[[2]][,1],na.rm = T),
                         mad_2 = mad(x[[2]][,1],na.rm = T),
                         mean_2 = mean(x[[2]][,1],na.rm = T),
                         sd_2 = sd(x[[2]][,1],na.rm = T)
        )
        # colnames(x2) <- c("taxon","p.value")
        return(x2)
      })
      wilc_test <- do.call(rbind,wilc_test)
      #Adjusting FDR
      wilc_test$FDR <- p.adjust(wilc_test$p.value,method = "fdr")
      #Ordering ascending
      wilc_test <- wilc_test[order(wilc_test$FDR,decreasing = F),]
      wilc_test$level <- level
    }
    
    
    wilc_test$method <- "RA"
    wilc_test$dataset <- paste0("ORIGINAL","-",wilc_test$level)
    wilc_test <- wilc_test[order(wilc_test$FDR,decreasing = F),]
    wilc_test <- wilc_test[which(wilc_test$FDR<0.05),]
    write.csv(wilc_test,paste0(csv_dir, "/WILCOXON", "_","ORIGINAL","_", level, ".csv"),row.names=F)
    
    ##############################################################################
    message("detecting taxa to show density differences")
    #taxa to use for density plots
    taxa_to_test_density <- wilc_test$taxon
    
    #Obtaining data from two datasetss to create density plots
    x <- lapply(1:length(taxa_to_test_density),function(m){
      # m <- 1
      y <- lapply(1:length(markers),function(n){
        x <- data.frame(value = as.numeric(x_list[[n]][taxa_to_test_density[[m]],sam_to_test]),
                        marker = markers[[n]])
        return(x)
      })
      #Getting data in a unique file
      y <- do.call(rbind,y)
      #Plotting in a density plot     
      x_dens <- ggdensity(y, x = "value",
                          add = "mean", rug = TRUE,
                          color = "marker", fill = "marker",
                          palette = c("#00AFBB", "#E7B800"),repel = T,
                          xlab = "Relative abundance",ylab = "Density",
                          title = paste0(level,": ",taxa_to_test_density[[m]])
      ) +
        
        labs(fill = NULL, color = NULL)
      
      # Save to high-resolution PDF
      ggsave(
        filename = paste0(dens_dir,"/",level,"_",taxa_to_test_density[[m]],"_original_RA", ".pdf"),
        plot = x_dens,
        width = 15,
        height = 10,
        units = "in",
        dpi = 600,
        device = cairo_pdf
      )
    })
    
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
      return(otu1)
    })
    
    names(phy_names) <- markers
    
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
    #WILCOXON
    
    message("Wilcoxon test for rarefied phyloseqs")
    
    #test differences in composition
    tax_to_test <- Reduce(intersect, phy_names)
    
    if(length(physeq_list2)>2){
      stop("NO IMPLEMENTED FOR KRUSKAL WALLIS YET")
    } else {
      #get otu table
      x_list <- lapply(1:length(physeq_list2),function(l){
        x <- data.frame(otu_table(physeq_list2[[l]]))
        # x$dataset <- markers[[l]]
        return(x)
      })
      #get samples
      sample_list <- lapply(length(x_list),function(l){
        x <- colnames(x_list[[l]])
        return(x) 
      })
      #get taxa to do
      tax_to_test <- Reduce(intersect, phy_names)
      #get sample intersected
      sam_to_test <- Reduce(intersect, sample_list)
      sam_to_test <- sam_to_test[sam_to_test!="dataset"]
      
      #get wilcoxon paired tests
      wilc_test <- lapply(1:length(tax_to_test),function(m){
        # m <- 1
        x <- lapply(1:length(x_list),function(n){
          x <- data.frame(value = as.numeric(x_list[[n]][tax_to_test[[m]],sam_to_test]),
                          marker = markers[[n]])
          return(x)
        })
        
        #x <- do.call(rbind, x)
        #Calculating wilcoxon paired test
        x2 <- data.frame(taxon=tax_to_test_original[[m]],
                         p.value=wilcox.test(x[[1]][,1],x[[2]][,1], exact = T,paired=T)$p.value,
                         ks=ks.test(as.numeric(x[[1]][,1]),as.numeric(x[[2]][,1]))$p.value,
                         median_1 = median(x[[1]][,1],na.rm = T),
                         mad_1 = mad(x[[1]][,1],na.rm = T),
                         mean_1 = mean(x[[1]][,1],na.rm = T),
                         sd_1 = sd(x[[1]][,1],na.rm = T),
                         median_2 = median(x[[2]][,1],na.rm = T),
                         mad_2 = mad(x[[2]][,1],na.rm = T),
                         mean_2 = mean(x[[2]][,1],na.rm = T),
                         sd_2 = sd(x[[2]][,1],na.rm = T)
        )
        # colnames(x2) <- c("taxon","p.value")
        return(x2)
      })
      wilc_test <- do.call(rbind,wilc_test)
      #Adjusting FDR
      wilc_test$FDR <- p.adjust(wilc_test$p.value,method = "fdr")
      #Ordering ascending
      wilc_test <- wilc_test[order(wilc_test$FDR,decreasing = F),]
      wilc_test$level <- level
    }
    
    if(method=="clr"){
      wilc_test$method <- "CLR"
      wilc_test$dataset <- paste0(wilc_test$method,"-",wilc_test$level)
      write.csv(wilc_test,paste0(csv_dir, "/WILCOXON", "_","CLR","_", level, ".csv"),row.names=F)
    } else if(method=="compositional"){
      wilc_test$method <- "RA"
      wilc_test$dataset <- paste0(wilc_test$method,"-",wilc_test$level)
      write.csv(wilc_test,paste0(csv_dir, "/WILCOXON", "_","RA","_", level, ".csv"),row.names=F)
    } else if(method=="counts"){
      wilc_test$method <- "COUNT"
      wilc_test$dataset <- paste0(wilc_test$method,"-",wilc_test$level)
      write.csv(wilc_test,paste0(csv_dir, "/WILCOXON", "_","COUNT","_", level, ".csv"),row.names=F)
    }
    
    ############################################################################
    message("detecting taxa to show density differences")
    #taxa to use for density plots
    taxa_to_test_density <- wilc_test$taxon
    
    #Obtaining data from two datasetss to create density plots
    x <- lapply(1:length(taxa_to_test_density),function(m){
      # m <- 1
      y <- lapply(1:length(markers),function(n){
        x <- data.frame(value = as.numeric(x_list[[n]][taxa_to_test_density[[m]],sam_to_test]),
                        marker = markers[[n]])
        return(x)
      })
      #Getting data in a unique file
      y <- do.call(rbind,y)
      #Plotting in a density plot     
      x_dens <- ggdensity(y, x = "value",
                          add = "mean", rug = TRUE,
                          color = "marker", fill = "marker",
                          palette = c("#00AFBB", "#E7B800"),repel = T,
                          xlab = "Relative abundance",ylab = "Density",
                          title = paste0(level,": ",taxa_to_test_density[[m]])
      ) +
        
        labs(fill = NULL, color = NULL)
      # Save to high-resolution PDF
      ggsave(
        filename = paste0(dens_dir,"/",level,"_",taxa_to_test_density[[m]],"_",method, ".pdf"),
        plot = x_dens,
        width = 15,
        height = 10,
        units = "in",
        dpi = 600,
        device = cairo_pdf
      )
    })
    
    ############################################################################
    
    return(physeq_list2)
  })
  
  message("Saving image")
  
  save.image(paste0(RDS_dir,"/","Comparisons_lotus_",method,".RData"))
  
}

# compare_phyloseq(dir, markers,out_comparison_dir,method="clr")
compare_phyloseq(dir, markers,out_comparison_dir,method="compositional")
# compare_phyloseq(dir, markers,out_comparison_dir,method="counts")
