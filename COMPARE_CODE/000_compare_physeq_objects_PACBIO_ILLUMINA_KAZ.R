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
source("~/SCRIPTS/000_Prepare_dirs.R")
source("~/SCRIPTS/000_wilcoxon_test_phyloseq_paired.R")
source("~/SCRIPTS/001_Accumulation_Species_threshold.R")
# source("~/SCRIPTS/002_regression_alpha_vars.R")
# source("~/SCRIPTS/003_Differentially_abundant_taxon.R")

################################################################################
set.seed(1000)
################################################################################
#Home folder
dir <- "~"
#MAKE DIR
out_comparison_dir <- paste0(dir,"/","comparisons_ITS2_PACBIO_ILLUMINA_KAZ")
if(!dir.exists(out_comparison_dir)){
  dir.create(out_comparison_dir)
}

#how to name the results
markers <- c(
  "ITS_ILLUMINA",
  "ITS_PACBIO"
)

#original_paths
paths <- c(
  "~/data/pipeline/kazakhstan-experimental/illumina/ITS2/ecm_physeq.Rdata",
  "~/data/pipeline/kazakhstan-experimental/pacbio/ITS9MUN/ITS2/ecm_physeq.Rdata"
)

compare_phyloseq <- function (dir, markers,out_comparison_dir){
  method="counts"
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
    #i <- 2
    phy_i <- 
      readRDS(paste0(dir, "/", markers[[i]], "/outcomes/rds_dir/001_filtered_PSD.RDS"))
    phy_i <- microViz::tax_fix(phy_i, unknowns = c("?"))
    phy_i@sam_data$file <- markers[[i]]
    sample_data(phy_i)$SampleID <- sample_names(phy_i)
    sample_data(phy_i)$SampleID <- sub(pattern = "_ITS",replacement = "",sample_data(phy_i)$SampleID)
    sample_data(phy_i)$SampleID <- sub(pattern = "_",replacement = "",sample_data(phy_i)$SampleID)
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
    sample_names(phy_i) <- sample_data(phy_i)$SampleID
    # 
    return(phy_i)
  })

  samp1  <- colnames(otu_table(physeq_list[[1]]));  samp2  <- colnames(otu_table(physeq_list[[2]]));
  samp2 <- sub(pattern = "_ITS",replacement = "",x = samp2)
  samp2 <- sub(pattern = "_",replacement = "",x = samp2)
  #get sample intersected
  sam_to_test <- Reduce(intersect, list(samp1,samp2))
  physeq_list <- lapply(1:length(markers), function(i) {
    x <- prune_samples(sample_names(physeq_list[[i]]) %in% sam_to_test, physeq_list[[i]])
    x <- prune_taxa(taxa_sums(x) > 0, x)
    return(x)
  })

  

  # #Aggregating rarefied object (physeq_rare)
  # physeq_rare <- lapply(1:length(markers), function(i) {
  #   phy_i <- 
  #     readRDS(paste0(dir, "/", markers[[i]], "/outcomes/rds_dir/002_rarefied_PSD_4.RDS"))
  #   phy_i <- microViz::tax_fix(phy_i, unknowns = c("?"))
  #   phy_i@sam_data$file <- markers[[i]]
  #   sample_data(phy_i)$SampleID <- sample_names(phy_i)
  #   sample_data(phy_i)$SampleID <- sub(pattern = "_ITS",replacement = "",sample_data(phy_i)$SampleID)
  #   sample_data(phy_i)$SampleID <- sub(pattern = "_",replacement = "",sample_data(phy_i)$SampleID)
  #   
  #   #get taxonomic table to be modified
  #   tax_t <- tax_table(phy_i)
  #   #get unique species
  #   UT <- unique(as.character(tax_t[,7]))
  #   #getting taxa as Russula_sp or Russula Genus and replace by "?"
  #   target_positions <- grep("(_sp$| Genus$)", UT)
  #   UT[target_positions] <- "?"
  #   UT <- UT[which(UT!="?")]
  #   tax_t[,7][!tax_t[,7] %in% UT] <- "?"
  #   tax_table(phy_i) <- tax_t
  #   phy_i <- microViz::tax_fix(phy_i,unknowns = c("?"))
  #   sample_names(phy_i) <- sample_data(phy_i)$SampleID
  #   return(phy_i)
  # })
  # 
  # samp1  <- colnames(otu_table(physeq_rare[[1]]));  samp2  <- colnames(otu_table(physeq_rare[[2]]));
  # samp2 <- sub(pattern = "_ITS",replacement = "",x = samp2)
  # samp2 <- sub(pattern = "_",replacement = "",x = samp2)
  # #get sample intersected
  # sam_to_test <- Reduce(intersect, list(samp1,samp2))
  # physeq_rare <- lapply(1:length(markers), function(i) {
  #   x <- subset_samples(physeq_rare[[i]],SampleID %in% sam_to_test)
  #   x <- prune_taxa(taxa_sums(x) > 0, x)
  #   return(x)
  # })
  
  
  #Loading original phyloseqs without rarefy
  physeq_original <- lapply(1:length(markers), function(i) {
    # i <- 1
    phy_o <- readRDS(paths[[i]])
    phy_o <- microViz::tax_fix(phy_o, unknowns = c("?"))
    sample_data(phy_o)$SampleID <- sample_names(phy_o)
    sample_data(phy_o)$SampleID <- sub(pattern = "_ITS",replacement = "",sample_data(phy_o)$SampleID)
    sample_data(phy_o)$SampleID <- sub(pattern = "_",replacement = "",sample_data(phy_o)$SampleID)
    
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
    sample_names(phy_o) <- sample_data(phy_o)$SampleID
    
    return(phy_o)
  })
  
  samp1  <- colnames(otu_table(physeq_original[[1]]));  samp2  <- colnames(otu_table(physeq_original[[2]]));
  samp2 <- sub(pattern = "_ITS",replacement = "",x = samp2)
  samp2 <- sub(pattern = "_",replacement = "",x = samp2)
  #get sample intersected
  sam_to_test <- Reduce(intersect, list(samp1,samp2))
  physeq_original <- lapply(1:length(markers), function(i) {
     x <- prune_samples(sample_names(physeq_original[[i]]) %in% sam_to_test, physeq_original[[i]])
    x <- prune_taxa(taxa_sums(x) > 0, x)
    return(x)
  })
  
  ##############################################################################
 
  ###############################################################################
  #First aggregate at species level for filtered object
  message("Aggregating phyloseq objects to species level and joining them")
  #aggregating taxonomic level using original
  phy_to_join <- lapply(1:length(physeq_list),function(j){
     # j <- 1
    phy1 <- physeq_list[[j]]
    phy1 <- microViz::tax_agg(ps = phy1,"Species")
    sample_names(phy1) <- paste0(sample_names(phy1),"_",markers[[j]])
    sample_data(phy1)$SampleID <- sample_names(phy1)
    return(phy1)
  })
################################################################################
  
  if(!file.exists(paste0(RDS_dir,"/","002_joined_PSD.RDS"))){
    phy_joined <- merge_phyloseq(phy_to_join[[1]],phy_to_join[[2]])
    samp_names <- sample_names(phy_joined)
    sample_data(phy_joined)$id <- sub(pattern = "_ITS",replacement = "",x = sample_data(phy_joined)$id )
    sample_data(phy_joined)$id <- sub(pattern = "_",replacement = "",x = sample_data(phy_joined)$id)
    un_samples <- tapply(sample_data(phy_joined)$id,sample_data(phy_joined)$id,length)
    un_samples <- un_samples[un_samples>1]
    phy_joined <- subset_samples(phy_joined, id %in% names(un_samples))
    phy_joined <- prune_taxa(taxa_sums(phy_joined) > 0, phy_joined)
    
    saveRDS(phy_joined,paste0(RDS_dir,"/","002_joined_PSD.RDS"))
    
  } else {
    phy_joined <- readRDS(paste0(RDS_dir, "/002_joined_PSD.RDS"))
    
  }

  if(!file.exists(paste0(RDS_dir,"/","002_rarefied_PSD_4joined.RDS"))){
    ps_rare <- accumulation_curve_threshold_function(phy_joined, option="SRS",threshold_label="joined",Cmin=min(sample_sums(phy_joined))) 
  } else {
    ps_rare <- readRDS(paste0(RDS_dir, "/002_rarefied_PSD_4joined.RDS"))
    
  }
  
  ##############################################################################
  ##############################################################################
  ##############################################################################
    message("Saving rarefaction curves plot (using ggrare function)")
    ggrare_plot <-
      ranacapa::ggrare(
        ps_rare,
        step = 100,
        se = TRUE,
        parallel = T,
        color = "file",
        label = "file"
      ) 
    rare_data <- ggrare_plot$data
    
    # Paso 2: Mediana y MAD por grupo
    summary_stats <- rare_data %>%
      group_by(file, Size) %>%
      summarise(
        median_richness = median(.S, na.rm = TRUE),
        mad_richness = mad(.S, na.rm = TRUE),
        .groups = "drop"
      )
    
    # Paso 3: Suavizar mediana y mad con loess por grupo
    smooth_data <- summary_stats %>%
      group_by(file) %>%
      arrange(Size) %>%
      mutate(
        smooth_median = predict(loess(median_richness ~ Size, span = 0.4)),
        smooth_mad = predict(loess(mad_richness ~ Size, span = 0.4))
      ) %>%
      ungroup()
    
    # Paso 4: Graficar
    ggrare_plot <- ggplot() +
      # Líneas individuales (una por muestra)
      geom_line(data = rare_data,
                aes(x = Size, y = .S, group = Sample, color = file),
                alpha = 0.3, size = 0.9) +
      
      # Banda: mediana ± MAD suavizada
      # geom_ribbon(data = smooth_data,
      #             aes(x = Size,
      #                 ymin = smooth_median - smooth_mad,
      #                 ymax = smooth_median + smooth_mad,
      #                 fill = file),
      #             alpha = 0.2, color = NA) +
      # 
      # # Línea de mediana suavizada
      # geom_line(data = smooth_data,
      #           aes(x = Size, y = smooth_median, color = file),
      #           size = 1.2) +
      ggsci::scale_color_jco() +
      ggsci::scale_fill_jco() +
      
      labs(x = "Number of Sequences", y = "Richness") +
      theme_minimal() +
      theme(
        legend.position = "right",
        legend.title = element_blank()
      )
    
  message("Plotting rarefaction curves")
    ggsave(
      paste0(graph_dir, "/", "001_rarefaction_curve_ggrare.pdf"),
      ggrare_plot,
      scale = 0.9,
      width = 25,
      height = 14,
      units = "in",
      dpi = 600
    )


  
  
  ##############################################################################

  
  p_r_p4_rare <-
    plot_richness(
      ps_rare,
      x = "file",
      color = "file",
      measures = c("Shannon", "Observed", "InvSimpson","Chao1"),
      title = ""#"Rarefied to twice the minimum size"
    ) +
    ylab("") +
    geom_boxplot(aes(fill = file), alpha = 0.7) +
    # stat_compare_means(
    #                    method = "wilcox.test",
    #                    size = 8,
    #                    aes(label = paste0("p = ", after_stat(p.format)))) +
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
  
  
    ggsave(
      paste0(graph_dir, "/", "002_Alpha_boxplots.pdf"),
      p_r_p4_rare,
      width = 20,
      height = 25,
      units = "in",
      dpi = 600
    )
  ################################################################################
  #Joining objects in one
  physeq_clr2 <- microbiome::transform(phy_joined, "compositional")
  
    ################################################################################
  # Perform PCoA ordination using Bray-Curtis distance
  psd5.mds.euc2 <-
    # ordinate(physeq_clr2, method = "PCoA", distance = "euclidean")
      ordinate(physeq_clr2, method = "PCoA", distance = "bray")
  evals <- psd5.mds.euc2$values$Eigenvalues
  
  
  # Extraer los datos de la gráfica
  ord_df <- as.data.frame(plot_ordination(physeq_clr2, psd5.mds.euc2, justDF = TRUE))
  ord_df$file <- sample_data(physeq_clr2)$file
  
  # Obtener los polígonos convexos
  hulls <- ord_df %>%
    group_by(file) %>%
    slice(chull(Axis.1, Axis.2))  # Convex hull para cada grupo
  # Create the ordination plot and apply the 'jco' color palette
  pord2 <-
    plot_ordination(physeq_clr2, psd5.mds.euc2, color = "file", shape = "file") +
    geom_polygon(data = hulls, aes(x = Axis.1, y = Axis.2, fill = file), alpha = 0.2, color = NA) +
    
    # labs(col = "Fungi System") +
    geom_point(size = 8) +  # Increase point size
    
    ggsci::scale_color_jco() +         # Color de puntos
    ggsci::scale_fill_jco() +          # Color de polígonos (igual que los puntos)
    coord_fixed(sqrt(evals[2] / evals[1])) +
    #theme_pubr(base_size = 20) +  # Increase all base text size
    theme(
      axis.title = element_text(size = 22),
      axis.text = element_text(size = 20),
      strip.text = element_text(size = 20),
      legend.title = element_text(size = 0),
      legend.text = element_text(size = 18)
    )
  # Show the plot
  message("Beta diversity...Saved in graphics subfolder")
  ggsave(
    paste0(graph_dir, "/", "004_Beta_PCoA_EUC.pdf"),
    pord2,
    width = 20,
    height = 25,
    units = "in",
    dpi = 600
  )
  ##############################################################################
  ##############################################################################
  Alpha_N <- phyloseq::estimate_richness(ps_rare)
  Alpha_N$pielou <- evenness(ps_rare, index = "pielou")[, 1]
  calc <- ChaoRichness(x = otu_table(ps_rare), datatype = "abundance", conf = 0.95)
  Alpha_N$Chao1 <- calc$Estimator
  Alpha_N$se.chao1 <- calc$Est_s.e.
  
  #Adding sample data
  Alpha_N <- cbind(sample_data(ps_rare), Alpha_N)
  
  message("Saving Alpha diversity indexes CSV files")
    write.csv(
    Alpha_N,
    paste0(csv_dir, "/", "002_Alpha.csv"),
    na = "",
    row.names = F,
    quote = T
  )


  ##############################################################################
  ##############################################################################
    #SUMMARY TABLE
    message("creating summary, be patient")
    summary <- list() 
      
#      lapply(1:length(markers), function(i) {
    for(i in 1:length(markers)){
      # phy_o <- physeq_original[[i]]#readRDS(paths[[i]])
      #rarifed
      # phy_i <-   physeq_rare[[i]] #readRDS(paste0(dir, "/", markers[[i]], "/outcomes/rds_dir/002_rarefied_PSD_4.RDS"))       
      counts <- sample_sums(physeq_list[[i]])
      ps_rare_joined <- subset_samples(ps_rare, file==markers[[i]])
      ps_rare_joined2 <- prune_taxa(taxa_sums(ps_rare_joined) > 0, ps_rare_joined)

      x <- data.frame(
        marker = markers[[i]],
        orig_nsamples = length(sample_sums(physeq_original[[i]])),
        orig_read_counts = sum(sample_sums(physeq_original[[i]])),
        orig_min = min(sample_sums(physeq_original[[i]]),na.rm=T),
        orig_max =  max(sample_sums(physeq_original[[i]]),na.rm=T),
        orig_min_sample = names(sample_sums(physeq_original[[i]])[which.min(sample_sums(physeq_original[[i]]))]),
        orig_max_sample = names(sample_sums(physeq_original[[i]])[which.max(sample_sums(physeq_original[[i]]))]),

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
        orig_sp_no_na = 
          length(unique(na.omit(unlist(lapply(1:nrow(tax_table(physeq_original[[i]])),function(j){
            x <- paste(tax_table(physeq_original[[i]])[j,7],collapse = "-")
            x2 <- strsplit(x," ")
            if(length(x2[[1]])>1){
              x <- NA    
            } else {
              x <- x
            }
            return(x)
          }))))),

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
        filtered_sp_no_na = 
          length(unique(na.omit(unlist(lapply(1:nrow(tax_table(physeq_list[[i]])),function(j){
            x <- paste(tax_table(physeq_list[[i]])[j,7],collapse = "-")
            x2 <- strsplit(x," ")
            if(length(x2[[1]])>1){
              x <- NA    
            } else {
              x <- x
            }
            return(x)
          }))))),
        rare_threshold=min(sample_sums(ps_rare_joined2)),
        rare_nsamples = length(Alpha_N$file[which(Alpha_N$file==markers[[i]])]),

        rare_families =
          length(unique(unlist(lapply(1:nrow(tax_table(ps_rare_joined2)),function(j){
            x <- paste(tax_table(ps_rare_joined2)[j,1:5],collapse = "-")
            return(x)
          })))),
        rare_genus =
          length(unique(unlist(lapply(1:nrow(tax_table(ps_rare_joined2)),function(j){
            x <- paste(tax_table(ps_rare_joined2)[j,1:6],collapse = "-")
            return(x)
          })))),
        rare_sp =
          length(unique(unlist(lapply(1:nrow(tax_table(ps_rare_joined2)),function(j){
            x <- paste(tax_table(ps_rare_joined2)[j,1:7],collapse = "-")
            return(x)
          })))),
        rare_sp_no_na = 
          length(unique(na.omit(unlist(lapply(1:nrow(tax_table(ps_rare_joined2)),function(j){
            x <- paste(tax_table(ps_rare_joined2)[j,7],collapse = "-")
            x2 <- strsplit(x," ")
            if(length(x2[[1]])>1){
              x <- NA    
            } else {
              x <- x
            }
            return(x)
          }))))),
        rare_obs_sp_av = mean(Alpha_N$Observed[Alpha_N$file==markers[[i]]],na.rm=T),
        rare_obs_sp_sd = sd(Alpha_N$Observed[Alpha_N$file==markers[[i]]],na.rm=T),
        rare_Shannon_av = mean(Alpha_N$Shannon[Alpha_N$file==markers[[i]]],na.rm=T),
        rare_Shannon_sd = sd(Alpha_N$Shannon[Alpha_N$file==markers[[i]]],na.rm=T),
        rare_InvSimpson_av = mean(Alpha_N$InvSimpson[Alpha_N$file==markers[[i]]],na.rm=T),
        rare_InvSimpson_sd = sd(Alpha_N$InvSimpson[Alpha_N$file==markers[[i]]],na.rm=T),
        rare_Pielou_av = mean(Alpha_N$pielou[Alpha_N$file==markers[[i]]],na.rm=T),
        rare_Pielou_sd = sd(Alpha_N$pielou[Alpha_N$file==markers[[i]]],na.rm=T),
        rare_Chao1_av = mean(Alpha_N$Chao1[Alpha_N$file==markers[[i]]],na.rm = T),
        rare_Chao1_sd = sd(Alpha_N$Chao1[Alpha_N$file==markers[[i]]],na.rm = T)

      )
      summary[[i]] <- x
    }

    summary <- do.call(rbind,summary)
    write.csv(summary,paste0(out_comparison_dir, "/comparison_summary.csv"),row.names = F)



  ##############################################################################
  ##############################################################################
     # library(MicrobiotaProcess)
     # mpse4 <- ps_rare %>% as.MPSE()
     # detach("package:MicrobiotaProcess", unload = TRUE)
    theme_custom <- theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(size = 14),
      strip.text = element_text(size = 18, face = "bold"),
      # Facet labels (AM, ECM)
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14),
      plot.title = element_text(size = 16, face = "bold")
    )
    
    p1 <-
      ps_rare  #phyloseq::merge_samples(ps_rare, group = "Fungi_system")
    p1 <- tax_fix(p1, unknowns = c("?"))
    
    #Order barplot
    bp1 <- microViz::comp_barplot(
      tax_fix(p1),
      tax_level = "Order",
      bar_outline_colour = NA,
      sample_order = "bray",
      facet_by = "file",
      bar_width = 0.7,
      taxon_renamer = toupper
    ) +
      coord_flip() +
      # facet_wrap("Fungi_system", nrow = 1, scales = "free") +
      labs(x = NULL, y = NULL) +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
      theme_custom
    
    #family barplot
    
    bp2 <- microViz::comp_barplot(
      tax_fix(p1),
      tax_level = "Family",
      n_taxa = 10,
      bar_outline_colour = NA,
      sample_order = "bray",
      facet_by = "file",
      bar_width = 0.7,
      taxon_renamer = toupper
    ) +
      coord_flip() +
      # facet_wrap("Fungi_system", nrow = 1, scales = "free") +
      labs(x = NULL, y = NULL) +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
      theme_custom
    
    
    #genus barplot
    
    bp3 <- microViz::comp_barplot(
      tax_fix(p1),
      tax_level = "Genus",
      n_taxa = 10,
      bar_outline_colour = NA,
      facet_by = "file",
      sample_order = "bray",
      bar_width = 0.7,
      taxon_renamer = toupper
    ) +
      coord_flip() +
      # facet_wrap("Fungi_system", nrow = 1, scales = "free") +
      labs(x = NULL, y = NULL) +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
      theme_custom
    
    #species barplot
    
    bp4 <- microViz::comp_barplot(
      tax_fix(p1),
      tax_level = "Species",
      n_taxa = 10,
      bar_outline_colour = NA,
      facet_by = "file",
      sample_order = "bray",
      bar_width = 0.7,
      taxon_renamer = toupper
    ) +
      coord_flip() +
      # facet_wrap("Fungi_system", nrow = 1, scales = "free") +
      labs(x = NULL, y = NULL) +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
      theme_custom
    
    
    bp_total <- gridExtra::grid.arrange(bp1,
                                        bp2, bp3, bp4)
    
    
    ggsave(
      paste0(graph_dir, "/", "005_BARPLOT_top10_Rare.pdf"),
      bp_total,
      width = 20,
      height = 25,
      units = "in",
      dpi = 600
    )
    
    ggsave(
      paste0(graph_dir, "/", "005_BARPLOT_top10_Genus.pdf"),
      bp3,
      width = 20,
      height = 25,
      units = "in",
      dpi = 600
    )
    
    ggsave(
      paste0(graph_dir, "/", "005_BARPLOT_top10_Species.pdf"),
      bp4,
      width = 20,
      height = 25,
      units = "in",
      dpi = 600
    )
    
    
  ##############################################################################
  ##############################################################################
    
  levels <-
    c("Family", "Genus", "Species")
  
  x_jacc <- lapply(1:length(levels), function(i) {
    # 
     # i <- 1
    level <- levels[[i]]
    message(paste0("Level: ",level))
    ##############################################################################
    ##############################################################################
    ##############################################################################
    #aggregating taxonomic level using original
    phy_original <- lapply(1:length(physeq_original),function(j){
      phy1 <- physeq_original[[j]]

        phy1 <- microViz::tax_agg(ps = phy1, level)
        if(method!="counts"){
          phy1 <- microbiome::transform(phy1, method)          
        } else {
          phy1 <- phy1
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
      
      sample_list[[2]] <- sub(pattern = "_ITS",replacement = "",x = sample_list[[2]])
      sample_list[[2]] <- sub(pattern = "_",replacement = "",x = sample_list[[2]])
      
      colnames(x_list[[2]]) <- sub(pattern = "_ITS",replacement = "",x = colnames(x_list[[2]]))
      colnames(x_list[[2]]) <- sub(pattern = "_",replacement = "",x = colnames(x_list[[2]]))
      
      #get sample intersected
      sam_to_test <- Reduce(intersect, sample_list)    
      sam_to_test <- sam_to_test[sam_to_test!="dataset"]
      # sam_to_test <- sam_to_test[!sam_to_test %in% c("NC_1","NC_2")]
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
  if(nrow(wilc_test)>0){
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
  }
    ##############################################################################
    #Obtaining data from two datasets to create density plots
    if(nrow(wilc_test)>0){
      plot_density_function(dens_dir,wilc_test,markers,x_list,sam_to_test,method,level,dataset="original")      
    } else {
      message("no Diff taxa for original Phyloseq object")
    }
    ##############################################################################
    #Aggregating taxa using rarefied objects
    message("Aggregating rarefied data to selected level")
    physeq_list2 <- lapply(1:length(physeq_list), function(j) {
      # print(j)
      phy1 <- physeq_list[[j]]
        phy1 <- microViz::tax_agg(ps = phy1, level)
      
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
      # j <- 2
      p1 <-physeq_list2[[j]]
      # p1 <- subset_samples(p1, id %in% names(un_samples))
      p1 <- prune_taxa(taxa_sums(p1) > 0, p1)
      tax1 <- tax_table(p1)
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

        phy1 <- microViz::tax_agg(ps = phy1, level)
        if(method!="counts"){
          phy1 <- microbiome::transform(phy1, method)          
        } else {
          phy1 <- phy1
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
      
      sample_list[[2]] <- sub(pattern = "_ITS",replacement = "",x = sample_list[[2]])
      sample_list[[2]] <- sub(pattern = "_",replacement = "",x = sample_list[[2]])
      
      colnames(x_list[[2]]) <- sub(pattern = "_ITS",replacement = "",x = colnames(x_list[[2]]))
      colnames(x_list[[2]]) <- sub(pattern = "_",replacement = "",x = colnames(x_list[[2]]))
      
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
  
    if(nrow(wilc_test)>0){
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
  }  
    #Obtaining data from two datasets to create density plots
    if(nrow(wilc_test)>0){
      plot_density_function(dens_dir,wilc_test,markers,x_list,sam_to_test,method,level,dataset="filtered")      
    } else {
      message("no Diff taxa for filtered Phyloseq object")
    }

    ##############################################################################
    ##############################################################################
    ##############################################################################
    #WILCOXON RAREFIED!
    message("Wilcoxon test for rarefied phyloseqs")
    
    
    physeq_rare_2 <- lapply(1:length(markers),function(j){
      # j <- 1
      phy1 <- subset_samples(ps_rare, file==markers[[i]])
      phy1 <- prune_taxa(taxa_sums(phy1) > 0, phy1)
      phy1 <- microViz::tax_agg(ps = phy1, level)
      sample_names(phy1) <- sample_data(phy1)$id
        if(method!="counts"){
          phy1 <- microbiome::transform(phy1, method)          
        } else {
          phy1 <- phy1
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
      
      sample_list[[2]] <- sub(pattern = "_ITS",replacement = "",x = sample_list[[2]])
      sample_list[[2]] <- sub(pattern = "_",replacement = "",x = sample_list[[2]])
      
      colnames(x_list[[2]]) <- sub(pattern = "_ITS",replacement = "",x = colnames(x_list[[2]]))
      colnames(x_list[[2]]) <- sub(pattern = "_",replacement = "",x = colnames(x_list[[2]]))
      
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
    
  if(nrow(wilc_test)>0){
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
  }  
    ############################################################################
    message("detecting taxa to show density differences")
    #taxa to use for density plots
    taxa_to_test_density <- wilc_test$taxon
    
    #Obtaining data from two datasets to create density plots
    #Obtaining data from two datasets to create density plots
    if(nrow(wilc_test)>0){
      plot_density_function(dens_dir,wilc_test,markers,x_list,sam_to_test,method,level,dataset="rarified")      
    } else {
      message("no Diff taxa for rarified Phyloseq object")
    }
    
    ############################################################################
    
    
    #obtain unique taxa per dataset and object and saving a CSV file
    if(level %in% levels[2:3]){
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
      
    # library(MicrobiotaProcess)
    # mpse4 <- esophagus %>% as.MPSE() 
    # detach("package:MicrobiotaProcess", unload = TRUE)
    }
    return(physeq_list2)
  })
  
  message("Saving image")
  
  save.image(paste0(RDS_dir,"/","Comparisons_lotus_",method,".RData"))
  
}

# compare_phyloseq(dir, markers,out_comparison_dir,method="clr")
compare_phyloseq(dir, markers,out_comparison_dir)
# compare_phyloseq(dir, markers,out_comparison_dir,method="counts")
