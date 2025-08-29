#' Generate and Save Rarefaction Curves with Optional Rarefaction
#'
#' This function generates rarefaction curves from a `phyloseq` object using either
#' the `iNEXT` or `ggrare` methods. It also performs sample rarefaction based on
#' estimated or user-defined thresholds and saves the output as an RDS object.
#'
#' @param physeq4 A `phyloseq` object with OTU/ASV abundance data.
#' @param option Character. Method to use for rarefaction:
#' `"iNEXT"` for interpolation/extrapolation curves,
#' or `"ggrare"` (from `ranacapa`) for a traditional approach.
#'
#' @return A rarefied `phyloseq` object.
#' @export
#'
#' @details
#' - If `option = "iNEXT"`, the function estimates diversity and thresholds using
#'   interpolation/extrapolation methods and automatically rarefies the samples.
#' - If `option = "ggrare"`, the user must input a rarefaction threshold interactively.
#'
#' Plots are saved as high-resolution PDFs under the folder defined by `graph_dir`, and
#' the rarefied `phyloseq` object is saved as `"002_rarefied_PSD.RDS"` under `RDS_dir`.
#'
#' @importFrom iNEXT iNEXT estimateD
#' @importFrom ggpubr ggsave
#' @importFrom ggplot2 geom_line geom_vline xlim theme
#' @importFrom dplyr filter group_by slice_max ungroup
#' @importFrom ggrepel geom_text_repel
#' @importFrom phyloseq sample_sums rarefy_even_depth otu_table
#' @importFrom ranacapa ggrare
#' @importFrom scales label_number
#'
#' @examples
#' \dontrun{
#' ps_rare <- accumulation_curve_function(physeq4, option = "iNEXT")
#' }
#'


# ################################################################################
# #Load pacman to load complementary libraries
# library(pacman);
# #Load libraries or install if it is the case from CRAN and github respetively
pacman::p_load(phyloseq, ggpubr, iNEXT,SRS)
################################################################################
################################################################################
accumulation_curve_function <- function(physeq4, option) {
  if (is.null(option) | !option %in% c("iNEXT", "ggrare","plateau","SRS")) {
    stop("Please use a valid option")
  }
  ##############################################################################
  #Adding a set seed to get reproducible results
  set.seed(1000)
  ##############################################################################
  if (option == "iNEXT") {
    ############################################################################
    ############################################################################
    ############################################################################
    ############################################################################
    message("Obtaining OTU table from the Phyloseq Object")
    #Obtaining OTU table to try iNEXT
    otu_tab <- as.data.frame(otu_table(physeq4))
    for (i in 1:ncol(otu_tab)) {
      otu_tab[, i] <- as.numeric(otu_tab[, i])
    }
    rm(i)
    
    #https://github.com/adlape95/ME/blob/main/Stats-and-Figures.Rmd
    # Only running q=0
    
    message("Using iNEXT to get Interpolation and extrapolation approach, be patient!")
    
    out <- iNEXT::iNEXT(otu_tab, q = 0,
                        #size=seq(1, 10000, by=500),
                        datatype = "abundance")
    message("Using estimateD with abundance data \n
          and size option to obtain threshold to rarefy")
    message(
      "This function computes the diversity estimates\n
          for the minimum among all doubled reference sample sizes"
    )
    out1 <-
      estimateD(otu_tab,
                q = c(0),
                datatype = "abundance",
                base = "size")
    
    ##############################################################################
    #Getting labels to add into the iNEXT plot
    message("Getting labels to add into the iNEXT plot")
    
    label_data <- out$iNextEst$size_based %>%
      filter(Order.q == 0) %>%
      group_by(Assemblage) %>%
      slice_max(qD, n = 1, with_ties = FALSE) %>%
      ungroup()
    
    message("Saving rarefaction curves plot (dashed lines are extrapolated values)")
    #Plotting
    Rareplot1 <- ggiNEXT(x = out, type = 1) +
      #theme_bw() +
      geom_line(size = 0.1) +
      #geom_point(size=0, na.rm = TRUE) +
      xlim(c(0, max(sample_sums(physeq4)))) +
      xlab("Number of sequences") + ylab("Richness (Observed OTUs)")  +
      geom_text_repel(
        data = label_data,
        aes(
          x = m,
          y = qD,
          label = Assemblage,
          color = Assemblage
        ),
        #hjust = -0.1,
        size = 3,
        fontface = "italic",
        show.legend = FALSE,
        max.overlaps = Inf,
        # Muestra todas las etiquetas
        direction = "y",
        # Solo se repelen en el eje Y
        nudge_x = 1000,
        # Corre las etiquetas un poco a la derecha
        min.segment.length = 0     # Siempre muestra líneas si hay nudges
      ) +
      
      # # Cuadrícula más densa
      # scale_x_continuous(
      #   breaks = seq(0, max(sample_sums(physeq4)), by = 20000),
      #   minor_breaks = seq(0, max(sample_sums(physeq4)), by = 5000),
      #   labels = scales::label_number()
      # ) +
      # scale_y_continuous(
      #   breaks = seq(0, max(out$iNextEst$size_based$qD), by = 50),
      #   minor_breaks = seq(0, max(out$iNextEst$size_based$qD), by = 10),
      #   labels = scales::label_number()
      # ) +
      geom_vline(
        xintercept = unique(out1$m),
        linetype = "dotted",
        color = "black",
        size = 0.5
      ) +
      
      
      theme(legend.position = "None")
    
    
    message(paste0("Plot saved as: 001_rarefaction_curve_iNEXT.pdf"))
    ggsave(
      paste0(graph_dir, "/", "001_rarefaction_curve_iNEXT.pdf"),
      Rareplot1,
      scale = 0.9,
      width = 25,
      height = 14,
      units = "in",
      dpi = 600
    )
    ##############################################################################
    #rarefying to selected value or suggested value
    message(paste("Rarefaction obtained using ", unique(out1$m), "value"))
    ps_rare <- phyloseq::rarefy_even_depth(
      physeq4,
      #sample.size = seqs_to[[i]],  # Or choose a fixed number
      sample.size = unique(out1$m),
      # Or choose a fixed number
      rngseed = 123,
      # For reproducibility
      replace = FALSE,
      # No replacement
      trimOTUs = TRUE,
      # Remove unobserved taxa
      verbose = FALSE
    )
    
    
    
    #Saving Phyloseq object
    message("Phyloseq object rarefied obtained...Saving to RDS subfolder")
    saveRDS(ps_rare, paste0(RDS_dir, "/002_rarefied_PSD_1.RDS"))
    message(paste0("Sample with: minimum read counts: ",min(sample_sums(physeq4)),
                   " /Threshold used: ",unique(sample_sums(ps_rare))))
    write.csv(data.frame(method="iNEXT",threshold=unique(sample_sums(ps_rare))),
              paste0(csv_dir,"/001_threshold_1.csv"),row.names = F)
  } else if (option == "ggrare") {
    ############################################################################
    ############################################################################
    ############################################################################
    ############################################################################
    message("Saving rarefaction curves plot (using ggrare function)")
    ggrare_plot <-
      ranacapa::ggrare(
        physeq4,
        step = 500,
        se = TRUE,
        parallel = T,
        color = "Sample",
        label = "Sample"
      )
    ggrare_plot  <- ggrare_plot + theme(legend.position = "none")
    ggrare_plot <- ggrare_plot +
      labs(title = "",
           x = "Number of Sequences",
           y = "Richness")
    ggsave(
      paste0(graph_dir, "/", "001_rarefaction_curve_ggrare.pdf"),
      ggrare_plot,
      scale = 0.9,
      width = 25,
      height = 14,
      units = "in",
      dpi = 600
    )
    
    message("Plotting rarefaction curves")
    ggrare_plot
    
    val_rare <-
      as.numeric(readline(prompt = "Choose a value to rarify: "))
    
    #Rarefy using the second mininum number of samples
    ps_rare <- phyloseq::rarefy_even_depth(
      physeq4,
      #sample.size = seqs_to[[i]],  # Or choose a fixed number
      sample.size = val_rare,
      # Or choose a fixed number
      rngseed = 123,
      # For reproducibility
      replace = FALSE,
      # No replacement
      trimOTUs = TRUE,
      # Remove unobserved taxa
      verbose = FALSE
    )
    #Saving Phyloseq object
    message("Phyloseq object rarefied obtained...Saving to RDS subfolder")
    saveRDS(ps_rare, paste0(RDS_dir, "/002_rarefied_PSD_2.RDS"))
    write.csv(data.frame(method="ggrare",threshold=val_rare),
              paste0(csv_dir,"/001_threshold_2.csv"),row.names = F)
    
    
  } else if(option == "plateau") {
    ############################################################################
    ############################################################################
    ############################################################################
    ############################################################################
    
    #Threshold
    threshold <- 2 # It means difference of two OTUs
    #Consecutive number of differences to test
    n_consec <- 3
    
    #Filtering to use only rarified data from iNEXT
    message("Using Plateau based approach")
    
    message("Obtaining OTU table from the Phyloseq Object")
    #Obtaining OTU table to try iNEXT
    otu_tab <- as.data.frame(otu_table(physeq4))
    for (i in 1:ncol(otu_tab)) {
      otu_tab[, i] <- as.numeric(otu_tab[, i])
    }
    rm(i)
    
    message("Using iNEXT to get Interpolation and extrapolation approach, be patient!")
    
    out <- iNEXT::iNEXT(otu_tab, q = 0,
                        #size=seq(1, 10000, by=500),
                        datatype = "abundance")
    message("Getting size based data frame and ommiting Extrapolation")
    df_raref <- out$iNextEst$size_based %>%
      filter(Method != "Extrapolation")  # focus on observed range
    df_raref$diff_qD <- NA
    df_raref$threshold <- NA
    df_raref$idx <- NA
    
    #Getting unique Assemblages
    assemblage <- unique(df_raref$Assemblage)
    #Obtaining differences < 2 OTUS as input to calculate threshold to samples
    
    
    message(paste0("Testing with: ", n_consec," differences"))
    
    ##################################
    #Obtaining the consecutive number
    d1_thr <- list()
    for(i in 1:length(assemblage)){
      
      #Filtering by assemblage
      assem <- assemblage[i]
      
      #Filtering data using the assemblage
      df_sub <- df_raref %>% filter(Assemblage == assem)
      
      #Calculate differences of qD (Observed species due to q=0)
      d_qD <- diff(df_sub$qD)
      #Adding differences to subset data frame
      df_sub$diff_qD <- c(NA,d_qD)
      df_sub$threshold <- df_sub$diff_qD  < threshold
      df_sub$threshold[which(is.na(df_sub$diff_qD))] <- FALSE
      
      #Obtaining the index to get the provided threshold
      df_sub$idx <- (df_sub$threshold*1) * 1:nrow(df_sub)
      #Getting optimal read counts from consecutive differences
      #If there are no n_consec then ommiting that sample for calculations
      if(length(df_sub$idx>0) >=n_consec){
        x <- df_sub[df_sub$idx>0,]
        d1_thr[[i]] <-x[3,"m"]
      } else {
        d1_thr[[i]] <- NA
      }
    }  
    #Joining in one data.frame
    d1_thr <- do.call(rbind,d1_thr)
    
    message(paste0("Calculating threshold as the Median of samples (Read counts where ",n_consec, "differences)"))
    threshold_final <- median(d1_thr,na.rm = T)
    
    
    
    ##############################################################################
    #Getting labels to add into the iNEXT plot
    message("Getting labels to add into the iNEXT plot")
    
    label_data <- out$iNextEst$size_based %>%
      filter(Order.q == 0) %>%
      group_by(Assemblage) %>%
      slice_max(qD, n = 1, with_ties = FALSE) %>%
      ungroup()
    
    message("Saving rarefaction curves plot (dashed lines are extrapolated values)")
    #Plotting
    Rareplot1 <- ggiNEXT(x = out, type = 1) +
      #theme_bw() +
      geom_line(size = 0.1) +
      #geom_point(size=0, na.rm = TRUE) +
      xlim(c(0, max(sample_sums(physeq4)))) +
      xlab("Number of sequences") + ylab("Richness (Observed OTUs)")  +
      geom_text_repel(
        data = label_data,
        aes(
          x = m,
          y = qD,
          label = Assemblage,
          color = Assemblage
        ),
        #hjust = -0.1,
        size = 3,
        fontface = "italic",
        show.legend = FALSE,
        max.overlaps = Inf,
        # Muestra todas las etiquetas
        direction = "y",
        # Solo se repelen en el eje Y
        nudge_x = 1000,
        # Corre las etiquetas un poco a la derecha
        min.segment.length = 0     # Siempre muestra líneas si hay nudges
      ) +
      
      # # Cuadrícula más densa
      # scale_x_continuous(
      #   breaks = seq(0, max(sample_sums(physeq4)), by = 20000),
      #   minor_breaks = seq(0, max(sample_sums(physeq4)), by = 5000),
      #   labels = scales::label_number()
      # ) +
      # scale_y_continuous(
      #   breaks = seq(0, max(out$iNextEst$size_based$qD), by = 50),
      #   minor_breaks = seq(0, max(out$iNextEst$size_based$qD), by = 10),
      #   labels = scales::label_number()
      # ) +
      geom_vline(
        xintercept = threshold_final,
        linetype = "dotted",
        color = "black",
        size = 0.5
      ) +
      
      theme(legend.position = "None")
    
    message(paste0("Plot saved as: 001_rarefaction_curve_iNEXT_Plateau.pdf"))
    ggsave(
      paste0(graph_dir, "/", "001_rarefaction_curve_iNEXT_Plateau.pdf"),
      Rareplot1,
      scale = 0.9,
      width = 25,
      height = 14,
      units = "in",
      dpi = 600
    )
    ##############################################################################
    #rarefying to selected value or suggested value
    message(paste("Rarefaction obtained using ", threshold_final, "value"))
    ps_rare <- phyloseq::rarefy_even_depth(
      physeq4,
      #sample.size = seqs_to[[i]],  # Or choose a fixed number
      sample.size =threshold_final,
      # Or choose a fixed number
      rngseed = 123,
      # For reproducibility
      replace = FALSE,
      # No replacement
      trimOTUs = TRUE,
      # Remove unobserved taxa
      verbose = FALSE
    )
    
    
    #Saving Phyloseq object
    message("Phyloseq object rarefied obtained...Saving to RDS subfolder")
    saveRDS(ps_rare, paste0(RDS_dir, "/002_rarefied_PSD_3.RDS"))
    message("returning rarefied file to calculate Alpha diversity")
    message(paste0("Sample with: minimum read counts: ",min(sample_sums(physeq4)),
                   " /Threshold used: ",unique(sample_sums(ps_rare))))
    write.csv(data.frame(method="plateau",threshold=unique(sample_sums(ps_rare))),
              paste0(csv_dir,"/001_threshold_3.csv"),row.names = F)
    
  } else if(option=="SRS"){
    ############################################################################
    ############################################################################
    ############################################################################
    ############################################################################
    #Filtering to use only rarified data from iNEXT
    message("Using SRS based approach")
    
    message("Obtaining OTU table from the Phyloseq Object")
    otu_tab <- as.data.frame(otu_table(physeq4))
    for (i in 1:ncol(otu_tab)) {
      otu_tab[, i] <- as.numeric(otu_tab[, i])
    }
    rm(i)
    #Obtaining otu mininum sample to minimum size sample
    Cmin <- min(colSums(otu_tab))*2
    message("Calculating SRS and assigning to a rarefied Phyloseq object")
    #Calculating SRS and assigning to a rarefied Phyloseq object
    SRS_output <- SRS(otu_tab, Cmin,seed = 1000)
    row.names(SRS_output) <- row.names(otu_table(physeq4))
    SRS_output <- as.data.frame(SRS_output)
    #Assigning rarefied object
    phy5 <- physeq4
    otu_table(phy5) <- otu_table(SRS_output, taxa_are_rows = T)
    
    ggrare_plot <-
      ranacapa::ggrare(
        physeq4,
        step = 500,
        se = TRUE,
        parallel = T,
        color = "Sample",
        label = "Sample"
      )
    ggrare_plot  <- ggrare_plot + theme(legend.position = "none")
    ggrare_plot <- ggrare_plot +
      labs(title = "",
           x = "Number of Sequences",
           y = "Richness")
    ggsave(
      paste0(graph_dir, "/", "001_rarefaction_curve_ggrare.pdf"),
      ggrare_plot,
      scale = 0.9,
      width = 25,
      height = 14,
      units = "in",
      dpi = 600
    )
    
    message("Plotting rarefaction curves")
    # ggrare_plot
    
    ##############################################################################
    
    SRSCurve_rare <- 
      SRScurve(otu_tab, metric = "richness", step = 50, sample = 200, max.sample.size = max(colSums(otu_tab)), 
             rarefy.comparison = TRUE, rarefy.repeats = 50,
             rarefy.comparison.legend = TRUE, ylab = "richness",
             col = c(rep(c("#000000", "#E69F00", "#56B4E9"),2)),
             lty = c(1,2))
    
    
    pdf(paste0(graph_dir, "/", "001_rarefaction_curve_SRS.pdf"))
    replayPlot(SRSCurve_rare)
    dev.off()

    #Saving Phyloseq object
    message("Phyloseq object rarefied obtained...Saving to RDS subfolder")
    ps_rare <- phy5
    ps_rare <- phyloseq::prune_taxa(taxa_sums(ps_rare) > 0, ps_rare)
    saveRDS(ps_rare, paste0(RDS_dir, "/002_rarefied_PSD_4.RDS"))
    message("returning rarefied file to calculate Alpha diversity")
    message(paste0("Sample with: minimum read counts: ",min(sample_sums(phy5)),
                   " /Threshold used: ",unique(sample_sums(phy5))))
    write.csv(data.frame(method="SRS",threshold=unique(sample_sums(ps_rare))),
              paste0(csv_dir,"/001_threshold_4.csv"),row.names = F)
    
    
  }
  

  
  ##############################################################################
  
  return(ps_rare)
}



  
  
  
#   
#   ps_rare <- phyloseq::rarefy_even_depth(physeq4,
#                                          #sample.size = seqs_to[[i]],  # Or choose a fixed number
#                                          sample.size = threshold_final,  # Or choose a fixed number
#                                          rngseed = 123,                        # For reproducibility
#                                          replace = FALSE,                      # No replacement
#                                          trimOTUs = TRUE,                      # Remove unobserved taxa
#                                          verbose = FALSE)
#   
#   
#   ##############################################################################
#   #Getting labels to add into the iNEXT plot
#   message("Getting labels to add into the iNEXT plot")
#   
#   label_data <- out$iNextEst$size_based %>%
#     filter(Order.q == 0) %>%
#     group_by(Assemblage) %>%
#     slice_max(qD, n = 1, with_ties = FALSE) %>%
#     ungroup()
#   
#   
#   
#   
#   #Saving Phyloseq object
#   message("Phyloseq object rarefied obtained...Saving to RDS subfolder")
#   saveRDS(ps_rare, paste0(RDS_dir, "/002_rarefied_PSD.RDS"))
#   message(paste0("Sample with: minimum read counts: ",min(sample_sums(physeq4)),
#                  " /Threshold used: ",unique(sample_sums(ps_rare))))
#   
#   
#   
#   ggrare_plot <-
#     ranacapa::ggrare(
#       ps_rare,
#       step = 500,
#       se = TRUE,
#       parallel = T,
#       color = "Sample",
#       label = "Sample"
#     )
#   ggrare_plot  <- ggrare_plot + theme(legend.position = "none")
#   ggrare_plot <- ggrare_plot +
#     labs(title = "",
#          x = "Number of Sequences",
#          y = "Richness")
#   ggsave(
#     paste0(graph_dir, "/", "001_rarefaction_curve.pdf"),
#     ggrare_plot,
#     scale = 0.9,
#     width = 25,
#     height = 14,
#     units = "in",
#     dpi = 600
#   )
#   
#   
#    
#   
# #   
# #   
# #   #Calculating Differences
# #   d1 <- diff(df_raref$qD[df_raref$Assemblage==assemblage[[i]]])
# #   #Obtaining input data to get m and qD for approach
# #   thr <- data.frame(Assemblage = assemblage[[i]] ,
# #                     seqs = as.numeric(sample_sums(physeq4)[assemblage[[i]]]),
# #                     m=df_raref$m[df_raref$Assemblage==assemblage[[i]]][-1],
# #                     qD = df_raref$qD[df_raref$Assemblage==assemblage[[i]]][-1],
# #                     diff=d1)
# #   x <- thr[which(thr$diff< threshold),]
# #   plot(sqrt(x$m),sqrt(x$qD))
# # # 
# # #   d1_thr[[i]] <- thr
# # #   
# #   qD_i= thr$qD[which.min(thr$m[thr$diff < threshold])]
# #   # m_i = min(thr$m[thr$diff < threshold],na.rm = T)
# #   
# #   if(length(qD_i)==0){
# #     d1_thr[[i]] <- data.frame(Assemblage = assemblage[[i]],
# #                               seqs = as.numeric(sample_sums(physeq4)[assemblage[[i]]]),
# #                               qD=NA,
# #                               m=NA)
# #   } else {
# #     d1_thr[[i]] <- data.frame(Assemblage = assemblage[[i]],
# #                               seqs = as.numeric(sample_sums(physeq4)[assemblage[[i]]]),
# #                               qD=qD_i,
# #                               m=m_i)
# #   }
# #   
# #   
# }
# 
# 
# # p1 <- ggpubr::ggline(out$iNextEst$size_based, x = "m", y = "qD", group  = "Assemblage",color = "Assemblage",point.size = 0.2,legend=NULL) + 
# #   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size = 2)) + 
# #   rremove("legend")
# # 
# # p2 <- ggpubr::ggline(d1_thr, x = "m", y = "diff", group  = "Assemblage",color = "Assemblage",point.size = 0.2) + 
# #   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size = 2)) + 
# #   rremove("legend")
# # 
# # grid.arrange(p1,p2)
# 
# 
#   #Getting mininum M (sequences)
# 
# # }
# #  d1_thr <- do.call(rbind,d1_thr)
# #  #Obtaining probable threshold using first decile
# #  thr_final <- quantile(d1_thr$m,probs = 0.1)
# #  
# #  plot(d1_thr$m,d1_thr$seqs)
# # 
# # #Using estimateD
# #  
# # 
# # 
# # #############################################################################
# # 
# # ################################################################################
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # # #Obtaining Library sizes per samples
# # # sum_taxa_samp <- sample_sums(ps_rare)
# # 
# # 
# # 
# # #Testing if different samples affect richness results
# # 
# # # #for(i in 1:length(seqs_to)){
# # # for(i in 1:3){
# # #
# # #   ps_rare <- phyloseq::rarefy_even_depth(physeq4,
# # #                                          #sample.size = seqs_to[[i]],  # Or choose a fixed number
# # #                                          sample.size = sum_taxa_samp2[[i]],  # Or choose a fixed number
# # #                                          rngseed = 123,                        # For reproducibility
# # #                                          replace = FALSE,                      # No replacement
# # #                                          trimOTUs = TRUE,                      # Remove unobserved taxa
# # #                                          verbose = FALSE)
# # #   #
# # #
# # #   Alpha_i <- phyloseq::estimate_richness(ps_rare,measures = c("Observed","InvSimpson"))
# # #   Alpha_i$library_size <-sample_sums(ps_rare)
# # #   #Alpha_i$simulated <- as.character(seqs_to[[i]])
# # #   Alpha_i$simulated <- as.character(sum_taxa_samp2[[i]])
# # #   Alpha_i$nsamples <- nsamples(ps_rare)
# # #   Alpha_i$id <- row.names(Alpha_i)
# # #   Alph_list[[i]] <- Alpha_i
# # # }
# # # #ADDING NO RARIFIED
# # # Alph_list[[i+1]]<- data.frame(Observed=Alpha[,c("Observed")],
# # #                               InvSimpson=Alpha[,c("InvSimpson")],
# # #                               library_size  = as.numeric(sample_sums(physeq4)),
# # #                               simulated="No rarefaction",
# # #                               nsamples= nsamples(physeq4),
# # #                               id=row.names(Alpha)
# # # )
# # # rm(ps_rare)
# # # Alph_list <- do.call(rbind,Alph_list)
# # # Alph_list$simulated <- as.character(Alph_list$simulated)
# # # Alph_list$Fungi_system <-  dplyr::left_join(x = samp_data,y = Alph_list,by = c("id"="id"))$Fungi_system
# # # Alph_list$simulated <- factor(Alph_list$simulated,levels = c(as.character(sum_taxa_samp2[1:3]),"No rarefaction"))
# # # # ################################################################################
# # # # ####Plotting trajectories with Observed species
# # # label_data <- Alph_list %>%
# # #   group_by(id) %>%
# # #   filter(simulated == "No rarefaction") %>%
# # #   ungroup()
# # #
# # # p <- ggscatter(
# # #   data=Alph_list,
# # #   x="simulated",
# # #   y="Observed",#"Height",
# # #   combine = FALSE,
# # #   merge = FALSE,
# # #   color = "Fungi_system",
# # #   fill = "Fungi_system",
# # #   palette = "jco",
# # #   shape = 19,
# # #   size = 2,
# # #   point = TRUE,
# # #   rug = FALSE,
# # #   title = NULL,
# # #   xlab = "Reads",
# # #   ylab = "Observed species",
# # #   facet.by = "Fungi_system",
# # #   panel.labs = NULL,
# # #   short.panel.labs = TRUE,
# # #   add = c("reg.line"),
# # #   conf.int = TRUE,
# # #   conf.int.level = 0.95,
# # #   fullrange = FALSE,
# #   font.label = c(12, "plain"),
# #   font.family = "",
# #   label.select = NULL,
# #   repel = T,
# #   label.rectangle = FALSE,
# #   parse = FALSE,
# #   cor.coef = FALSE,
# #   cor.coeff.args = list(),
# #   cor.method = "pearson",
# #   cor.coef.coord = c(NULL, NULL),
# #   cor.coef.size = 4,
# #   ggp = NULL,
# #   show.legend.text = NA,
# #   ggtheme = theme_pubr()
# # ) +
# #   stat_cor(aes(color = Fungi_system)) +
# #   geom_line(aes(group = id)) +
# #   geom_text_repel(data = label_data,
# #                   aes(x = simulated, y = Observed, label = id, color = Fungi_system),
# #                   hjust = -0.1,
# #                   size = 3.2,
# #                   fontface = "italic",
# #                   show.legend = FALSE,max.overlaps = 10000)
# #
# # p <- ggpar(p, x.text.angle = 90)
# # p
# ################################################################################
# 
# ################################################################################
# 
# 
# 
# # sum_taxa_samp2 <- sum_taxa_samp[order(sum_taxa_samp,decreasing = F)]
# # #Trying some values
# # # seqs_to <- seq(100,20000,500)
# # ################################################################################
# # #Rarefaction curve
# # ps_rare <- ranacapa::ggrare(physeq4, step = 100, se = TRUE,parallel = T,color = "Sample",label = "Sample")
# # ps_rare  <- ps_rare + theme(legend.position = "none")
# # ps_rare <- ps_rare +
# #   labs(
# #     title = "Rarefaction Curves",
# #     x = "Number of Sequences",
# #     y = "Richness"
# #   )
# # ggsave(
# #   paste0(graph_dir, "/", "001_rarefaction_curve.pdf"), ps_rare,scale = 0.9,
# #   width = 25, height = 14, units = "in")
# 
# 
# 
# # physeq
# # plot_accum_sp <- function(physeq){
# # # Extract OTU table (taxa as columns)
# # otu_mat <- as(otu_table(physeq), "matrix")
# #
# # # Transpose if taxa_are_rows = TRUE
# # if (taxa_are_rows(physeq)) {
# #   otu_mat <- t(otu_mat)
# # }
# # tot <- rowSums(otu_mat)
# # S <- rowSums(otu_mat > 0)
# # nr <- nrow(otu_mat)
# #
# #
# # spec_accum <- vegan::specaccum(otu_mat, method = "exact", permutations = 1000)
# # plot(spec_accum)
# #   accum_df <- with(spec_accum, data.frame(Sites = sites,
# #                                         Samples = phyloseq::sample_names(physeq),
# #                                         NSeqs = tot,
# #                                         Richness = richness,
# #                                         SD = sd,
# #                                         Upper = richness + sd,
# #                                         Lower = richness - sd))
# #
# # accum_df <- accum_df[order(accum_df$Richness,decreasing = T),]
# # p <- ggpubr::ggline(accum_df,
# #                     x = "NSeqs",
# #                     y = "Richness",
# #                     color = "Samples",
# #                     xlab = "Number of sequences",ylab = "Richness") +
# #   geom_line(size = 1.2) +
# #   geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.4) +
# #   labs(title = "Species Accumulation Curve",
# #        x = "Sample size",
# #        y = "Observed Richness") +
# #   theme_minimal()
# # #return(p)
# # p
# # }