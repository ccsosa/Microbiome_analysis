library(ggrepel)
library(phyloseq)
library(iNEXT)
library(ggpubr)

Alpha_function <- function(physeq,name){
  ##Getting outliers
  findoutlier <- function(x) {
    return(x < quantile(x, .25) - 1.5*IQR(x) | x > quantile(x, .75) + 1.5*IQR(x))
  }
  # Now add labels to your plot
  p_r_p4_rare <-
    plot_richness(
      physeq,
      x = "file",
      color = "file",
      measures = c("Observed","Shannon", "InvSimpson","Chao1"),
      title = ""
    )
  #getting data
  out_df <- p_r_p4_rare$data
  #obtaining Alpha measures to get outliers
  measures = c("Observed","Shannon", "InvSimpson","Chao1")
  
  #Getting outliers per Alpha measure
  i_list <- list()
  for(i in 1:length(markers_name)){
    #per sequencing tech
    j_list <- list()
    #getting outliers
    for(j in 1:length(measures)){
      x <- out_df$SampleID[which(out_df$variable==measures[[j]] & out_df$file==markers_name[[i]])]
      x_i <- findoutlier(
        out_df$value[which(out_df$variable==measures[[j]] & out_df$file==markers_name[[i]])]
      )
      #subsetting values
      x <- x[x_i]
      outlier <- out_df$value[which(out_df$variable==measures[[j]] & out_df$file==markers_name[[i]])][x_i]
      #saving only 
      if(length(x)>0){
        #saving outliers
        j_list[[j]] <- data.frame(SampleID = x,
                                  variable=measures[[j]] ,
                                  value=outlier,
                                  file=markers_name[[i]])
      }
    }
    #saving outliers
    j_list <- do.call(rbind,j_list)
    i_list[[i]] <- j_list
  }
  i_list <- do.call(rbind,i_list)
  
  #outlier data.frame
  out_df <- i_list
  
  p_r_p4_rare <- p_r_p4_rare +
    ylab("") +
    geom_boxplot(aes(fill = file), alpha = 0.7) +
    stat_compare_means(
      method = "wilcox.test",
      size = 8,
      aes(label = paste0("p = ", after_stat(p.format)))
    ) +
    geom_text_repel(
      data = out_df,
      aes(x = file, y = value, label = SampleID), # change `sample_id` to your ID variable
      size = 5,max.overlaps = 10000
    ) +
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
  
  #
  ggsave(
    paste0(graph_dir, "/", name,"_002_Alpha_boxplots.pdf"),
    p_r_p4_rare,
    width = 20,
    height = 25,
    units = "in",
    dpi = 600
  )
  
  ##############################################################################
  Alpha_N <- phyloseq::estimate_richness(physeq)
  Alpha_N$pielou <- evenness(physeq, index = "pielou")[, 1]
  calc <- ChaoRichness(x = otu_table(physeq), datatype = "abundance", conf = 0.95)
  Alpha_N$Chao1 <- calc$Estimator
  Alpha_N$se.chao1 <- calc$Est_s.e.
  
  #Adding sample data
  Alpha_N <- cbind(sample_data(physeq), Alpha_N)
  
  message("Saving Alpha diversity indexes CSV files")
  write.csv(
    Alpha_N,
    paste0(csv_dir, "/", name,"_002_Alpha.csv"),
    na = "",
    row.names = F,
    quote = T
  )
  
}


