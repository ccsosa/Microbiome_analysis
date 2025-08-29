
wilcoxon_function_microb <- function(tax_to_test_original,x_list,sam_to_test){


#get wilcoxon paired tests
wilc_test <- lapply(1:length(tax_to_test_original),function(m){
     # print(m)
      # m <- 14
  x <- lapply(1:length(x_list),function(n){
      # print(n)
     # n <- 2
    # if(tax_to_test_original[[m]] %in% colnames(x_list[[n]])==FALSE){
    #   x <- data.frame(value = NA,
    #                   marker = markers[[n]])
    #   
    # } else {
      x <- data.frame(value = as.numeric(x_list[[n]][tax_to_test_original[[m]],sam_to_test]),
                      marker = markers_name[[n]])
    #   
    # }
    return(x)
  })
  # if(sum(x[[1]][,1],na.rm = T)>0  & sum(x[[1]][,1],na.rm = T)>0){
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
  
  # }else {
  #   x2 <- data.frame(taxon=tax_to_test_original[[m]],
  #                    p.value=NA,
  #                    ks=NA,
  #                    median_1 = NA,
  #                    mad_1 = NA,
  #                    mean_1 = NA,
  #                    sd_1 = NA,
  #                    median_2 = NA,
  #                    mad_2 = NA,
  #                    mean_2 = NA,
  #                    sd_2 = NA
  #   )
  #   
  # }
  # colnames(x2) <- c("taxon","p.value")
  return(x2)
})
return(wilc_test)
}

library(ggpubr)

plot_density_function <- function(dens_dir,wilc_test,markers_name,x_list,sam_to_test,method,level,dataset){
  message("detecting taxa to show density differences")
  #taxa to use for density plots
  taxa_to_test_density <- wilc_test$taxon
  
  #Obtaining data from two datasetss to create density plots
  x <- lapply(1:length(taxa_to_test_density),function(m){
    # m <- 1
    y <- lapply(1:length(markers_name),function(n){
      x <- data.frame(value = as.numeric(x_list[[n]][taxa_to_test_density[[m]],sam_to_test]),
                      marker = markers_name[[n]])
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
      filename = paste0(dens_dir,"/",level,"_",taxa_to_test_density[[m]],"_",dataset,"_",method, ".pdf"),
      plot = x_dens,
      width = 15,
      height = 10,
      units = "in",
      dpi = 600,
      device = cairo_pdf
    )
  })
}
