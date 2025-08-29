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

################################################################################
remove_outliers <- function(vec, na.rm = TRUE) {
  qnt <- quantile(vec, probs = c(.25, .75), na.rm = na.rm)
  H <- 1.5 * IQR(vec, na.rm = na.rm)
  vec[vec < (qnt[1] - H) | vec > (qnt[2] + H)] <- NA
  return(vec)
}
################################################################################
cor_var_function <- function(ps,padj_threshold){
# Trying to get correlation per taxonomic level
  
levels_tax <- c("Genus","Species","OTU")
# i <- 1
#Correlation function
cor_levels <- lapply(1:length(levels_tax),function(i){

  message(paste0("correlation table for ",levels_tax[[i]]))
  if(levels_tax[[i]]=="OTU"){
    pF <- ps
    message("transforming using compositional approach")
    pF <- microbiome::transform(pF, "compositional")
    tax_table(pF)[,7] <- paste0(tax_table(pF)[,7],"-",row.names(tax_table(pF)))
    taxa_names(pF) <- tax_table(pF)[,7]   } else {
    #Aggregating taxa
    pF <- microViz::tax_agg(ps = ps, levels_tax[[i]])
    message("transforming using compositional approach")
    pF <- microbiome::transform(pF, "compositional")
  }
  
  #assigning ps to pF
  # ps <- pF
  meta_pF <- data.frame(meta(pF))
  #soil data
  sample_pF <- meta_pF[,c(34:45,47:73)]
  #removing zeros
  x = data.frame(t(otu_table(pF)))
  tol <- 1e-12
  x[abs(x) < tol] <- NA
  #Performing correlations with FDR <0.05

  cor_list <- list()
  for(j in 1:ncol(x)){
    # j <- 1
    xy_i_list <- list()
    for(k in 1:ncol(sample_pF)){
      # k <- 1
      xy_i <- data.frame(level = levels_tax[i],
                         taxa = colnames(x)[j],
                         taxa_freq = x[,j],
                         var= colnames(sample_pF)[k],
                         var_freq = sample_pF[,k])
      #Only using complete cases
      xy_i <- xy_i[complete.cases(xy_i),]
      #remove outliers
      xy_i <- xy_i[xy_i$var_freq %in% na.omit(remove_outliers(xy_i$var_freq)),]
      #adding range of values used in the correlation
      range_x <- range(xy_i$taxa_freq)
      range_y <- range(xy_i$var_freq)
      #using only 15 or more full observations
      if(nrow(xy_i)>14){
        x_cor <- cor.test(xy_i$taxa_freq,xy_i$var_freq,method="spearman")
        x_cor <- data.frame(level = levels_tax[i],
                            taxa = colnames(x)[j],
                            mean_taxa = mean(xy_i$taxa_freq),
                            sd_taxa = sd(xy_i$taxa_freq),
                            var = colnames(sample_pF)[k],
                            mean_var= mean(xy_i$var_freq),
                            sd_var = sd(xy_i$var_freq),
                            range_taxa = paste(range_x,collapse = " - "),
                            range_var = paste(range_y,collapse = " - "),
                            correlation =x_cor$estimate,
                            n = nrow(xy_i),
                            p.value = x_cor$p.value,
                            p.adj = NA
        )
      } else {
        x_cor <- data.frame(level = levels_tax[i],
                            taxa = colnames(x)[j],
                            mean_taxa =NA,
                            sd_taxa = NA,
                            var = colnames(sample_pF)[k],
                            mean_var = NA,
                            sd_var = NA,
                            range_taxa = NA,
                            range_var = NA,
                            correlation =NA,
                            n = NA,
                            p.value = NA,
                            p.adj = NA
        )
      }
      
      rm(xy_i)
      xy_i_list[[k]] <- x_cor
    }
    xy_i_list <- do.call(rbind,xy_i_list)
    if(length(na.omit(xy_i_list$p.value))>2){
    xy_i_list$p.adj <- p.adjust(xy_i_list$p.value,
                                method = "fdr"
    )
                                
    } else {
      xy_i_list$p.adj <- NA
    }
    cor_list[[j]] <- xy_i_list
  }
  correlation.table <- do.call(rbind,cor_list)

  # correlation.table <-
  #   microbiome::associate(
  #     x = x,
  #     y = sample_pF,
  #     method = "spearman",
  #     mode = "table",
  #     order = T,
  #     # p.adj.threshold = padj_threshold,
  #     n.signif = NULL,
  #     p.adj.method = "fdr",
  #     filter.self.correlations = T
  #   )
  #filtering at p-value desired
  correlation.table <- correlation.table[which(correlation.table$p.adj<padj_threshold),]

  #Returning only tables with correlations with FDR < 0.01
  if(nrow(correlation.table)>0){
    # colnames(correlation.table) <- c("id","var","correlation","p.adj")
    # #adding 
    # #tax table to be added to correlations
    # tax <- data.frame(tax_table(pF))
    # tax$id <- row.names(tax)
    # tax <- tax[,c(ncol(tax),1:(ncol(tax)-1))]
    # 
    # xx <-
    #   dplyr::left_join(x = correlation.table,
    #                    y = tax,
    #                    by = c("id" = "id"))
    # if(levels_tax[[i]]=="OTU"){
    #   xx <- data.frame(id = xx[,"id"],
    #                    var = xx[,"var"],
    #                    correlation = xx[,"correlation"],
    #                    p.adj = xx[,"p.adj"],
    #                    taxonomy = xx[,"id"]
    #                  )
    # 
    # } else {
    #   xx <- data.frame(id = xx[,"id"],
    #                    var = xx[,"var"],
    #                    correlation = xx[,"correlation"],
    #                    p.adj = xx[,"p.adj"],
    #                    taxonomy = xx[,levels_tax[[i]]]
    #   )
    # }
    # xx$level <- levels_tax[[i]]
    xx <- correlation.table
    return(xx)
  } else {
    message(paste0("No results for level: ",levels_tax[[i]]))
    return(NULL)
  }
})

cor_levels <- do.call(rbind,cor_levels)
return(cor_levels)
}
