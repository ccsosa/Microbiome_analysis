################################################################################
remove_outliers <- function(vec, na.rm = TRUE) {
  qnt <- quantile(vec, probs = c(.25, .75), na.rm = na.rm)
  H <- 1.5 * IQR(vec, na.rm = na.rm)
  vec[vec < (qnt[1] - H) | vec > (qnt[2] + H)] <- NA
  return(vec)
}
################################################################################


cor_alpha_function <- function(data,indexes,vars,padj_threshold){
  #list for correaltion per index
  cor_list_final <- list()
  x <- data[,vars]

  #performing correlation per Alpha indexes  
  for(i in 1:length(indexes)){
    print(i)
    #subsetting to index
    index <- indexes[[i]]
    cor_list <- list()
    #loop for correlation
    for(j in 1:ncol(x)){
      # print(j)
      # j <- 1
      #subsetting to index and soil variable
      xy_i <- data.frame(index = index,
                         val = data[,index],
                         soil_var = colnames(x)[j],
                         soil_var_val = x[,j]
      )
      #Only using complete cases
      xy_i <- xy_i[complete.cases(xy_i),]
      #remove outliers
      xy_i <- xy_i[xy_i$soil_var_val %in% na.omit(remove_outliers(xy_i$soil_var_val)),]
      #adding range of values used in the correlation
      range_x <- range(xy_i$val)
      range_y <- range(xy_i$soil_var_val)
      #using only 15 or more full observations
      if(nrow(xy_i)>14){
        x_cor <- cor.test(xy_i$val,xy_i$soil_var_val,method="spearman")
        x_cor <- data.frame(
          index = index,
          soil_var = colnames(x)[j],
          range_index = paste(range_x,collapse = " - "),
          range_soil_var = paste(range_y,collapse = " - "),
          correlation =x_cor$estimate,
          n = nrow(xy_i),
          p.value = x_cor$p.value,
          p.adj = NA
        )
      } else {
        
        x_cor <- data.frame(
          index = index,
          soil_var = colnames(x)[j],
          range_index = NA,
          range_soil_var = NA,
          correlation =NA,
          n = NA,
          p.value = NA,
          p.adj = NA
        )
      }
      
      rm(xy_i)
      
      cor_list[[j]] <- x_cor
    };rm(j)
    
    cor_list <- do.call(rbind,cor_list)
    cor_list$p.adj <- p.adjust(cor_list$p.value,
                               method = "fdr"
    )
    
    cor_list_final[[i]] <- cor_list
  };rm(i)
  
  #final correlation
  cor_list_final <- do.call(rbind,cor_list_final)
  #removing p value < padj_threshold
  cor_list_final <- cor_list_final[cor_list_final$p.adj<padj_threshold,]
  return(cor_list_final)
}

