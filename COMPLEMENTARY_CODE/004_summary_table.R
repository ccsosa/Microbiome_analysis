library(ggrepel)
library(phyloseq)
library(iNEXT)
library(ggpubr)



summary_table <- function(markers_name,physeq,physeq_original,physeq_list,name){
  #SUMMARY TABLE
  message("creating summary, be patient")
  summary <- list()
  #Loading Alpha diversity results
  Alpha_N <- read.csv(paste0(csv_dir,"/",name,"_002_Alpha.csv"))
  
  for(i in 1:length(markers_name)){
    #counts
    counts <- sample_sums(physeq_list[[i]])
    #loading joined phyloseq object
    ps_rare_joined <- subset_samples(physeq, file==markers_name[[i]])
    ps_rare_joined2 <- prune_taxa(taxa_sums(ps_rare_joined) > 0, ps_rare_joined)
    # message(paste0("processing ",markers_name[[i]]))
    
    orig_families = length(unique(unlist(lapply(1:nrow(tax_table(physeq_original[[i]])),function(j){
      x <- paste(tax_table(physeq_original[[i]])[j,1:5],collapse = "-")
      return(x)
    }))))
    orig_genus = length(unique(unlist(lapply(1:nrow(tax_table(physeq_original[[i]])),function(j){
      x <- paste(tax_table(physeq_original[[i]])[j,1:6],collapse = "-")
      return(x)
    }))))
    
    orig_sp = length(unique(unlist(lapply(1:nrow(tax_table(physeq_original[[i]])),function(j){
      x <- paste(tax_table(physeq_original[[i]])[j,1:7],collapse = "-")
      return(x)
    }))))
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
      })))))
    
    
    filtered_families =
      length(unique(unlist(lapply(1:nrow(tax_table(physeq_list[[i]])),function(j){
        x <- paste(tax_table(physeq_list[[i]])[j,1:5],collapse = "-")
        return(x)
      }))))
    
    filtered_genus =
      length(unique(unlist(lapply(1:nrow(tax_table(physeq_list[[i]])),function(j){
        x <- paste(tax_table(physeq_list[[i]])[j,1:6],collapse = "-")
        return(x)
      }))))
    
    filtered_sp =
      length(unique(unlist(lapply(1:nrow(tax_table(physeq_list[[i]])),function(j){
        x <- paste(tax_table(physeq_list[[i]])[j,1:7],collapse = "-")
        return(x)
      }))))
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
      })))))
    
    Alpha_N_i <- Alpha_N[Alpha_N$file==markers_name[[i]],]
    
    
    marker = markers_name[[i]]
    
    # message("original phyloseq objects without filtering summary metrics")
      orig_nsamples = length(sample_sums(physeq_original[[i]]))
      orig_read_counts = sum(sample_sums(physeq_original[[i]]))
      orig_min = min(sample_sums(physeq_original[[i]]),na.rm=T)
      orig_max =  max(sample_sums(physeq_original[[i]]),na.rm=T)
      orig_min_sample = names(sample_sums(physeq_original[[i]])[which.min(sample_sums(physeq_original[[i]]))])
      orig_max_sample = names(sample_sums(physeq_original[[i]])[which.max(sample_sums(physeq_original[[i]]))])

      # message("filtered phyloseq objects filtering summary metrics")
      filtered_nsamples = length(counts)
      filtered_read_counts = sum(sample_sums(physeq_list[[i]]))
      filtered_min = min(counts,na.rm=T)
      filtered_max =  max(counts,na.rm=T)
      filtered_min_sample = names(counts[which.min(counts)])
      filtered_max_sample = names(counts[which.max(counts)])
      filtered_families = filtered_families
      filtered_genus = filtered_genus
      filtered_sp = filtered_sp
      filtered_sp_no_na = filtered_sp_no_na
      
      
      #Alpha
      joined_obs_sp_av = mean(Alpha_N_i$Observed,na.rm=T)
      joined_obs_sp_sd = sd(Alpha_N_i$Observed,na.rm=T)
      joined_Shannon_av = mean(Alpha_N_i$Shannon,na.rm=T)
      joined_Shannon_sd = sd(Alpha_N_i$Shannon,na.rm=T)
      joined_InvSimpson_av = mean(Alpha_N_i$InvSimpson,na.rm=T)
      joined_InvSimpson_sd = sd(Alpha_N_i$InvSimpson,na.rm=T)
      joined_Pielou_av = mean(Alpha_N_i$pielou,na.rm=T)
      joined_Pielou_sd = sd(Alpha_N_i$pielou,na.rm=T)
      joined_Chao1_av = mean(Alpha_N_i$Chao1,na.rm = T)
      joined_Chao1_sd = sd(Alpha_N_i$Chao1,na.rm = T)
      
    x <- data.frame(
      #marker
      marker = marker,
      # message("original phyloseq objects without filtering summary metrics")
      orig_nsamples = orig_nsamples,
      orig_read_counts = orig_read_counts,
      orig_min = orig_min,
      orig_max =  orig_max,
      orig_min_sample = orig_min_sample,
      orig_max_sample = orig_max_sample,
      orig_families = orig_families,
      orig_genus = orig_genus,
      orig_sp = orig_sp,
      orig_sp_no_na = orig_sp_no_na,

      #filtered phyloseq objects summary metrics
      filtered_nsamples = filtered_nsamples,
      filtered_read_counts = filtered_read_counts,
      filtered_min = filtered_min,
      filtered_max =  filtered_max,
      filtered_min_sample = filtered_min_sample,
      filtered_max_sample = filtered_max_sample,
      filtered_families = filtered_families,
      filtered_genus = filtered_genus,
      filtered_sp = filtered_sp,
      filtered_sp_no_na = filtered_sp_no_na,

      joined_obs_sp_av = joined_obs_sp_av,
      joined_obs_sp_sd = joined_obs_sp_sd,
      joined_Shannon_av = joined_Shannon_av,
      joined_Shannon_sd = joined_Shannon_sd,
      joined_InvSimpson_av = joined_InvSimpson_av,
      joined_InvSimpson_sd = joined_InvSimpson_sd,
      joined_Pielou_av = joined_Pielou_av,
      joined_Pielou_sd = joined_Pielou_sd,
      joined_Chao1_av = joined_Chao1_av,
      joined_Chao1_sd = joined_Chao1_sd
      
    )
    summary[[i]] <- x
  }
  summary <- do.call(rbind,summary)
  write.csv(summary,paste0(csv_dir, "/",name,"_comparison_summary.csv"),row.names = F)
}
