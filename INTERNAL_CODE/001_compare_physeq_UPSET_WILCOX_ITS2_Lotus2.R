

wilcoxon_summary <- function(markers,dir, out_dir_name){
#MAKE DIR
out_comparison_dir <- paste0(dir,"/",out_dir_name)
if(!dir.exists(out_comparison_dir)){
  dir.create(out_comparison_dir)
}

#MAKE graphic dir
graph_dir <- paste0(out_comparison_dir,"/","graphics")
if(!dir.exists(graph_dir)){
  dir.create(graph_dir)
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

methods <- c("counts")
# levels <-c("Class","Order","Family", "Genus", "Species")
levels <-c("Genus", "Species")
#levels
x <- lapply(1:length(levels),function(i){
   # i <- 1
  #list files Wilcoxon 
  x <- list.files(path = csv_dir, pattern = "WILCOXON_")
  #obtain wilcoxon files per taxonomic level
  x_files <- grep(levels[[i]], x, value = TRUE)
  
  #filtered by FDR <0.05
  y <- lapply(1:length(x_files),function(j){
    y <- read.csv(paste0(csv_dir,"/",x_files[[j]]),header = T)
    # y$dataset <- x_files[[j]]
    y$context <- NA
    y$context[which(y$FDR<0.05)] <- 1 
    y$context[which(y$FDR>=0.05)] <- 0
    y$context[which(is.na(y$FDR))] <- 0
    return(y)
  })
  y <- do.call(rbind,y)
  y_only <- y[which(y$context==1),]
  y_only <-y_only[which(y_only$method!="CLR"),]
  # y_only <- y_only[order(y_only$FDR,decreasing = F),]
  # tapply(y_only$dataset,y_only$dataset,length)
  # y_only_s <- y_only %>%   count(taxon, name = "n") 
  return(y_only)
})

#creating a unique file
x <- do.call(rbind,x)
x$context <- NULL
#Obtaining what dataset has greater values for median
x$major_dataset_median <- NA
x$major_dataset_median[which(x$median_1 > x$median_2)] <- markers[[1]]
x$major_dataset_median[which(x$median_1 < x$median_2)] <- markers[[2]]

#Obtaining what dataset has greater values for mean
x$major_dataset_mean <- NA
x$major_dataset_mean[which(x$mean_1 > x$mean_2)] <- markers[[1]]
x$major_dataset_mean[which(x$mean_1 < x$mean_2)] <- markers[[2]]

#Creating subsets
# x_COUNT <- x[which(x$method=="COUNT"),]
x_RA <- x[which(x$method=="COUNT"),]

# #Saving results
# write.csv(x_COUNT,paste0(out_comparison_dir,"/","Wilcoxon_summary_COUNT.csv"),row.names=F)
write.csv(x_RA,paste0(out_comparison_dir,"/","Wilcoxon_summary_COUNT.csv"),row.names=F)
}


#Home folder
dir <- "~"

# #markers
# markers1 <- c(
#   "ITS2_chimecheckseparate.lotus2tax",
#   "ITS2_chimechecktogether.lotus2tax"
# )
# 
# markers2 <- c(
#   "ITS2_chimechecktogether.lotus2tax",
#   "ITS2_chimechecktogether.derepmin.lotus2tax"
# )
# 
# 
# markers3 <- c(
#   "ITS2_1st_batch",
#   "ITS2_chimechecktogether.derepmin.lotus2tax"
# )

# markers4 <- c(
#   "ITS2_1st_batch",
#   "chimechecktogether.derepmin.withdiscardsmapped.lotus2tax"
# )


# markers5 <- c(
#   "ITS2_1st_batch",
#   "ITS2lulu"
# )

#how to name the results
markers6 <- c(
  "ITS_ILLUMINA",
  "ITS_PACBIO"
)
# out_dir_name5 <- "comp_ITS2_20250725_FILT"out_dir_name1 <- "comparisons_lotus2"
# out_dir_name2 <- "comp_lotus_TOG_DEP"
# out_dir_name3 <- "comp_ITS2_LOT_DEREP_0714"
# out_dir_name4 <- "comp_ITS2_20250725"
# out_dir_name5 <- "comp_ITS2_20250725_FILT"
# out_dir_name5 <- "comparisons_ITS2_LULU"
out_dir_name6 <- "comparisons_ITS2_PACBIO_ILLUMINA_KAZ"
 # wilcoxon_summary(markers1,dir, out_dir_name1)
# wilcoxon_summary(markers2,dir, out_dir_name2)
# wilcoxon_summary(markers3,dir, out_dir_name3)
# wilcoxon_summary(markers4,dir, out_dir_name4)
# wilcoxon_summary(markers4,dir, out_dir_name5)
# wilcoxon_summary(markers5,dir, out_dir_name5)
wilcoxon_summary(markers = markers6,dir = dir, out_dir_name = out_dir_name6)
