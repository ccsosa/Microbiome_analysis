


#MAKE DIR
creat_dir_function <- function(dir,comparison_name){
  out_comparison_dir <- paste0(dir,"/",comparison_name)
  if(!dir.exists(out_comparison_dir)){
    dir.create(out_comparison_dir)
  }
  
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
  
  return(list(
    out_comparison_dir = out_comparison_dir,
    graph_dir = graph_dir,
    dens_dir = dens_dir,
    csv_dir = csv_dir,
    RDS_dir = RDS_dir
  ))
}