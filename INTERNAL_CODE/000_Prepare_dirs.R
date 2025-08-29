#' Initialize Directory Structure for Analysis Project
#'
#' Creates a standard set of directories under a given base directory:
#' - `SCRIPTS/` for R scripts
#' - `outcomes/` with subfolders for `graphics/`, `csv_dir/`, and `rds_dir/`
#'
#' This function is useful for organizing project outputs and code files.
#'
#' @param dir Character string. Path to the base directory. Defaults to the user's home directory `"~"`.
#'
#' @return A named list containing the full paths of all created (or existing) directories:
#' `SCRIPTS_dir`, `out_dir`, `graph_dir`, `csv_dir`, `RDS_dir`.
#' @export
#'
#' @examples
#' \dontrun{
#' paths <- init_directories("~")
#' print(paths$csv_dir)
#' }

dir_function <- function(dir,marker){
  # Load or create main directory
  message("Loading or creating scripts dir")
  
  SCRIPTS_dir <- paste0(dir, "/SCRIPTS")
  if (!dir.exists(SCRIPTS_dir)) dir.create(SCRIPTS_dir)
  
  marker_dir <- marker
  ("Loading or creating marker dir")
  marker_dir <- paste0(dir, "/",marker_dir)
  if (!dir.exists(marker_dir)) dir.create(marker_dir)
  
  
  # Create Outcomes dir
  message("Loading or creating outcomes dir")
  out_dir <- paste0(marker_dir, "/outcomes")
  if (!dir.exists(out_dir)) dir.create(out_dir)
  
  # Create subfolders inside Outcomes
  message("Loading or creating outcomes subfolders")
  graph_dir <- paste0(out_dir, "/graphics")
  if (!dir.exists(graph_dir)) dir.create(graph_dir)
  
  csv_dir <- paste0(out_dir, "/csv_dir")
  if (!dir.exists(csv_dir)) dir.create(csv_dir)
  
  RDS_dir <- paste0(out_dir, "/rds_dir")
  if (!dir.exists(RDS_dir)) dir.create(RDS_dir)
  
  # # Return paths as a named list
  return(list(
    SCRIPTS_dir = SCRIPTS_dir,
    out_dir = out_dir,
    graph_dir = graph_dir,
    csv_dir = csv_dir,
    RDS_dir = RDS_dir
  ))

}
 