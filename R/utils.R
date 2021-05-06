
#'Helper function to get data set names from the path of csv file
#'
#'returns filenames with .csv ending
#'
#'@param path_list List of paths to csv files (one for each dataset)
#'@param string The parts of the file name that doesn't belong to the name, e.g. "abundance_cleaned"
#'
#'@return a list of data set file names
#'
#'@export

get_names <- function(path_list, string){
  n = length(path_list)
  names <- unlist(strsplit(path_list, string))[seq(2,n*2,2)]
  return(names)
}

#' Helper functions to get data set names 
#' 
#' returns filenames without .csv ending
#' 
#' @param path_list List of paths to csv files 
#' @param string The parts of the file name that doesn't belong to the name, e.g. "Abundance_cleaned"
#' 
#' @return a list of data set names
#' 
#' @export
#' 
get_names_sans <- function(path_list, string){
  n = length(path_list)
  names <- unlist(strsplit(path_list, string))[seq(2,n*2,2)]
  names_sans <- tools::file_path_sans_ext(names)
  return(names_sans)
}