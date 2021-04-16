
#'Helper function to get data set names from the path of csv files
#'
#'@param path_list List of paths to csv files (one for each dataset)
#'@param string The parts of the file name that doesn't belong to the name, e.g. "abundance_cleaned"
#'
#'@export

get_names <- function(path_list, string){
  #split into two function? -> get_names and get_names_sans
  n = length(path_list)
  names <- unlist(strsplit(path_list, string))[seq(2,n*2,2)]
  names_sans <- tools::file_path_sans_ext(names)
  return(list(names, names_sans))
}