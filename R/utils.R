
#define a function to get the names from one of the path lists
get_names <- function(path_list, string){
  n = length(path_list)
  names <- unlist(strsplit(path_list, string))[seq(2,n*2,2)]
  names_sans <- file_path_sans_ext(names)
  return(list(names, names_sans))
}