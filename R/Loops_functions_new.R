#Functions to calculate loop strengths and weights


#' Calculate the loop strength and weight of a given feedback loop
#' 
#' Used internally by loops(). 
#' 
#' A loop is defined by a list of species e.g. c(1,2,3) is the 3-link loop 1 -> 2 -> 3
#' these species numbers are used as indices in order to pick the corresponding 
#' interaction strengths from the Jacobian \eqn{\mathbf{A}}
#' 
#' @param loop a vector containing the species that are in the loop, e.g. c(1,2,3) 
#' for the 3-link loop 1 -> 2 -> 3
#' @param A A Jacobian matrix
#' 
#' @return A vector containing the loop strength and the loop weight 
#' 
#' 
loop_weight <- function(loop, A){
  #returns the strength and weights of a single loop
  #a loop is defines by the sequence of species
  
  l = length(loop) #length of the loop
  
  #select the right coefficients from the matrix A 
  #first species -> index of the row we need 
  #second species -> index of the column we need 
  
  coefficients <- vector('numeric', length = l)
  
  row_indices <- loop
  col_indices <- c(loop[2:l], loop[1])
  
  for (i in 1:l){
    coefficients[i] <- A[row_indices[i], col_indices[i]]
  }
  
  #calculate loop strength and loop weight
  #loop strength = absolute product of all coefficients 
  strength = prod(coefficients)
  weight = abs(strength)^(1/l) * sign(strength)
  
  return(data.frame(species_sequence = paste(as.character(loop), collapse = " "),
                    loop_strength = strength, loop_weight= weight))
}


#' Calculate the strengths and weights of all loops of length *n* within the matrix *A* 
#'
#' To do this, I simply assume that all possible loops of length n (all possible 
#' combinations of n-species) exist. Then, the function goes through all possible
#' loops and calculates their strengths using loop_weights().
#' 
#' If the loop weight is 0, at least one of the links is missing and the loop
#' is removed from the list.
#'
#'@param n loop length
#'@param A a Jacobian matrix in which all loops of length should be identified
#'
#'@return A dataframe containing the loops, their strengths and weights
#'
#'@export
#'

loops <- function(n, A){
  
  #find all loops of length n in the matrix
  #and calculate their strengths as loop weights
  
  N <- nrow(A)  #number of species in the community 
  
  #combn returns a list of all unique n-species subsets of the network
  comb <- utils::combn(c(1:N), n, simplify = FALSE)
  
  #initialise a data frame to store the results for all subsets 
  #r <- data.frame(species_sequence = character(),
  #                loop_strength = numeric(), 
  #                loop_weight = numeric())
  datalist <- list()
  
  #loop through all n-species subsets and calculate loop weights
  for (i in seq_along(comb)){
    
    subset_i <- comb[[i]] #pick one n-species subset 
    start_species <- subset_i[1] #what is the first species in the list?
    
    all_loops <- gtools::permutations(n,n,subset_i)
    
    #only those permutations that begin with the start_species are unique feedback loops
    all_loops <- all_loops[all_loops[,1]==start_species,]
    
    #for 2-link loops, there is only one loop per subset
      if(n == 2){
        all_weights <- loop_weight(all_loops, A)
      } else {
      
        #for longer loops there are several loops per subset
        #purrr is used to apply the loop weights function to each of them 
      
        #turn all loops into a list of vectors first
        all_loops_list <- as.list(as.data.frame(t(all_loops)))
      
        #apply loop weights function to each element of the list
        all_weights <- purrr::map_dfr(all_loops_list, loop_weight, A)
      }
    
    datalist[[i]] <- all_weights
    #add loop weights to the data frame 
    #r <- rbind(r, all_weights)
  }
  
  r <- dplyr::bind_rows(datalist)
  #remove zeros from the list 
  r <- r[r$loop_strength != 0,]
  
  return(r)
}


#'Calculates the strengths and weights of all loops of length *n* within the matrix *A*
#'
#'Same function as loops() but calculations can be run on multiple cores
#'(setup of parallel backend only works on Linux!!!)
#'
#'@param n loop length
#'@param A a Jacobian matrix in which all loops of length should be identified
#'@param n_cores number of cores used for computation
#'
#'@return A dataframe containing the loops, their strengths and weights
#'
#'@importFrom foreach %dopar%
#'
#'@export
#'

loops_parallel <- function(n, A, n_cores){
  
  i <- NULL #bind variable to avoid warnings 
  
  #set up parallel backend
  doMC::registerDoMC(n_cores)
  
  #find all loops of length n in the matrix
  #and calculate their strengths as loop weights
  
  N <- nrow(A)  #number of species in the community 
  
  #combn returns a list of all unique n-species subsets of the network
  comb <- utils::combn(c(1:N), n, simplify = FALSE)
  
  datalist <- list()
  
  #loop through all n-species subsets and calculate loop weights
  datalist <- foreach::foreach (i = seq_along(comb)) %dopar% {
    
    subset_i <- comb[[i]] #pick one n-species subset 
    start_species <- subset_i[1] #what is the first species in the list?
    
    all_loops <- gtools::permutations(n,n,subset_i)
    
    #only those permutations that begin with the start_species are unique feedback loops
    all_loops <- all_loops[all_loops[,1]==start_species,]
    
    #for 2-link loops, there is only one loop per subset
    if(n == 2){
      all_weights <- loop_weight(all_loops, A)
    } else {
      
      #for longer loops there are several loops per subset
      #purrr is used to apply the loop weights function to each of them 
      
      #turn all loops into a list of vectors first
      all_loops_list <- as.list(as.data.frame(t(all_loops)))
      
      #apply loop weights function to each element of the list
      all_weights <- purrr::map_dfr(all_loops_list, loop_weight, A)
    }
    
  }
  
  r <- dplyr::bind_rows(datalist)
  #remove zeros from the list 
  r <- r[r$loop_strength != 0,]
  
  return(r)
}



