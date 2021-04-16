#Functions to calculate loop strengths and weights


#' Calculates the loop strength and weight of a given feedback loop
#' 
#' A loop is defined by a list of species e.g. c(1,2,3) is the 3-link loop 1 -> 2 -> 3
#' these species numbers are used as indices in order to pick the corresponding 
#' interaction strengths from the Jacobian A
#' 
#' @param loop a vector containing the species that are in the loop, e.g. c(1,2,3) 
#' for the 3-link loop 1 -> 2 -> 3
#' @param A a Jacobian matrix
#' 
#' @return a vector containing the loop strength and the loop weight 
#' 
loop_weight <- function(loop, A){
  #returns the strength and weights of a single loops (given by the species sequence loop)
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
  strength = abs(prod(coefficients))
  weight = strength^(1/l)
  
  return(c(strength, weight))
}

#'Calculates the strengths and weights of all loops of length n within the matrix A 
#'
#'To do this, I simply assume that all possible loops of length n (all possible combinations of n-species) exist. 
#'Then, the function goes through all possible loops and calculates their strengths using loop_weights()
#'If the loop weight is 0, at least one of the links is missing and the loop is removed from the list.
#'
#'@param n loop length
#'@param A a Jacobian matrix in which all loops of length should be identified
#'
#'@return A list containing the strengths and weights of all loops of length n within A 
#'
loops <- function(n, A){
  
  N <- nrow(A)  #number of species in the community 
  
  #combn returns a list of all unique n-species subsets of the network
  comb <- utils::combn(c(1:N), n, simplify = FALSE)
  
  #prepare lists to store loop weights and strengths 
  #(we know that there are twice as many possible loops as subsets in the list )
  #unless its 2-link loops -> maybe put in an if-condition here?
  loop_number = length(comb)*2
  strengths = vector('numeric', length = loop_number)
  weights = vector('numeric', length = loop_number)
  
  index = 1 #used to index l3_strengths and l3_weights
  
  for (i in comb){
    #first loop -> right order
    l1 <- loop_weight(i, A)
    
    #store the values in the lists prepared above
    strengths[index] <- l1[1]
    weights[index] <- l1[2]
    index = index + 1
    
    #same for the second loop -> same species but in reversed order
    l2 <- loop_weight(rev(i), A)
    strengths[index] <- l2[1]
    weights[index] <- l2[2]
    
    index = index+1
  }
  
  #loops that have a weight of zero do not exist because one(or more) of their links is 0
  #remove them from the list
  strengths <-strengths[strengths > 0]
  weights <- weights[weights > 0]
  
  #return the list of loop strengths and weights 
  return(list(strengths, weights))
}