#'Calculates total feedback at a given level *k*
#'
#'F_k is calculated via the coefficients of the characteristic polynomial 
#'
#'@param k level for which total feedback should be calculated
#'@param A a Jacobian matrix 
#'
#'@return F_k total feedback at level k 
#'
#'@export
#'

get_Fk <- function(k, A){
  
  #get the coefficients of the characteristic polynomial using the charpoly function
  charpoly <- pracma::charpoly(A)
  
  #we have to get the (k+1)th element since there is also F0, which is not 
  #interesting for us
  F_k <- charpoly[(k+1)]*-1
  
  return(F_k)
}