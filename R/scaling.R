#Functions used to apply Anjes Scaling procedure 

usethis::use_package("dplyr")
#see: Neutel and Thorne, 2016 + 2021

#'Replace missing diagonal values 
#'
#'The scaling procedure as well as the calculation of s* require that all 
#'diagonal values $a_{ii} != 0$. Therefore, missing values need to be replaced by
#'some small value.
#'
#'This replacement value is calculated as the mean strengths of all off-diagonal 
#'matrix elements, mutltiplied with a factor f.
#'
#'@param Jac A Jacobian matrix
#'@param f factor multiplied with mean interaction strength to choose replacement value
#'
#'@return A new version of the Jacobian in which missing diagonal values have been replaced
#'
#'@export
#'

replace_zeros <- function(Jac, f){
  
  #calculate mean of off-diagonal interaction strengths (that are not zero)
  Jac_copy = Jac
  diag(Jac_copy) = 0
  mean_aij = mean(Jac_copy[Jac_copy !=0])
  
  #replace zeros in the Jacobian's diagonal with mean of interaction strengths * factor
  replacement <- mean_aij * f
  
  original_diagonal <- diag(Jac) #get original diagonal 
  zeros <- original_diagonal == 0 #records the indices of zeros
  
  new_diagonal <- replace(original_diagonal, zeros, replacement)#replace zeros
  
  diag(Jac) <- new_diagonal 
  
  return(Jac)
}


#' Scale interaction strengths within a Jacobian matrix
#'
#' Off-diagonal values are scaled by dividing them by the absolute value of 
#' their corresponding diagonal element. The resulting scaled interaction strengths
#' are thus defined as: \eqn{\bar{a_{ij}} = \frac{a_{ij}}{|a_{ii}|}} (Neutel and Thorne
#' 2016).
#' 
#' This does not work if there is a 0 on the diagonal ! Replace missing elements
#' with replace_zeros() before scaling !. 
#'
#'@param A A Jacobian matrix, with non-zero diagonal elements
#'@param aii determines the value of diagonal elements 
#'
#'@return The same Jacobian matrix with scaled interaction strengths and -1 on the diagonal 
#'
#'@export
#'

scale_matrix <- function(A, aii = 0){
  
  #extract diagonal values
  diagonal <- diag(A)
  #get inverse of those 
  d_inv <- -1/diagonal
  
  #divide matrix elements with corresponding diagonal values:
  A_scaled <- A * d_inv
  
  #set diagonal to 0 
  diag(A_scaled) <- aii
  
  return(A_scaled)
}

#' Scale interaction strengths within an interaction table
#'
#' Off-diagonal values are scaled by dividing them by the absolute value of 
#' their corresponding diagonal element. The resulting scaled interaction strengths
#' are thus defined as: \eqn{\bar{a_{ij}} = \frac{a_{ij}}{|a_{ii}|}} (Neutel and Thorne
#' 2016).  
#' 
#' For each pairwise interaction, the function identifies the corresponding 
#' self-regulation term needed to scale an interspecific interaction strength.
#' If the self-regulation value is missing, it is replaced by the mean 
#' interaction strength multiplied  with a specified factor. 
#' 
#' @param interaction_table table of pairwise interactions
#' @param r_factor factor that is mulitplied with the mean interaction strength 
#' to replace missing self-regulation values
#' 
#' @return interaction table data.frame with two new columns a_ij_scaled and a_ji_scaled
#' 
#' @export
#' @importFrom rlang .data
#' 
scale_interaction_table <- function(interaction_table, r_factor){
  
  Species_i <- Species_j <- . <- NULL
  
  #missing diagonal values are replaced by mean_aij* replacement
  
  #calculate mean interspecific interaction strengths
  mean_aij <- interaction_table %>% 
    dplyr::filter(.data$Species_i != .data$Species_j) %>% 
    dplyr::select(.data$a_ij, .data$a_ji) %>% 
    unlist() %>% mean()
  
  a_ij_scaled <- vector("numeric", length = nrow(interaction_table))
  a_ji_scaled <- vector("numeric", length = nrow(interaction_table))
  
  #scale interactions strengths in interaction table
  for (i in 1:nrow(interaction_table)){
    
    species_i <- interaction_table[i,]$Species_i
    species_j <- interaction_table[i,]$Species_j
    
    #only continue for interspecific interactions
    
    if (species_i != species_j){
      
      #pick a_ii
      a_ii <- interaction_table %>% 
        dplyr::filter(Species_i == species_i & Species_j == species_i) %>% .$a_ii
      
      a_jj <- interaction_table %>% 
        dplyr::filter(Species_i == species_j & Species_j == species_j) %>% .$a_ii
      
      
      #check whether a_ii or a_jj are missing, 
      #set them to mean_aij * replacement factor
      if (length(a_ii) == 0){a_ii <- mean_aij * r_factor}
      if (length(a_jj) == 0){a_jj <- mean_aij * r_factor}
      
      #scale F_ij by dividing it by a_ii
      a_ij_scaled[i] <- (interaction_table[i,]$a_ij / a_ii)*-1
      
      #scale F_ji by dividing it by a_jj
      a_ji_scaled[i] <- (interaction_table[i,]$a_ji/ a_jj)*-1
    }
  }
  
  interaction_table$a_ij_scaled <- a_ij_scaled
  interaction_table$a_ji_scaled <- a_ji_scaled
  
  return(interaction_table)
}