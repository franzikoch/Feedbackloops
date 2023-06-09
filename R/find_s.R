usethis::use_pipe(export = TRUE)

#'Find the relative amount of self-regulation needed for stability s*
#'
#'For unstable matrices, diagonal values are incrementally increased by multiplying it 
#'with a factor *s*, until its dominant eigenvalue becomes negative (=the matrix becomes stable). 
#'
#'This will not work, if the diagonal contains zeros! 
#'Replace missing values before calculating s*!
#'
#'For stable matrices, diagonal values are incrementally decrease, until the 
#'dominant eigenvalue becomes positive (= the matrix is unstable)
#'
#'@param Jacobian A Jacobian matrix
#'@param step_size Amount by which *s* is increases in each step 
#'@param start_s Starting value of *s*
#'@param max_s Maximum possible *s* value. If the matrix is still not stable, NA is returned
#'
#'@export

find_s <- function(Jacobian, step_size = 0.01, start_s = 1, max_s = 1000){
  
  original_diagonal = diag(Jacobian)
  s = start_s
  
  #Step 1: Check if the original matrix is stable
  lambda_d <- eigen(Jacobian)$values %>% Re() %>% max()
  
  #if it is unstable: decrease s until lambda_d becomes negative 
  if (lambda_d > 0){
    while(TRUE){ 
      
      #set diagonal to original diagonal values multiplied with s
      diag(Jacobian) <- original_diagonal * s
      #calculate lambda_d again
      lambda_d <- eigen(Jacobian)$values %>% Re() %>% max()
      
      #if it is negative , the loop can be stopped
      if (lambda_d < 0){
        return(s)
        break
      }
      #if the system is not stable yet, increase s and repeat 
      s = s+step_size
      
      if (s > max_s){
        return(NA)
        break
      }#end if 
    }#end while
  }#end if 
  
  #if it is stable: increase s until lambda_d becomes negative 
  if (lambda_d < 0){
    while(TRUE){ 
      
      #set diagonal to original diagonal values multiplied with s
      diag(Jacobian) <- original_diagonal * s
      #calculate lambda_d again
      lambda_d <- eigen(Jacobian)$values %>% Re() %>% max()
      
      #if it is positive , the loop can be stopped
      if (lambda_d > 0){
        return(s)
        break
      }
      #if the system is not stable yet, increase s and repeat 
      s = s-step_size
      
      if (s > max_s){
        return(NA)
        break
      }#end if 
    }#end while
  }#end if 
  
}#end function

#'Find the relative amount of self-regulation needed to fullfill Levins I 
#'stability criterium
#'
#'Diagonal values are incrementally increases by a factor s until all total feedback 
#'values (calculated via the coefficients of the characteristic polynomial). 
#'This will not work, if the diagonal contains zeros. 
#'Replace missing values before calculating s*!
#'
#'@param Jacobian A Jacobian matrix
#'@param step_size Amount by which s is increases in each step 
#'@param max_s Maximum possible s value. If the matrix is still not stable, NA is returned
#'
#'@export

find_s_F <- function(Jacobian, step_size = 0.01, max_s = 1000){
  
  original_diagonal = diag(Jacobian)
  s = 0
  while(TRUE){ #the loop increases s until lambda_d becomes negative
    
    #set diagonal to original diagonal values multiplied with s
    diag(Jacobian) <- original_diagonal * s
    #calculate lambda_d again
    lambda_d <- eigen(Jacobian)$values %>% Re() %>% max()
    
    F_k <- pracma::charpoly(Jacobian)[c(-1)]*-1
    
    #if all FK-values are negative, the loop can be stopped
    if (sum(F_k>0) == 0){
      return(s)
      break
    }
    #if the system is not stable yet, increase s and repeat 
    s = s+step_size
    
    if (s > max_s){
      return(NA)
      break
    }#end if 
  }#end while
}#end function