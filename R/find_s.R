#'Calculates the relative amount of self-regulation needed for stability s*
#'
#'Diagonal values are incrementally increases by a factor s until its dominant
#'eigenvalue becomes negative (=the matrix becomes stable). This will not work, 
#'if the diagonal contains zeros! Replace missing values before calculating s*!
#'
#'@param Jacobian A Jacobian matrix
#'@param step_size Amount by which s is increases in each step 
#'@param max_s Maximum possible s value. If the matrix is still not stable, NA is returned

find_s <- function(Jacobian, step_size = 0.01, max_s = 1000){
  
  original_diagonal = diag(Jacobian)
  s = 0
  while(TRUE){ #the loop increases s until lambda_d becomes negative
    
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
}#end function