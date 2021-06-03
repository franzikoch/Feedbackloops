library(dplyr)
library(tidyr)

scale_interaction_table <- function(interaction_table, replacement){
  
  #missing diagonal values are replaced by mean_aij* replacement
  
  #calculate mean interspecific interaction strengths
  mean_aij <- interaction_table %>% 
    filter(Species_i != Species_j) %>% 
    select(F_ij_B, F_ji_B) %>% 
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
        filter(Species_i == species_i & Species_j == species_i) %>% .$F_ii_B
    
      a_jj <- interaction_table %>% 
        filter(Species_i == species_j & Species_j == species_j) %>% .$F_ii_B
  
    
      #check whether a_ii or a_jj are missing, 
      #set them to mean_aij * replacement factor
      if (length(a_ii) == 0){a_ii <- mean_aij * replacement}
      if (length(a_jj) == 0){a_jj <- mean_aij * replacement}
    
      #scale F_ij by dividing it by a_ii
      a_ij_scaled[i] <- (interaction_table[i,]$F_ij_B / a_ii)*-1
    
      #scale F_ji by dividing it by a_jj
      a_ji_scaled[i] <- (interaction_table[i,]$F_ji_B/ a_jj)*-1
    }
  }

  interaction_table$a_ij_scaled <- a_ij_scaled
  interaction_table$a_ji_scaled <- a_ji_scaled
  
  return(interaction_table)
}

scaled_interaction_table <- scale_interaction_table(interaction_table = interaction_table,
                                                    replacement = 0.1)
#assembled the scaled Jacobian
Jacobian_scaled <- assemble_jacobian(interaction_table = scaled_interaction_table, 
                                     species_list = abundance[,1],
                                     ij_col = "a_ij_scaled",
                                     ji_col = "a_ji_scaled")
