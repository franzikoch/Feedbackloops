
#' Full randomisation of an interaction table 
#' 
#' During the full randomisation procedure, ALL interspecific interaction strengths are randomly 
#' reshuffled within the network. Intraspecific interactions/ Diagonal values are not affected.
#. Interaction strengths are reshuffled within the interaction table. They are also reshuffled across the two 
#' interaction columns so that an F_ij value can become an F_ji value and 
#' vice-vers. The function returns a new interaction table
#' that contains two columns with randomised interaction strengths. 
#'  Use this table as input
#' to assemble_jacobian_randomised() to get a fully randomised Jacobian matrix.
#' Note: values are only randomised across existing links! (the interaction table 
#' does not contain non-existent links). Thus, if the Jacobian is reconstructed from 
#' the randomised interaction table, the network will have the same topology as
#' the empirical one. 
#' 
#' @param df interaction table (created by interaction_strenghts())
#' 
#' @return the same interaction table but with two additional columns $ F_ij_B_rand 
#' and $F_ji_B_rand that contain the same interaction strengths in a randomised order

randomize_all <- function(df){
  ###Implements full randomization of all interaction strenghts
  
  #add two new columns to df to store randomized interaction pairs
  z = length(df$Species_i)
  df$F_ij_B_rand <- vector("numeric", z)
  df$F_ji_B_rand <- vector("numeric", z)
  
  #pick all interspecific interactions from df 
  #(intraspecific interactions are not randomized)
  df_inter <-  df[df$Species_i != df$Species_j,]
  n = length(df_inter$Species_i)
  #add the two columns of interactions strengths together and randomise their order
  vec <- sample(c(df_inter$F_ij_B, df_inter$F_ji_B))
  
  #fill randomized values back into the new columns in df 
  df_inter$F_ij_B_rand <- vec[1:n]
  df_inter$F_ji_B_rand <- vec[(1+n):(n*2)]
  
  #sort df_inter values back into the original df 
  df[df$Species_i != df$Species_j,] <- df_inter
  
  return(df)
}



#' Pairiwise randomisation of an interaction table 
#' 
#' During pairwise randomisation, pairs of interaction strengths are kept intact but their 
#' location is reshuffled across the network. In the interaction table this means the following:
#' F_ij and F_ji values that appear in the same row in the original table, will also be in the same
#' row in the pairwise randomised table(all though it is possible that they switch columns). 
#' The function returns a new interaction table that contains two columns with randomised interaction strengths,
#' Use this table as input
#' to assemble_jacobian_randomised() to get a fully randomised Jacobian matrix.
#' Note: values are only randomised across existing links! (the interaction table 
#' does not contain non-existent links). Thus, if the Jacobian is reconstructed from 
#' the randomised interaction table, the network will have the same topology as
#' the empirical one!
#' 
#' @param df interactiontable (created by interaction_strenghts())
#' 
#' @return the same interaction table with additional columns $F_ij_B_pw and $F_ji_B_pw, 
#' containing randomised interaction strengths
#' 
randomize_pw <- function(df){
  ###Implements pairwise randomizations of interaction strengths in df 
  
  #add two new columns to df to store randomized interaction pairs
  z = length(df$Species_i)
  df$F_ij_B_pw <- vector("numeric", z)
  df$F_ji_B_pw <- vector("numeric", z)
  
  #pick all interspecific interactions from df 
  #(intraspecific interactions are not randomized)
  df_inter <-  df[df$Species_i != df$Species_j,]
  
  n = length(df_inter$Species_i)
  #each entry of pw_list is one pair of interaction strengths
  pw_list <- vector("list", length = n)
  
  for (i in 1:n){
    #use sample so the order within the interaction pair can be exchanged
    pw_list[[i]] <- sample(c(df_inter$F_ij_B[i], df_inter$F_ji_B[i]))
  }
  #randomize the order of list items
  pw_list <- sample(pw_list)
  
  #fill the shuffeled items into the df_inter
  for (i in 1:n){
    df_inter$F_ij_B_pw[i] <- pw_list[[i]][1]
    df_inter$F_ji_B_pw[i] <- pw_list[[i]][2]
  }
  #sort df_inter values back into the original df 
  df[df$Species_i != df$Species_j,] <- df_inter
  
  return(df)
}



#' Second version of pairwise randomisation procedure
#' 
#' In this version, pairs are randomly reshuffled but F_ji and F_ij values cannot 
#' be switched!! Thus values always stay on one side of the matrix diagonal, which 
#' should avoid the creation of strong intransitive loops ? Not really well tested though
#' 
#' @param df interaction table (created by interaction_strengths)
#' 
#' @return the same interaction table but with two additional columns $F_ij_B_pw2 
#' and $F_ji_B_pw2, that contain randomised interaction strengths
#' 
randomize_pw2 <- function(df){
  ###Implements pairwise randomizations of interaction strengths in df 
  
  #add two new columns to df to store randomized interaction pairs
  z = length(df$Species_i)
  df$F_ij_B_pw2 <- vector("numeric", z)
  df$F_ji_B_pw2 <- vector("numeric", z)
  
  #pick all interspecific interactions from df 
  #(intraspecific interactions are not randomized)
  df_inter <-  df[df$Species_i != df$Species_j,]
  
  n = length(df_inter$Species_i)
  #each entry of pw_list is one pair of interaction strengths
  pw2_list <- vector("list", length = n)
  
  for (i in 1:n){
    #in this version we do not randomly exchange F_ji and F_ij
    #so that all F_ji stay below the diagonal and all F_ij stay above it
    pw2_list[[i]] <- c(df_inter$F_ij_B[i], df_inter$F_ji_B[i])
  }
  #randomize the order of list items
  pw2_list <- sample(pw2_list)
  
  #fill the shuffeled items into the df_inter
  for (i in 1:n){
    df_inter$F_ij_B_pw2[i] <- pw2_list[[i]][1]
    df_inter$F_ji_B_pw2[i] <- pw2_list[[i]][2]
  }
  #sort df_inter values back into the original df 
  df[df$Species_i != df$Species_j,] <- df_inter
  
  return(df)
}


assemble_jacobian_randomized <- function(df, species_list, column){
  
  ##Used to generate a randomized version of the Jacobian matrix 
  ##Use the parameter column to specify the randomization type:
  ##if column = "random", all interaction strengths are exchanged (except for intraspecific ones)
  ##if column = "pw", pairwise interactions are randomized
  ##in both cases the topology stays the same -> 0's in the original Jacobian 
  ##stay 0 in the randomized versions 
  
  ##prepare a dataframe to fill:
  n <- length(species_list)
  Jacobian<- as.data.frame(matrix(nrow = n, ncol = n), stringsAsFactors = FALSE)
  rownames(Jacobian)<- species_list
  colnames(Jacobian)<- species_list
  
  for (i in 1:n){#loops through the rows of the matrix
    row_species <- species_list[i]
    for (j in 1:n){#loops through the columns of the matrix
      column_species <- species_list[j]
      #select row in df with the corresponding species_i and species_j
      r<- df[(df$Species_i == row_species & df$Species_j == column_species),]
      #check if r contains data (if not, the two species did not interact)
      if (nrow(r) > 0){
        if (i == j){ #when on the diagonal, use intraspecific coefficient
          Jacobian[i,j] <- r$F_ii_B
        } else if (column == "random"){
          #F_ij_B_rand values can be used to fill the matrix above the diagonal
          Jacobian[i,j] <- r$F_ij_B_rand
          #F_ji_B_rand values can be used to fill the matrix below the diagonal
          Jacobian[j,i]<- r$F_ji_B_rand
        } else if (column == "pw"){
          #F_ij_B_pw values can be used to fill the matrix above the diagonal
          Jacobian[i,j] <- r$F_ij_B_pw
          #F_ji_B_pw values can be used to fill the matrix below the diagonal
          Jacobian[j,i]<- r$F_ji_B_pw 
        } else if (column == "pw2"){
          #F_ij_B_pw values can be used to fill the matrix above the diagonal
          Jacobian[i,j] <- r$F_ij_B_pw2
          #F_ji_B_pw values can be used to fill the matrix below the diagonal
          Jacobian[j,i]<- r$F_ji_B_pw2 
      }else{print("Invalid column specified")}
    }#end of nrow(r) > 0 condition
    }#end of inner loop
    }#end of outer loop
  #turn NAs into zeros 
  Jacobian[is.na(Jacobian)]<-0
  return(Jacobian)
}

#functioons to calculate loop weights and strengths

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


loops <- function(n, A){
  #calculates the strengths and weights of all loops of length n within the matrix A
  N <- nrow(A)  #number of species in the community 
  #combn returns a list of all unique n-species subsets of the network
  comb <- combn(c(1:N), n, simplify = FALSE)
  
  #prepare lists to store loop weights and strengths 
  #(we know that there are twice as many possible loops as subsets in the list )
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

#define a function to get the names from one of the path lists
get_names <- function(path_list, string){
  n = length(path_list)
  names <- unlist(strsplit(path_list, string))[seq(2,n*2,2)]
  names_sans <- file_path_sans_ext(names)
  return(list(names, names_sans))
}
