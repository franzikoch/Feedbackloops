
#' Full randomisation of an interaction table 
#' 
#' During the full randomisation procedure, ALL interspecific interaction strengths are randomly 
#' reshuffled within the network. Intraspecific interactions/ Diagonal values are not affected.
#. Interaction strengths are reshuffled within the interaction table. They are also reshuffled across the two 
#' interaction columns so that an a_ij value can become an a_ji value and 
#' vice-versa. The function returns a new interaction table
#' that contains two columns with randomised interaction strengths.
#'  
#' Use this table as inputto assemble_jacobian() and specify the new columns to get a fully randomised Jacobian matrix.
#' Note: values are only randomised across existing links! (the interaction table 
#' does not contain non-existent links). Thus, if the Jacobian is reconstructed from 
#' the randomised interaction table, the network will have the same topology as
#' the empirical one. 
#' 
#' @param df interaction table (created by interaction_strenghts())
#' @param ij_col column of a_ij values to randomise (choose scaled or unscaled)
#' @param ji_col column of a_ji values to randomise (scaled or unscaled)
#' 
#' @return the same interaction table but with two additional columns $a_ij_B_rand 
#' and $a_ji_B_rand that contain the same interaction strengths in a randomised order
#' @export

randomize_all <- function(df, ij_col, ji_col){
  
  ###Implements full randomization of all interaction strenghts
  
  #some defensive programming: 
  #check if specified columns exist, raise an error if not
  if((ij_col %in% colnames(df))== FALSE){stop('ij_col does not exist')}
  if((ji_col %in% colnames(df))==FALSE){stop('ij_col does not exist')}
  
  #print a warning if the two specified columns are the same:
  if(ij_col == ji_col){warning("ij_col and ji_col are identical!")}
  
  
  #add two new columns to df to store randomized interaction pairs
  z = length(df$Species_i)
  df$a_ij_rand <- vector("numeric", z)
  df$a_ji_rand <- vector("numeric", z)
  
  #pick all interspecific interactions from df 
  #(intraspecific interactions are not randomized)
  df_inter <-  df[df$Species_i != df$Species_j,]
  n = length(df_inter$Species_i)
  #add the two chosen columns of interactions strengths together and randomise their order
  vec <- sample(c(df_inter[[ij_col]], df_inter[[ji_col]]))
  
  #fill randomized values back into the new columns in df 
  df_inter$a_ij_rand <- vec[1:n]
  df_inter$a_ji_rand <- vec[(1+n):(n*2)]
  
  #sort df_inter values back into the original df 
  df[df$Species_i != df$Species_j,] <- df_inter
  
  return(df)
}



#' Pairiwise randomisation of an interaction table 
#' 
#' During pairwise randomisation, pairs of interaction strengths are kept intact but their 
#' location is reshuffled across the network. In the interaction table this means the following:
#' a_ij and a_ji values that appear in the same row in the original table, will also be in the same
#' row in the pairwise randomised table(all though it is possible that they switch columns). 
#' The function returns a new interaction table that contains two columns with randomised interaction strengths.
#' 
#' Use this table as input to assemble_jacobian() and specify the two new columns
#' to get a pairwise randomised Jacobian matrix.
#' 
#' Note: values are only randomised across existing links! (the interaction table 
#' does not contain non-existent links). Thus, if the Jacobian is reconstructed from 
#' the randomised interaction table, the network will have the same topology as
#' the empirical one!
#' 
#' @param df interaction table (created by interaction_strenghts())
#' @param ij_col column of a_ij values to randomise (choose scaled or unscaled)
#' @param ji_col column of a_ji values to randomise (scaled or unscaled)
#' 
#' @return the same interaction table with additional columns $a_ij_B_pw and $a_ji_B_pw, 
#' containing randomised interaction strengths
#' 
#' @export
#' 
randomize_pw <- function(df, ij_col, ji_col){
  ###Implements pairwise randomizations of interaction strengths in df 
  
  #some defensive programming: 
  #check if specified columns exist, raise an error if not
  if((ij_col %in% colnames(df))== FALSE){stop('ij_col does not exist')}
  if((ji_col %in% colnames(df))==FALSE){stop('ij_col does not exist')}
  
  #print a warning if the two specified columns are the same:
  if(ij_col == ji_col){warning("ij_col and ji_col are identical!")}
  
  #add two new columns to df to store randomized interaction pairs
  z = length(df$Species_i)
  df$a_ij_pw <- vector("numeric", z)
  df$a_ji_pw <- vector("numeric", z)
  
  #pick all interspecific interactions from df 
  #(intraspecific interactions are not randomized)
  df_inter <-  df[df$Species_i != df$Species_j,]
  
  n = length(df_inter$Species_i)
  #each entry of pw_list is one pair of interaction strengths
  pw_list <- vector("list", length = n)
  
  for (i in 1:n){
    #use sample so the order within the interaction pair can be exchanged
    col1 <- df_inter[[ij_col]]
    col2 <- df_inter[[ji_col]]
    
    pw_list[[i]] <- sample(c(col1[i], col2[i]))
  }
  #randomize the order of list items
  pw_list <- sample(pw_list)
  
  #fill the shuffeled items into the df_inter
  for (i in 1:n){
    df_inter$a_ij_pw[i] <- pw_list[[i]][1]
    df_inter$a_ji_pw[i] <- pw_list[[i]][2]
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
#' @export
randomize_pw2 <- function(df){
  ###Implements pairwise randomizations of interaction strengths in df 
  
  #some defensive programming: 
  #check if specified columns exist, raise an error if not
  if((ij_col %in% colnames(df))== FALSE){stop('ij_col does not exist')}
  if((ji_col %in% colnames(df))==FALSE){stop('ij_col does not exist')}
  
  #print a warning if the two specified columns are the same:
  if(ij_col == ji_col){warning("ij_col and ji_col are identical!")}
  
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


#' Constructs a Jacobian matrix with randomised interaction strengths
#' 
#' This is not needed anymore! assemble_jacobian now works with randomised dataframes as well!!
#' 
#' @param df interaction table, containing columns with randomised interaction strengths
#' @param species_list list of species names (can be taken from abundance table)
#' @param column name specifying which columns of interaction strengths should be used to assemble the matrix
#' 
#' @return a Jacobian matrix 

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


#' Randomise a community matrix to have perfect asymmetry 
#' 
#' During the randomisation procedure, all interaction strengths are reordered to 
#' make the 2-link loops in the randomised system as asymmetric as possible. To do this, 
#' all links are ordered by size. Then, the very strongest link is paired with the
#' weakest one, the second strongest with the second weakest etc. The location of 
#' pairwise interactions in the network is chosen randomly. Also, links are randomised 
#' between the ij and ji column, so that the strong links can appear both above and below
#' the diagonal. 
#'
#' 
#' The function returns a new interaction table that contains two columns with randomised interaction strengths,
#' Use this table as input to assemble_jacobian() and specify the new columns to 
#' build the Jacobian matrix.
#' 
#' Note: values are only randomised across existing links! (the interaction table 
#' does not contain non-existent links). Thus, if the Jacobian is reconstructed from 
#' the randomised interaction table, the network will have the same topology as
#' the empirical one!
#' 
#' @param it interaction table (created by interaction_strengths())
#' @param ij_col column of a_ij values to randomise (choose scaled or unscaled)
#' @param ji_col column of a_ji values to randomise (scaled or unscaled)
#' 
#' @return the same interaction table with additional columns $a_ij_asym and $a_ji_asym, 
#' containing randomised interaction strengths
#' 
#' @export

randomise_asymmetric <- function(it, ij_col, ji_col){
  
  #some defensive programming: 
  #check if specified columns exist, raise an error if not
  if((ij_col %in% colnames(it))== FALSE){stop('ij_col does not exist')}
  if((ji_col %in% colnames(it))==FALSE){stop('ij_col does not exist')}
  
  #print a warning if the two specified columns are the same:
  if(ij_col == ji_col){warning("ij_col and ji_col are identical!")}
  
  #add two new columns to df to store randomized interaction pairs
  z = nrow(it)
  it$a_ij_asym <- vector("numeric", z)
  it$a_ji_asym <- vector("numeric", z)
  
  #pick all interspecific interactions from df 
  #(intraspecific interactions are not randomized)
  it_inter <-  it[it$Species_i != it$Species_j,]
  n = nrow(it_inter)
  #add the two chosen columns of interactions strengths together and sort them by size 
  all_aij <- c(it_inter[[ij_col]], it_inter[[ji_col]])
  all_aij_ordered <- all_aij[order(all_aij)]
  
  #each entry of the pairwise list is one pair of interaction strengths
  #that belong together
  pw_list <- vector("list", length = n)
  
  index1 = 1 #starts at the first list element
  index2 = length(all_aij_ordered) #starts at the last list element
  
  #fill the list of pairwise interactions from the list of ordered interaction strengths
  for (i in 1:n){
    
    col1 <- all_aij_ordered[index1] #stronger link
    col2 <- all_aij_ordered[index2] #weaker link
    
    #use sample, so that the column can be exchanged -> whether the
    #stronger value appears below or above the diagonal is randomised 
    pw_list[[i]] <- sample(c(col1, col2))
    
    #update indices 
    index1 <- index1 + 1
    index2 <- index2 - 1
  }
  
  #randomize the order of list items
  pw_list <- sample(pw_list)
  
  #fill the shuffeled items back into the interaction table
  for (i in 1:n){
    it_inter$a_ij_asym[i] <- pw_list[[i]][1]
    it_inter$a_ji_asym[i] <- pw_list[[i]][2]
  }
  #sort df_inter values back into the original interaction table
  it[it$Species_i != it$Species_j,] <- it_inter
  
  return(it)
}


#' Randomise a community matrix to have perfect asymmetry and hierarchy
#' 
#' During the randomisation procedure, all interaction strengths are reordered to 
#' make the 2-link loops in the randomised system as asymmetric as possible. To do this, 
#' all links are ordered by size. Then, the very strongest link is paired with the
#' weakest one, the second strongest with the second weakest etc. The location of 
#' pairwise interactions in the network is chosen randomly. The stronger link of each 
#' interaction is placed below the diagonal of the matrix, while all weaker links are
#' placed above the diagonal. This makes sure that the matrix is also perfectly hierarchical, 
#' meaning the the weight of all loops > 2 is reduced as well. 
#' 
#' The function returns a new interaction table that contains two columns with randomised interaction strengths,
#' Use this table as input to assemble_jacobian() and specify the new columns to 
#' build the Jacobian matrix.
#' 
#' Note: values are only randomised across existing links! (the interaction table 
#' does not contain non-existent links). Thus, if the Jacobian is reconstructed from 
#' the randomised interaction table, the network will have the same topology as
#' the empirical one!
#' 
#' @param it interaction table (created by interaction_strengths())
#' @param ij_col column of a_ij values to randomise (choose scaled or unscaled)
#' @param ji_col column of a_ji values to randomise (scaled or unscaled)
#' 
#' @return the same interaction table with additional columns $a_ij_asym_h and $a_ji_asym_h, 
#' containing randomised interaction strengths
#' 
#' @export

randomise_asymmetric_hierarchical <- function(it, ij_col, ji_col){
 
  #some defensive programming: 
  #check if specified columns exist, raise an error if not
  if((ij_col %in% colnames(it))== FALSE){stop('ij_col does not exist')}
  if((ji_col %in% colnames(it))==FALSE){stop('ij_col does not exist')}
  
  #print a warning if the two specified columns are the same:
  if(ij_col == ji_col){warning("ij_col and ji_col are identical!")}
  
  #add two new columns to the interaction table to store randomized interaction pairs
  z = nrow(it)
  it$a_ij_asym_h <- vector("numeric", z)
  it$a_ji_asym_h <- vector("numeric", z)
  
  #pick all interspecific interactions from it
  #(intraspecific interactions are not randomized)
  it_inter <-  it[it$Species_i != it$Species_j,]
  n = nrow(it_inter)
  #add the two chosen columns of interactions strengths together and sort them by size 
  all_aij <- c(it_inter[[ij_col]], it_inter[[ji_col]])
  all_aij_ordered <- all_aij[order(all_aij)]
  
  #each entry of the pairwise list is one pair of interaction strengths
  #that belong together
  pw_list <- vector("list", length = n)
  
  index1 = 1 #starts at the first list element
  index2 = length(all_aij_ordered) #starts at the last list element
  
  #fill the list of pairwise interactions from the list of ordered interaction strengths
  for (i in 1:n){
  
    col1 <- all_aij_ordered[index1] #stronger link
    col2 <- all_aij_ordered[index2] #weaker link
    
    pw_list[[i]] <- c(col1, col2)
    
    #update indices 
    index1 <- index1 + 1
    index2 <- index2 - 1
    
  }
  
  #randomize the order of list items/ the location of links in the matrix
  pw_list <- sample(pw_list)
  
  #fill the shuffeled items back into the interaction table
  for (i in 1:n){
    #the stronger link is filled below the diagonal 
    it_inter$a_ji_asym_h[i] <- pw_list[[i]][1]
    #the weaker link is filled above the diagonal 
    it_inter$a_ij_asym_h[i] <- pw_list[[i]][2]
    
  }
  #sort df_inter values back into the original interaction table
  it[it$Species_i != it$Species_j,] <- it_inter
  
  return(it)
}








