
#' Full randomisation of an interaction table 
#' 
#' During the full randomisation procedure, ALL interspecific interaction strengths are randomly 
#' reshuffled within the network. Intraspecific interactions, the diagonal matrix elements 
#' are not affects and network topology is preserved (zeros remain in place). 
#' 
#' The function returns a new interaction table that contains two new columns $a_ij_rand and
#' $a_ji_rand that contains the same values as the original columns but in a different order.
#' 
#' To get a fully randomised Jacobian matrix, use assemble_jacobian() and specify
#' the new columns. 
#' 
#' @param df interaction table (created by interaction_strenghts())
#' @param ij_col column of \eqn{a_{ij}} values to randomise (choose scaled or unscaled)
#' @param ji_col column of \eqn{a_{ji}} values to randomise (scaled or unscaled)
#' 
#' @return interaction table with two additional columns *$a_ij_rand* 
#' and *$a_ji_rand* that contain the interaction strengths in a randomised order
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



#' Pairiwise (Weak) randomisation of an interaction table 
#' 
#' During pairwise randomisation, pairs of interaction strengths are kept intact but their 
#' location is reshuffled across the network. Intraspecific links and network topology are 
#' preserved.
#' 
#' In the interaction table this means the following: \eqn{a_{ij}} and \eqn{a_{ji}}
#' values that appear in the same row in the original table, will also be in the same
#' row in the pairwise randomised table. However, they can switch from  \eqn{a_{ij}} to 
#' \eqn{a_{ji}} and vice versa, which means above-below diagonal orientation can be 
#' reversed. 
#' 
#' The function returns a new interaction table that contains two new columns $a_ij_pw and
#' $a_ji_pw that contains the same values as the original columns but in a different order.
#' 
#' To get a pairwise randomised Jacobian matrix, use assemble_jacobian() and specify
#' the new columns. 
#' 
#' @param df interaction table (created by interaction_strenghts())
#' @param ij_col column of \eqn{a_{ij}} values to randomise (choose scaled or unscaled)
#' @param ji_col column of \eqn{a_{ji}} values to randomise (scaled or unscaled)
#' 
#' @return Interaction table with two additional columns *$a_ij_pw* and *$a_ji_pw*, 
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



#' Minimal randomisation procedure
#' 
#' In this version, pairs are kept together and their above-below diagonal orientation 
#' is also preserved. 
#' 
#' In the interaction table this means the following: \eqn{a_{ij}} and \eqn{a_{ji}}
#' values that appear in the same row in the original table, will also be in the same
#' row in the pairwise randomised table. In contrast to the pairwise (weak) randomisation, 
#' they cannot switch from  \eqn{a_{ij}} to \eqn{a_{ji}} or vice versa. Only the order of rows
#' in the table is randomised. 
#' 
#' The function returns a new interaction table that contains two new columns $a_ij_min and
#' $a_ji_min that contains the same values as the original columns but in a different order.
#' 
#' To get a minimally randomised Jacobian matrix, use assemble_jacobian() and specify
#' the new columns. 
#' 
#' @param df interaction table (created by interaction_strenghts())
#' @param ij_col column of \eqn{a_{ij}} values to randomise (choose scaled or unscaled)
#' @param ji_col column of \eqn{a_{ji}} values to randomise (scaled or unscaled)
#' 
#' @return Interaction table with two additional columns *$a_ij_min* and *$a_ji_min*, 
#' containing randomised interaction strengths
#' 
#' @export
#' 
randomize_minimal <- function(df, ij_col, ji_col){
  
  #some defensive programming: 
  #check if specified columns exist, raise an error if not
  if((ij_col %in% colnames(df))== FALSE){stop('ij_col does not exist')}
  if((ji_col %in% colnames(df))==FALSE){stop('ij_col does not exist')}
  
  #print a warning if the two specified columns are the same:
  if(ij_col == ji_col){warning("ij_col and ji_col are identical!")}
  
  #add two new columns to df to store randomized interaction pairs
  z = length(df$Species_i)
  df$a_ij_min <- vector("numeric", z)
  df$a_ji_min <- vector("numeric", z)
  
  #pick all interspecific interactions from df 
  #(intraspecific interactions are not randomized)
  df_inter <-  df[df$Species_i != df$Species_j,]
  
  n = length(df_inter$Species_i)
  #each entry of pw_list is one pair of interaction strengths
  min_list <- vector("list", length = n)
  
  for (i in 1:n){
    
    col1 <- df_inter[[ij_col]]
    col2 <- df_inter[[ji_col]]
    
    #in contrast to the pairwise randomisation, values are 
    #not randomised between columns
    min_list[[i]] <- c(col1[i], col2[i])
  }
  #randomize the order of list items -> where in the matrix the pairs are located
  min_list <- sample(min_list)
  
  #fill the shuffeled items into the df_inter
  for (i in 1:n){
    df_inter$a_ij_min[i] <- min_list[[i]][1]
    df_inter$a_ji_min[i] <- min_list[[i]][2]
  }
  #sort df_inter values back into the original df 
  df[df$Species_i != df$Species_j,] <- df_inter
  
  return(df)
}





#' Maximise pairwise asymmetry 
#' 
#' During this randomisation procedure, all interaction strengths are reordered to 
#' make the 2-link loops in the randomised system as asymmetric as possible. Diagonal values
#' and network topology are preserved. 
#' 
#' To do this, all links are ordered by size. Then, the very strongest link is paired with the
#' weakest one, the second strongest with the second weakest etc. The location of 
#' pairwise interactions in the network is chosen randomly. Also, links are randomised 
#' between the ij and ji column, so that the strong links can appear both above and below
#' the diagonal. 
#' 
#' The function returns a new interaction table that contains two new columns $a_ij_asym and
#' $a_ji_asym that contains the same values as the original columns but in a different order.
#' 
#' To get a Jacobian matrix with maximised pairwise asymmetry, use assemble_jacobian() and specify
#' the new columns. 
#' 
#' @param it interaction table (created by interaction_strengths())
#' @param ij_col column of a_ij values to randomise (choose scaled or unscaled)
#' @param ji_col column of a_ji values to randomise (scaled or unscaled)
#' 
#' @return Interaction table with two additional columns *$a_ij_asym* and *$a_ji_asym*, 
#' containing randomised interaction strengths
#' 
#' @export

randomize_asymmetric <- function(it, ij_col, ji_col){
  
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


#' Maximise pairwise and community asymmetry 
#' 
#' During this randomisation procedure, all interaction strengths are reordered to 
#' make the 2-link loops in the randomised system as asymmetric as possible. To do this, 
#' all links are ordered by size. Then, the very strongest link is paired with the
#' weakest one, the second strongest with the second weakest etc. The location of 
#' pairwise interactions in the network is chosen randomly (but network topology 
#' is preserved- off-diagonal zeros remain in place).
#' 
#' To also maximise community asymmetry, the stronger link of each interaction 
#' is placed below the diagonal of the matrix (aji column), while all weaker links
#' are placed above the diagonal (aji column).  
#' 
#' The function returns a new interaction table that contains two new columns 
#' *$a_ij_asym_h* and *$a_ji_asym_h* that contains the same values as the original 
#' columns but in a different order.
#' 
#' To get a Jacobian matrix with maximised pairwise asymmetry, use assemble_jacobian() and specify
#' the new columns. 
#' 
#' @param it interaction table (created by interaction_strengths())
#' @param ij_col column of a_ij values to randomise (choose scaled or unscaled)
#' @param ji_col column of a_ji values to randomise (scaled or unscaled)
#' 
#' @return Interaction table with additional columns *$a_ij_asym_h* and *$a_ji_asym_h*, 
#' containing randomised interaction strengths
#' 
#' @export

randomize_asymmetric_hierarchical <- function(it, ij_col, ji_col){
 
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








