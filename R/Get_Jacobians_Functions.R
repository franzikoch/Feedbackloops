#contains functions that are used to calculate Jacobian matrices
#by calculating interaction strengths from species-contact matrices and abundances


#' Read in the raw data files
#'
#'
#' Competition and abundance data must be read in together, to make sure that species match.
#' Some species may appear in the abundance table but not in the species contact matrix, 
#' because they do not interact with anyone yet). These species are dropped from the list. 
#' Species in the abundance table are also reordered to match the order of species
#' in the species contact matrix.
#'
#' @param path_competition path to species contact matrix .csv file
#' @param path_abundance path to abundance .csv file
#' 
#' @return A list containing two data frames: the
#' species contact matrix and a cleaned abundance table 
#' 
#'@export
#'
read_data <- function(path_competition, path_abundance){

  #load competition csv file
  competition <- utils::read.csv(path_competition)
  #column 1 can be removed bc we dont need the species names here
  competition[,1] <- NULL

  #load abundance csv file, (file has no header line therefore FALSE)
  abundance <- utils::read.csv(path_abundance, header = FALSE)
  #drop rows that contain no data
  abundance <- stats::na.omit(abundance)

  #check whether the dataset is from Spitzbergen or from another side
  #Spitzbergen data sets have slightly different formats so they must be treated differently
  if (grepl("Spitzbergen", path_competition, fixed = TRUE)== TRUE){

    #fix the column names in the competition matrix
    cols <- colnames(competition)
    i = 1
    for (n in 1:(length(cols)/2)){
      cols[i+1] <- paste(cols[i],".1", sep= "")
      i = i+2
    }
    colnames(competition) <- cols
    rownames(competition) <- cols

    #spaces in the abundance names list must be replaced with . to match the colnames
    abundance[,1] <- gsub(" ", ".", abundance[,1])

  }else{ #if not use the normal read_data function

    #use colnames to set the row names
    rownames(competition) <- colnames(competition)

  }

  #remove species that dont appear in the competition matrix from the abundance list
  abundance <- abundance[abundance[,1] %in% colnames(competition),]

  #reorder the abundance list so that it matches the order in the competition matrix
  original_order <- rownames(competition)[seq(1,length(competition[,1]), 2)]

  #create a new abundance dataframe
  abundance_new <- data.frame(species = original_order, abundance = 0)

  #fill the new data frame with the right abundances from the old one
  for (i in 1:length(original_order)){
    current_species = as.character(abundance_new$species[i])
    abundance_new$abundance[i] <- abundance[abundance$V1 == current_species,2]
  }

  #return competition and abundance data frames
  return(list(competition, abundance_new))
}


#' Calculate interaction strengths
#'
#' Turns the species-contact matrix into a tabular format with one row per 
#' pairwise interaction. Biomass loss rates are calculated with a given list 
#' of cost-values. Interaction strengths are calculated by dividing biomass loss
#' rate \eqn{B_{ij}} by the abundance of species \eqn{B_j}. 
#' 
#' @param competition species-contact matrix
#' @param abundance abundance table
#' @param cost_list list containing win_cost, loss_cost, draw_cost
#'
#' @return A dataframe containing biomass loss rates and interaction strengths
#'
#' @export 
#' 
#' 
interaction_strengths <- function(competition, abundance, cost_list){
  #calculates interaction strengths from competition and abundance data frames

  #PART1 ##################################################################
  #Turn the competition matrix into a table of confrontations with
  #the number of wins/draws per confrontation

  #prepare a data frame that is filled in the for loops
  df <- data.frame(Species_i = character(), #row species
                   Species_j = character(), #column species
                   Interactions = numeric(), #total number of interactions
                   Wins_i = numeric(), #number of wins by the row species
                   Wins_j = numeric(), #number of wins by the column species
                   Draws = numeric(), #number of draws
                   stringsAsFactors = FALSE
  )

  z = 1 #used to index the dataframe

  for(i in seq(from = 1, to = dim(competition)[2], by = 2)){ #loops through the columns in the competition matrix
    #select row i (-> two rows in the competition matrix)
    current_row <- competition[i:(i+1),]
    for (j in seq(from = 1, to = dim(competition)[2], by = 2)){#loop trough columns in row i
      #select the next 2x2 matrix
      sub_mat <- current_row[,j:(j+1)]
      if ((is.na(sub_mat[2,2])== FALSE) && (sub_mat[2,2] != 0))   #check if there were confrontations, if not #interactions = NA
      {
        #if the submatrix contains data, it is added to the data frame

        #create a new row by adding an empty vector of length 6 to df
        df[z,] <- vector("numeric", length = 6)

        #fill the row with data from the submatrix
        df[z,1] <- rownames(sub_mat)[1] #name of rowspecies is taken from the row selected above
        df[z,2] <- colnames(sub_mat)[1] #name of colspecies is taken from the first col of the submatrix
        df[z,3] <- sub_mat[2,2] #number of interactions -> lower right corner
        df[z,4] <- sub_mat[2,1] #number of wins by column species i -> upper right corner
        df[z,5] <- sub_mat[1,2] #number of wins by row species j -> lower left corner
        df[z,6] <- sub_mat[1,1] #number of draws -> upper left corner

        z = z+1
      }
    }
  }

  #remove empty rows from df
  df <- stats::na.omit(df)


  #PART2 #############################################################
  #Use the table created above to calculate interaction coefficients for each row in the table
  #For interspecific interactions, two coefficients are calculated:
  #Impact of j on i(F_ji) and impact of i on j (F_ij)
  #For intraspecific interactions, only one coefficient is calculated: F_ii
  #Draws need to be counted double in those cases
  #All Fs are then divided by the abundance of the corresponding species
  #----------------------------------------------------------------------

  #define assumed costs for confrontations
  win_cost = cost_list[1]
  loss_cost = cost_list[2]
  draw_cost = cost_list[3]

  #initialize vectors to save interactions strengths
  n = dim(df)[1] #number of confrontations
  F_ii <- vector("numeric", length = n)
  a_ii <- vector("numeric", length = n)
  F_ij <- vector("numeric", length = n)
  a_ij <- vector("numeric", length = n)
  F_ji <- vector("numeric", length = n)
  a_ji <- vector("numeric", length = n)

  for(i in 1:n){#go through all confrontations
    #check whether the confrontation is intra- or interspecific
    if (df$Species_i[i] == df$Species_j[i]){ #if competition is intraspecific
      #calculate F_ii -> impact of species i on itself
      #draws are doubled here
      F_ii[i] <- win_cost*df$Wins_i[i]+ loss_cost*df$Wins_j[i]+ 2*draw_cost*df$Draws[i]
      #divide F_ii by abundance of species i, which is taken from the abundance table
      a_ii[i] <- F_ii[i]/abundance[abundance[,1] == df$Species_i[i],2]
    }else{ #if competition is interspecific
      #impact of species j on species i
      F_ij[i] <- win_cost*df$Wins_i[i]+ loss_cost*df$Wins_j[i] + draw_cost*df$Draws[i]
      #divide F_ij by the abundance of species j
      a_ij[i] <- F_ij[i]/abundance[abundance[,1] == df$Species_j[i],2]
      #impact of species i on species j
      F_ji[i] <- win_cost*df$Wins_j[i]+ loss_cost*df$Wins_i[i] + draw_cost*df$Draws[i]
      #divide F_ji by the abundance of species i
      a_ji[i] <- F_ji[i]/abundance[abundance[,1] == df$Species_i[i],2]
    }
  }
  #add vectors to the df
  df <- cbind(df, F_ii, a_ii,F_ij,a_ij, F_ji, a_ji)

  return(df)
}

#' Construct a community matrix
#'
#' Turns interaction strengths from the tabular format (output of *interaction_strengths()*) 
#' into a community matrix. 
#' 
#' To do this, species names from the abundance table is used to name the 
#' rows and columns of the matrix. The corresponding interaction strengths 
#' are then picked from the dataframe, based on the species name.
#'
#' @param interaction_table dataframe created by the interaction strengths function
#' @param species_list list of species names in the community, can be taken from the abundance table
#' @param ij_col name of the column that contains effects of species i on species j
#' @param ji_col name of the column that contains effects of species j on species j 
#'
#' @return A community matrix 
#'
#' @export 
#' 
#' 
assemble_jacobian <- function(interaction_table, species_list, ij_col, ji_col){
  
  #some defensive programming: 
  #check if specified columns exist, raise an error if not
  if((ij_col %in% colnames(interaction_table))== FALSE){stop('ij_col does not exist')}
  if((ji_col %in% colnames(interaction_table))==FALSE){stop('ij_col does not exist')}
  
  #print a warning if the two specified columns are the same:
  if(ij_col == ji_col){warning("ij_col and ji_col are identical!")}
  
  n <- length(species_list) #length of species list gives us the matrix size
  
  Jacobian <- matrix(0, nrow = n, ncol = n) #intialise a matrix in the correct size,
  #filled with zeros
  
  #create a named vector to match the location with species names
  location= 1:n
  names(location) <- species_list
  
  #iterate over the rows of the data frame
  for (i in seq_along(interaction_table[,1])){
    
    row <- interaction_table[i,] #pick the current row
    
    #finde their location in the species list
    index_i <- location[[row$Species_i]]
    index_j <- location[[row$Species_j]]
    
    if (index_i == index_j){#for intraspecific interactions
      
      #digaonal values are not affected by randomisations, so the same column does not change
      Jacobian[index_i, index_j] <- row[['a_ii']] #fill F_ii value into the matrix
      
    }else{ #for interspecific interactions
      
      #which column is used to pick values from depends on randomisation type
      
      Jacobian[index_i, index_j] <- row[[ij_col]] #fill F_ij value into the matrix
      
      Jacobian[index_j, index_i] <- row[[ji_col]] #fill F_ji value into the matrix
      
    }#end if condition
  }#end for-loop
  
  return(Jacobian)
}#end of function