#Functions to read in raw data files 
#this used to be a single function that return a list of [contact matrix, abundance_table]
#now split into two seperate functions
#to get abundance table, both raw data files are needes

#' Reads in the species-contact matrix and returns a nicely formatted data frame
#' that is more easy to work with in R.
#'
#'@param path_competition path to species contact matrix .csv file
#' 
#'@return a data frame containing the species contact matrix
#'
#'@export
read_contact_matrix <- function(path_competition){
  
  #load competition csv file
  competition <- utils::read.csv(path_competition)
  #column 1 can be removed bc we dont need the species names here
  competition[,1] <- NULL
 
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
    #abundance[,1] <- gsub(" ", ".", abundance[,1])
    
  }else{ #if not use the normal read_data function
    
    #use colnames to set the row names
    rownames(competition) <- colnames(competition)
    
  }
  #return competition and abundance data frames
  return(competition)
}


#'Reads in the abundance raw data file
#'
#'Needs the contact matrix as an input, to crosscheck the species names and 
#'their order
#'
#'@param path_abundance path to abundance data .csv file
#'@param contact_matrix species contact matrix data frame
#'
#'@return a dataframe containing species names and abundances (in number of colonies)
#'
#'@export
read_abundance <- function(path_abundance, contact_matrix){
  
  #load abundance csv file, (file has no header line therefore FALSE)
  abundance <- utils::read.csv(path_abundance, header = FALSE)
  #drop rows that contain no data
  abundance <- stats::na.omit(abundance)
  
  #for Spitzbergen: spaces in the abundance names list must be replaced with . 
  #to match the colnames
  if (grepl("Spitzbergen", path_abundance, fixed = TRUE)== TRUE){
    abundance[,1] <- gsub(" ", ".", abundance[,1])
  }
  
  #remove species that dont appear in the competition matrix from the abundance list
  abundance <- abundance[abundance[,1] %in% colnames(contact_matrix),]
  
  #reorder the abundance list so that it matches the order in the competition matrix
  original_order <- rownames(contact_matrix)[seq(1,length(contact_matrix[,1]), 2)]
  
  #create a new abundance dataframe
  abundance_new <- data.frame(species = original_order, abundance = 0)
  
  #fill the new data frame with the right abundances from the old one
  for (i in 1:length(original_order)){
    current_species = as.character(abundance_new$species[i])
    abundance_new$abundance[i] <- abundance[abundance$V1 == current_species,2]
  }
  
  return(abundance_new)
}
