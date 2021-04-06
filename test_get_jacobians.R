#testing get_jacobians_functions
competition_path <- "~/Documents/UNI/Bryozoans_git/Rohdaten/Data_cleaned/Competition_Signy_1.csv"
abundance_path <- "~/Documents/UNI/Bryozoans_git/Rohdaten/Data_cleaned/Abundance_Signy_1.csv"


#format species contact matrix and abundance list
contact_matrix <- read_data(competition_path, abundance_path)[[1]]
abundance <- read_data(competition_path, abundance_path)[[2]]

#calculate interaction strengths
interaction_table <- interaction_strengths(contact_matrix, abundance, c(-0.1, -0.9, -0.2))

#assemble the Jacobian matrix
Jacobian <- assemble_jacobian(interaction_table, abundance[,1])
