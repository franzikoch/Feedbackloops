#testing get_jacobians_functions
competition_path <- "vignettes/Competition_Signy_1.csv"
abundance_path <- "vignettes/Abundance_Signy_1.csv"


#format species contact matrix and abundance list
#contact_matrix1 <- read_data(competition_path, abundance_path)[[1]]
#abundance1 <- read_data(competition_path, abundance_path)[[2]]

contact_matrix <- read_contact_matrix(competition_path)
abundance <- read_abundance(abundance_path, contact_matrix)

#calculate interaction strengths
interaction_table <- interaction_strengths(contact_matrix, abundance, c(-0.1, -0.9, -0.2))

#assemble the Jacobian matrix
Jacobian <- assemble_jacobian(interaction_table, abundance[,1], "a_ij", "a_ji")
