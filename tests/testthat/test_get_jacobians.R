#format species contact matrix and abundance list
testthat::test_that("Read in data functions", {
  competition_path <- "Competition_Signy_1.csv"
  abundance_path <- "Abundance_Signy_1.csv"
  
  expect_snapshot(read_contact_matrix(competition_path))
  expect_snapshot(read_abundance(abundance_path, contact_matrix))
})

#calculation of interaction strengths
testthat::test_that("Calculation interaction strength table", {
  
  competition_path <- "Competition_Signy_1.csv"
  abundance_path <- "Abundance_Signy_1.csv"
  
  contact_matrix <- read_contact_matrix(competition_path)
  abundance <- read_abundance(abundance_path, contact_matrix)
  
  expect_snapshot(interaction_strengths(contact_matrix, abundance, c(-0.1, -0.9, -0.2)))
})
  
#assemble the Jacobian matrix
testthat::test_that("Assemble a Jacobian matrix", {
  
  competition_path <- "Competition_Signy_1.csv"
  abundance_path <- "Abundance_Signy_1.csv"
  
  contact_matrix <- read_contact_matrix(competition_path)
  abundance <- read_abundance(abundance_path, contact_matrix)
  
  interaction_table <- interaction_strengths(contact_matrix, abundance, c(-0.1, -0.9, -0.2))
  
  #output should no change
  expect_snapshot(assemble_jacobian(interaction_table, abundance[,1], "a_ij", "a_ji"))
  
  #error should be thrown when aij/ aji col does not exist
  expect_error(assemble_jacobian(interaction_table, abundance[,1], "aij", "a_ji"), "ij_col does not exist")
  
  #error should be thrown when aij == aji 
  expect_warning(assemble_jacobian(interaction_table, abundance[,1], "a_ij", "a_ij"), "ij_col and ji_col are identical!")
})

