testthat::test_that("Randomisation functions", {
  
  competition_path <- "Competition_Signy_1.csv"
  abundance_path <- "Abundance_Signy_1.csv"
  
  contact_matrix <- read_contact_matrix(competition_path)
  abundance <- read_abundance(abundance_path, contact_matrix)
  
  interaction_table <- interaction_strengths(contact_matrix, abundance, c(-0.1, -0.9, -0.2))
  
  set.seed(123458) # so that randomisations always return the same thing
  
  expect_snapshot(randomize_all(interaction_table, "a_ij", "a_ji"))
  expect_snapshot(randomize_pw(interaction_table, "a_ij", "a_ji"))
  expect_snapshot(randomize_asymmetric(interaction_table, "a_ij", "a_ji"))
  expect_snapshot(randomize_asymmetric_hierarchical(interaction_table, "a_ij", "a_ji"))
  
  #error should be thrown when aij/ aji col does not exist
  expect_error(randomize_all(interaction_table, "aij", "a_ji"), "ij_col does not exist")
  expect_error(randomize_pw(interaction_table, "aij", "a_ji"), "ij_col does not exist")
  expect_error(randomize_asymmetric(interaction_table, "aij", "a_ji"), "ij_col does not exist")
  expect_error(randomize_asymmetric_hierarchical(interaction_table, "aij", "a_ji"), "ij_col does not exist")
  
  #error should be thrown when aij == aji 
  expect_warning(randomize_all(interaction_table, "a_ij", "a_ij"), "ij_col and ji_col are identical!")
  expect_warning(randomize_pw(interaction_table, "a_ij", "a_ij"), "ij_col and ji_col are identical!")
  expect_warning(randomize_asymmetric(interaction_table, "a_ij", "a_ij"), "ij_col and ji_col are identical!")
  expect_warning(randomize_asymmetric_hierarchical(interaction_table, "a_ij", "a_ij"), "ij_col and ji_col are identical!")
})




