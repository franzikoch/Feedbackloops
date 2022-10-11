testthat::test_that("Scaling procedure", {
  
  competition_path <- test_path("Competition_Signy_1.csv")
  abundance_path <- test_path("Abundance_Signy_1.csv")
  
  contact_matrix <- read_contact_matrix(competition_path)
  abundance <- read_abundance(abundance_path, contact_matrix)
  
  interaction_table <- interaction_strengths(contact_matrix, abundance, c(-0.1, -0.9, -0.2))

  #scale within interaction table
  expect_snapshot(scale_interaction_table(interaction_table, 0.1))
  
  #scaled the matrix directly
  Jac <- assemble_jacobian(interaction_table, abundance[,1], "a_ij", "a_ji")
  Jac <- replace_zeros(Jac, 0.1)
  
  expect_snapshot(scale_matrix(Jac))
  
})



testthat::test_that("Calculation of loop weights", {
  
  Jacobian <- as.matrix(read.csv(test_path("Jacobian_Signy_1.csv")))
  
  list_l2 <- loops(2, Jacobian)[[2]]
  list_l3 <- loops(3, Jacobian)[[2]]
  
  expect_equal(max(list_l2), 0.005826275, tolerance = 0.0001)
  expect_equal(max(abs(list_l3)), 0.0002187744, tolerance = 0.0001)
})

testthat::test_that("Calculation of total feedback", {
  
  Jacobian <- as.matrix(read.csv(test_path("Jacobian_Signy_1.csv")))
  
  #get total feedback at levels 2 and 4
  expect_equal(get_Fk(2, Jacobian), 0.001777698, tolerance = 0.0001)
  expect_equal(get_Fk(3, Jacobian), 0.0008329313, tolerance = 0.0001)
  
  #check if total feedback F2 and F3 of the scaled matrix 
  #corresponds to the sum of 2/ 3 link loops of the scaled matrix
  
  Jacobian <- replace_zeros(Jacobian, 0.1)
  Jacobian_scaled <- scale_matrix(Jacobian)
  
  expect_equal(sum(loops(2, Jacobian_scaled)$loop_strength),  get_Fk(2, Jacobian_scaled))
  expect_equal(sum(loops(3, Jacobian_scaled)$loop_strength),  get_Fk(3, Jacobian_scaled))
  
})

testthat::test_that("Calculation of loop weights (parallelised)", {
  
  Jacobian <- as.matrix(read.csv(test_path("Jacobian_Signy_1.csv")))
  
  list_l2 <- loops_parallel(2, Jacobian, 2)[[2]]
  list_l3 <- loops_parallel(3, Jacobian, 2)[[2]]
  
  expect_equal(max(list_l2), 0.005826275, tolerance = 0.0001)
  expect_equal(max(abs(list_l3)), 0.0002187744, tolerance = 0.0001)
})

testthat::test_that("Calculation of total feedback (parallelised)", {
  
  Jacobian <- as.matrix(read.csv(test_path("Jacobian_Signy_1.csv")))
  
  #get total feedback at levels 2 and 4
  expect_equal(get_Fk(2, Jacobian), 0.001777698, tolerance = 0.0001)
  expect_equal(get_Fk(3, Jacobian), 0.0008329313, tolerance = 0.0001)
  
  #check if total feedback F2 and F3 of the scaled matrix 
  #corresponds to the sum of 2/ 3 link loops of the scaled matrix
  
  Jacobian <- replace_zeros(Jacobian, 0.1)
  Jacobian_scaled <- scale_matrix(Jacobian)
  
  expect_equal(sum(loops_parallel(2, Jacobian_scaled, 2)$loop_strength),  get_Fk(2, Jacobian_scaled))
  expect_equal(sum(loops_parallel(3, Jacobian_scaled, 2)$loop_strength),  get_Fk(3, Jacobian_scaled))
  
})