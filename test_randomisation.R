#testing randomisation functions 
#run test_get_jacobians first to get a Jacobian 

#scale interaction strengths
scaled_interaction_table <- scale_interaction_table(interaction_table = interaction_table,
                                                    replacement = 0.1)

#add new randomised columns to the interaction table data frame 
interaction_table <- interaction_table %>% 
  randomize_all(ij_col = "F_ij_B", ji_col = "F_ji_B") %>% 
  randomize_pw(ij_col = "F_ij_B", ji_col = "F_ji_B")

#assemble randomised Jacobians
Jacobian_rand <- assemble_jacobian(interaction_table, abundance[,1], "a_ij_rand", "a_ij_rand")
Jacobian_pw <- assemble_jacobian(interaction_table, abundance[,1], "a_ij_pw", "a_ji_pw")

#randomise scaled interaction_table
scaled_table <- scale_interaction_table(interaction_table = interaction_table,
                                                    replacement = 0.1)
scaled_table <- scaled_table %>% 
  randomize_all(ij_col = "a_ij_scaled", ji_col = "a_ji_scaled") %>% 
  randomize_pw(ij_col = "a_ij_scaled", ji_col = "a_ji_scaled")

#assembled scaled randomised Jacobians
Jacobian_rand <- assemble_jacobian(scaled_table, abundance[,1], "a_ij_rand", "a_ij_rand")
Jacobian_pw <- assemble_jacobian(scaled_table, abundance[,1], "a_ij_pw", "a_ji_pw")
