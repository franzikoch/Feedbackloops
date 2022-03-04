#testing randomisation functions 
#run test_get_jacobians first to get a Jacobian 
source("test_get_jacobians.R")

#-WITHOUT SCALING-------------------------------------------------------------
#add new randomised columns to the interaction table data frame 
interaction_table <- interaction_table %>% 
  randomize_all(ij_col = "a_ij", ji_col = "a_ji") %>% 
  randomize_pw(ij_col = "a_ij", ji_col = "a_ji")

#assemble randomised Jacobians
Jacobian_rand <- assemble_jacobian(interaction_table, abundance[,1], "a_ij_rand", "a_ij_rand")
Jacobian_pw <- assemble_jacobian(interaction_table, abundance[,1], "a_ij_pw", "a_ji_pw")


#-WITH SCALING---------------------------------------------------------------
#scale interaction strengths
scaled_table <- scale_interaction_table(interaction_table = interaction_table,
                                                    r_factor = 0.1)

scaled_table <- scaled_table %>% 
  randomize_all(ij_col = "a_ij_scaled", ji_col = "a_ji_scaled") %>% 
  randomize_pw(ij_col = "a_ij_scaled", ji_col = "a_ji_scaled") %>% 
  randomise_asymmetric(ij_col = "a_ij_scaled", ji_col = "a_ji_scaled") %>% 
  randomise_asymmetric_hierarchical(ij_col = "a_ij_scaled", ji_col = "a_ji_scaled")

#assembled scaled randomised Jacobians
Jacobian_rand <- assemble_jacobian(scaled_table, abundance[,1], "a_ij_rand", "a_ij_rand")
Jacobian_pw <- assemble_jacobian(scaled_table, abundance[,1], "a_ij_pw", "a_ji_pw")
Jacobian_asymmetric <- assemble_jacobian(scaled_table, abundance[,1], "a_ij_asym", "a_ji_asym")
Jacobian_asymmetric_hierarchical <- assemble_jacobian(scaled_table, abundance[,1], "a_ij_asym_h", "a_ji_asym_h")
