#test function to calculate loop weights

source("test_get_jacobians.R")

list_l2 <- loops(2, Jacobian)[[2]]
list_l3 <- loops(3, Jacobian)[[2]]

#get total feedback at levels 2 and 4
F2 <- get_Fk(1, Jacobian)
F3 <- get_Fk(3, Jacobian)

#check the calculations: 

#For a scaled Jacobian, is F2 = the sum of 2-link loops & F3 = the sum of 3-link loops
Jacobian <- replace_zeros(Jacobian, 0.1)
Jacobian_scaled <- scale_matrix(Jacobian)


#is F2 equal to the sum of 2-link loops?
suml2 <- sum(loops(2, Jacobian_scaled)$loop_strength)
F2 <- get_Fk(2, Jacobian_scaled)

if (suml2 != F2){
  print("F2 not equal to suml2! Something is wrong!")
}

#is F3 equal to the sum of 3-link loops?
suml3 <- sum(loops(3, Jacobian_scaled)$loop_strength)
F3 <- get_Fk(3, Jacobian_scaled)

if (suml3 != abs(F3)){
  print("F3 not equal to suml3! Something is wrong!")
}
