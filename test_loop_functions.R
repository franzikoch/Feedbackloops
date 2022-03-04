#test function to calculate loop weights

source("test_get_jacobians.R")

list_l2 <- loops(2, Jacobian)[[2]]
list_l3 <- loops(3, Jacobian)[[2]]

#get total feedback at levels 2 and 4
F2 <- get_Fk(2, Jacobian)
F3 <- get_Fk(3, Jacobian)
