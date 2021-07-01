#test function to calculate loop weights

source("test_get_jacobians.R")

list_l2 <- loops(2, Jacobian)[[2]]
list_l3 <- loops(3, Jacobian)[[2]]
