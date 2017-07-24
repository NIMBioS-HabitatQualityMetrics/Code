

# This code calculates contribution metrics for migrating Metapopulation Model. Followes the work of:
# C.M. Taylor and R.J. Hall 2012 Metapopulation models for seasonally migratory animals
# O. Ivaskainen and I. Hanski 2001 Spatially Structured Metapopulation Models: Global and Local Assessment of Metapopulation Capacity

# REQUIRED INPUTS
# EXP a matrix size NUMxNUM that contains exp(-alpha*dij)
# AzetaEX and AzetaEN both vectors size 1xNUM that contains node areas raised to the zeta
# gstar a vector that contains the equlibrium solutions

b_numerator <- matrix(0,NUM,NUM)
for (i in 1:NUM){
  for (j in 1:NUM){
    if (j <= NB){ # j is a breeding node
      b_numerator[i,j]<-gstar[j]*cb*AzetaEX[i]*AzetaEX[j]*EXP[i,j]
    }
    if (j > NB){ # j is a non-breeding node
      b_numerator[i,j]<-gstar[j]*cw*AzetaEX[i]*AzetaEX[j]*EXP[i,j]
    }
  }
}

b <- matrix(0,NUM,NUM)
for (i in 1:NUM){
  for (j in 1:NUM){
    SUM <- sum(b_numerator[i,])
    b[i,j] <- b_numerator[i,j]/SUM
  }
}
b_eigen <- eigen(b)