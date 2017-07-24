
# This code constructs the Landscape Matrix as defined in:
# C.M. Taylor and R.J. Hall 2012 Metapopulation models for seasonally migratory animals
# It creates M and then finds the dominant eigenvalue lambdaMM

# REQUIRED INPUTS
# dN a matrix size NUMxNUM that contains migration distances
# AN a vector size 1xNUM that contains node areas
# Parameters: zetaEX, zetaEM




# Average migration distance to calculate alpha
COUNT <- 0
TOT <- 0
for (i in 1:NUM){
  for (j in 1:NUM){
    if (dN[i,j]<DNE){
      COUNT <- COUNT +1
      TOT <- TOT + dN[i,j]
    }
  }
}
AveD <- TOT/COUNT
if(AveD == 0){ AveD <- 1*10^(-10)}
alpha <- 1/AveD # 1/alpha = average migration distance


# Exponential needed in A submatrix calculation
EXP <- exp(-alpha*dN)


# Construct A submatrix
AzetaEX <- AN^zetaEX
AzetaEM <- AN^zetaEM
A <- matrix(0,NUM,NUM)
for (i in 1:NUM){
  for (j in 1:NUM){
    A[i,j] <- AzetaEX[i]*AzetaEM[j]*EXP[i,j]
  }
}



# Build W matrix and find eigenvalues and eigenvectors
AT <- t(A)
ZERO <- matrix(0,NUM,NUM)
W <- rbind( cbind(ZERO,A), cbind(AT,ZERO))
W_eigen <- eigen(W)
lambdaMM <- max(abs(W_eigen[[1]]))
# Check that the dominant eigenvalue is positive
test <- max(W_eigen[[1]])
if (test != lambdaMM){lambdaMM <- -lambdaMM}

