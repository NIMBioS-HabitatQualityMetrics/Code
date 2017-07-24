
# This code iterates to find the occupancy probabilities following:
# C.M. Taylor and R.J. Hall 2012 Metapopulation models for seasonally migratory animals
# It uses equations (3.2) and (3.3) to find equilibrium solutions if they exist.

# REQUIRED INPUTS
# dN a matrix size NUMxNUM that contains migration distances
# AN a vector size 1xNUM that contains node areas
# Parameters: NB, NW, zetaEX, zetaEM, cb, cw, mub, muw

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


# Exponential needed in summation calculation
EXP <- exp(-alpha*dN)
AzetaEX <- AN^zetaEX
AzetaEM <- AN^zetaEM

# Choose number of time steps
END <- 200

for (i in 1:END){
  # One step forward breeding nodes
  for (l in 1:NB){
    SUMq <- 0
    for (j in 1:NW){
      SUMq <- SUMq + n[[NB+j]][i]*AzetaEM[NB+j]*EXP[l,NB+j]
    }
    
    n[[l]][i+1] <- n[[l]][i]+cb*SUMq*(1-n[[l]][i])-(mub/AzetaEX[l])*n[[l]][i]
  }
  
  # One step forward non-breeding nodes
  for (l in 1:NW){
    SUMp <- 0
    for (j in 1:NB){
      SUMp <- SUMp + n[[j]][i]*AzetaEM[j]*EXP[NB+l,j]
    }
    n[[NB+l]][i+1] <- n[[NB+l]][i]+cw*SUMp*(1-n[[NB+l]][i])-(muw/AzetaEX[NB+l])*n[[NB+l]][i]
  }
  
}


