# This code calculates metrics for migrating Metapopulation Model. Followes the work of:
# C.M. Taylor and R.J. Hall 2012 Metapopulation models for seasonally migratory animals
# O. Ivaskainen and I. Hanski 2001 Spatially Structured Metapopulation Models: Global and Local Assessment of Metapopulation Capacity

#############
# IMPORTANT #
#############
# Habitats must be ordered with breeding Habitats listed first and non-breeding Habitats listed second.
# Breeding Habitats are from i=1..NB and non breeding Habitats are from i=NB+1...NW
#############


# Depends on:
# LandscapeMatrix.R
# OccupancyProbabilities.R
# Set the working directory
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)


# Habitat PARAMETERS #

NB <- 2 #number of breeding Habitats
NW <- 2 #number of non-breeding Habitats


# Enter Habitat area and distance data
A1 <- 718 # breeding-only Habitat 1
A2 <-  1093 # year-round breeding Habitat 3
A3 <- 315 # nonbreeding-only Habitat 2
A4 <- A2 # year-round nonbreeding Habitat 4
# Construct vector containing area measurements
AN <- matrix(c(A1,A2,A3,A4), nrow=1)

# Enter intrinsic rates and area dependencies
cb <- 0.75
cw <- 0.75
mub <- 0.25
muw <- 0.25
zetaEM <- .1
zetaEX <- .1

# DNE is a large number implying that the migration route does not exist
DNE <- 1*10^10

# Define migration distances
d11 <- DNE
d12 <- DNE
d13 <- 40
d14 <- 60

d21 <- DNE
d22 <- DNE
d23 <- DNE
d24 <- 0

d31 <- d13
d32 <- DNE
d33 <- DNE
d34 <- DNE

d41 <- d14
d42 <- d24
d43 <- DNE
d44 <- DNE

# Migration distances in matrix 
dN <- rbind(c(d11,d12,d13,d14), c(d21,d22,d23,d24), c(d31,d32,d33,d34), c(d41,d42,d43,d44))
NUM <- ncol(dN) # Total number of Habitats defined
NODES <- c("Habitat 1", "Habitat 3", "Habitat 2", "Habitat 4")

# Initial occupancy probabilities
n <- list()
n[[1]] <- 1 # Habitat 1 / Yellowstone Breeding 
n[[2]] <- 1 # Habitat 2 / Cody Resident Breeding
n[[3]] <- 1 # Habitat 3 / Migrant Non-breeding
n[[4]] <- 1 # Habitat 4 / Cody Non-breeding
ninit <- matrix(c(n[[1]],n[[2]],n[[3]],n[[4]]))


### RUN THE BASELINE MODEL ###

source("OccupancyProbabilities.R")
gstar <- matrix(0,1,NUM)
for (k in 1:NUM){
  gstar[k] <- n[[k]][END+1] # Strong assumption here! -> we assume that by t+1 our solutions are in fact at equilibrium
}
cat("BASELINE \n")
cat("P 1 - Habitat 1 / Yellowstone Breeding = ", n[[1]][END+1],"\n")
cat("P 2 - Habitat 3 / Cody Resident Breeding = ", n[[2]][END+1],"\n")
cat("Q 1 - Habitat 2 /  Migrant Non-breeding = ", n[[3]][END+1],"\n")
cat("Q 2 - Habitat 4 / Cody Mixed Population Non-breeding = ", n[[4]][END+1],"\n\n")
rm(n)

# Calculate metrics based on the landscape matrix M, Metapop Size and Habitat removal
# Variables used to store data for perturbations to M
LAMBDA <- matrix(0,1,NUM+1)
WSTORE <- list()
SIZESTORE <- matrix(0,NUM+1,2)

# Calculate metapopulatin capacity
source("LandscapeMatrix.R")
LAMBDA[1]<-lambdaMM
WSTORE[[1]] <- W

cat("LambdaMM =", lambdaMM, "\n")
persist <- sqrt((mub*muw)/(cb*cw))
cat("Persistence Condition = ", persist, "\n\n" )

# Calculate Metapopulation Size
# For this example we will consider s_i=Ai or number of individuals assuming constant population density
source("PopSize.R")
SIZESTORE[1,] <- matrix(c(SIZEB,SIZEW),1,2)


# Store Baseline Variables
ANbase <- AN
dNbase <- dN
EXPbase <- EXP
NUMbase <- NUM
NBbase <- NB
NWbase <- NW
ninitbase <- ninit
gstarbase <- gstar
AzetaEXbase <- AzetaEX
AzetaEMbase <- AzetaEM
NODEbase <- NODES

rm(AN,dN,NUM, AzetaEX, AzetaEM,gstar,NB,NW)

## NODE REMOVAL PERTURBATION EXAMPLE ##

# Now remove each of the nodes and reclaculate lambda MM and Metapopsize
cat("PERTURBATIONS \n")
for (k in 1:NUMbase){
  
  # Remove each Habitat and redo the calculation
  AN <- ANbase[,-k]
  dN <- dNbase[-k,-k]
  ninit <- ninitbase[-k]
  NODES <- NODEbase[-k]
  NUM <- ncol(dN)
  
  source("LandscapeMatrix.R")
  LAMBDA[k+1]<-lambdaMM
  WSTORE[[k+1]] <- W
  
  # Re-run the Occupancy Probabilities
  if (k <= NBbase) { 
    NB <- NBbase-1 
    NW <- NWbase
  }
  if (k > NBbase) { 
    NB <- NBbase 
    NW <- NWbase-1
  }
  n <- list()
  for (i in 1:NUM){
    n[[i]] <- ninit[i]
  }

  source("OccupancyProbabilities.R")
  gstar <- matrix(0,1,NUM)
  for (l in 1:NUM){
    gstar[l] <- n[[l]][END+1] # Strong assumption here! -> we assume that by t+1 our solutions are in fact at equilibrium
  }
  cat("Remove:", NODEbase[k], "\n Population Persistence of remaining habitats:", NODES, "\n")
  cat(gstar)
  cat("\n")
  
  source("PopSize.R")
  SIZESTORE[k+1,] <- matrix(c(SIZEB,SIZEW),1,2)
  rm(AN,dN,NUM, AzetaEX, AzetaEM,gstar)
}

cat("\n")
cat("METRICS \n")
# Calculate relative patch value - Habitat removal
V_rR <- matrix(0,1,NUMbase+1)
# Use the eigenvalues that were just calculated
for (k in 1:NUMbase+1){
  V_rR[k] <- (LAMBDA[1]-LAMBDA[k]) / LAMBDA[k]
}
cat("V_rR = ", V_rR, "\n")

# Calculate contribution to size - Habitat removal
UB <- matrix(0,1,NUMbase+1)
UL <- matrix(0,1,NUMbase+1)
for (k in 1:NUMbase+1){
  UB[k] <- SIZESTORE[1,1]-SIZESTORE[k,1]
  UL[k] <- SIZESTORE[1,2]-SIZESTORE[k,2]
}
cat("UB = ", UB , "\n")
cat("UL = ", UL , "\n\n")

cat("Metapopulation Size - Baseline = Run 1 - Remove Habiats in Order \n")
cat("Breeding S_B     Nonbreeding S_W\n")
print(SIZESTORE)
cat("\n\n")


# Calculate long term contribution
NUM <- NUMbase
EXP <- EXPbase
AzetaEX <- AzetaEXbase
AzetaEM <- AzetaEMbase
gstar <- gstarbase
source("Contribution.R")
cat("Eigenvalues and Eigenvectors after Perturbation \n Lambda_MM and Long term contribution Wi: \n")
print(b_eigen)
