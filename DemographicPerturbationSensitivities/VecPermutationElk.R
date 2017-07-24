# # # # # # # # # # # # # # # # # # # # # # #
# # ELK EXAMPLE FOR SENSITIVITY ANALYSIS  # #
# # # # # # # # # # # # # # # # # # # # # # #
library(magic) # Need for creation of block diagonals

## ENTER PARAMETERS ##
# Enter Demography parameters
B1w <- rbind(c(1,0),c(0,1)) #NO POPULATION in node 1
B2w <- rbind(c(0.72,0),c(0,0.9)) #Winter survival node 2
B3w <- rbind(c(0.72,0),c(0,0.9)) #Winter survival node 3

B1s <- rbind(c(0,0.6*0.68*0.9),c(0.65,0.9)) #Summer reproduction and summer survival node 1
B2s <- rbind(c(1,0),c(0,1)) #NON POPULATION in node 2
B3s <- rbind(c(0,0.6*0.86*0.9),c(0.65,0.9)) #Summer reproduction and summer survival node 3

# Enter Dispersal parameters
M1w <-rbind(c(0,1*1,0.11*1),c(0,0,0),c(0,0,0.89*1)) #Winter transitions juveniles
M2w <-rbind(c(0,1*1,0.11*1),c(0,0,0),c(0,0,0.89*1)) #Winter transitions adults

M1s <-rbind(c(0,0,0),c(0.87*1,0,0),c(0.13*1,0,1*1)) #Summer transitions juveniles
M2s <-rbind(c(0,0,0),c(0.87*1,0,0),c(0.13*1,0,1*1)) #Summer transitions adults
## END ENTER PARAMETERS ##




# Calculate vec permutation matrix
E11 <- rbind(c(1,0,0), c(0,0,0))
E12 <- rbind(c(0,1,0), c(0,0,0))
E13 <- rbind(c(0,0,1), c(0,0,0))
E21 <- rbind(c(0,0,0), c(1,0,0))
E22 <- rbind(c(0,0,0), c(0,1,0))
E23 <- rbind(c(0,0,0), c(0,0,1))

K1 <- kronecker(E11,t(E11))
K2 <- kronecker(E12,t(E12))
K3 <- kronecker(E13,t(E13))
K4 <- kronecker(E21,t(E21))
K5 <- kronecker(E22,t(E22))
K6 <- kronecker(E23,t(E23))

P <- K1+K2+K3+K4+K5+K6

# Create block diagonal matrices
BW <- adiag(B1w,B2w,B3w)
BS <- adiag(B1s,B2s,B3s)
MW <- adiag(M1w,M2w)
MS <- adiag(M1s,M2s)

LW <- t(P) %*% MW %*% P %*% BW
LS <- t(P) %*% MS %*% P %*% BS
L <- LW %*% LS

# Find eigenvalues and vectors
EigenR <- eigen(L)
EigenL <- eigen(t(L))
MAX <- which.max(abs(EigenR$values))

w <- matrix(EigenR$vectors[,MAX],nrow=nrow(EigenR$vectors))
v <- matrix(EigenL$vectors[,MAX],nrow=nrow(EigenL$vectors))

# Sensitivity to full anual cycle
SensitivityA <- v%*%t(w)/as.numeric(t(v)%*%w)

# Sensitivity to summer demography
SummerFK <- t(P) %*% MW %*% P %*% BW %*% t(P) %*% MS %*% P
SummerGK <- diag(1,nrow(EigenR$vectors))
SensitivityBS <- t(SummerFK)%*%SensitivityA%*%t(SummerGK)

# Sensitivity to winter demography
WinterFK <- t(P) %*% MW %*% P 
WinterGK <- t(P) %*% MS %*% P %*% BS
SensitivityBW <- t(WinterFK)%*%SensitivityA%*%t(WinterGK)


cat("Growth rate sensitivities with respect to:\n\n")
cat("Winter s_2^j", SensitivityBW[3,3], "\n")
cat("Winter s_2^a", SensitivityBW[4,4], "\n")
cat("Winter s_3^j", SensitivityBW[5,5], "\n")
cat("Winter s_3^a", SensitivityBW[6,6], "\n")
cat("Summer s_1^j", SensitivityBS[2,1], "\n")
cat("Summer s_1^a", SensitivityBS[2,2] + 0.6*0.68*SensitivityBS[1,2], "\n")
cat("Summer s_3^j", SensitivityBS[6,5], "\n")
cat("Summer s_3^a", SensitivityBS[6,6] + 0.6*0.86*SensitivityBS[5,6], "\n")
