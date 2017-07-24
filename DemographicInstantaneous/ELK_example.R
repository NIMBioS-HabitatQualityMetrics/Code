# # # # # # # # # # # # # # # # # # # # # # #
# #   ELK EXAMPLE FOR CR CALCULATION      # #
# # # # # # # # # # # # # # # # # # # # # # #

#######################################
#######################################

## ENTER PARAMETER VALUES ##

nSites <- 3
seasons <- 2
breedSeason <- 1 # summer/fall = season 1 # winter/spring = season 2

## HABITAT CHARACTERISTICS ##
# Reproductive rate 
R0 <- matrix(0,nSites,seasons)
SR <- 0.60 # Sex Ratio
R0[,1] <- c(SR*0.68,0,SR*0.86)

# habitat survival
sA_node <- matrix(0,nSites,seasons)
sJ_node <- matrix(0,nSites,seasons)

#summer season survival - adults
sA_node[,1] <- c(.9,0,.9)
#juveniles
sJ_node[,1] <- c(0.65,0,0.65)

#winter season survival - adults
sA_node[,2] <-c(0,.9,.9)
#juveniles
sJ_node[,2] <- c(0,0.72,0.72)

## PATHWAY CHARACTERISTICS ##
# migration probabilities
p <- array(0,c(nSites,nSites,seasons))
#transition probabilities from habitat o to i
#13
p[1,3,1] <- 0.13
#12
p[1,2,1] <- 0.87
#33
p[3,3,1] <- 1

# transition probabilties from habitat i to d
#21
p[2,1,2] <- 1
#31
p[3,1,2] <- 0.11
#33
p[3,3,2] <- 0.89


# Migratory survival probability
sA_edge <- array(0,c(nSites,nSites,seasons))
sJ_edge <- array(0,c(nSites,nSites,seasons))
#Adult survival from  habitat o to i
sA_edge[1,2:3,1] <- 1
sA_edge[3,3,1] <- 1
#Juvenile survival from  habitat o to i
sJ_edge[1,2:3,1] <- 1
sJ_edge[3,3,1] <- 1
#Adult survival from habitat i to d
sA_edge[2:3,1,2] <- 1
sA_edge[3,3,2] <- 1
#Juvenile survival from habitat i to d
sJ_edge[2:3,1,2] <- 1
sJ_edge[3,3,2] <- 1


# Population proportions used in lambda calculation
wo <- matrix(0,nSites,seasons)
wo[,1] <- c(0.475, 0, 0.525) #summer/fall population proportions
wo[,2] <- c(0,0.41325,0.58675) #winter/spring population proportions
## END ENTER PARAMETERS ##




## CALCULATE Co ##

## Source code CR_ONEEQ should be saved in same file/location as CR_example
source("CR_ONEEQ.R") 

# Lambda contains one contribution value per node. 
# It is defined as the contribution for the first time the node is occupied
LAMBDA <- matrix(0,1,seasons)
for(i in 1:nSites){
    LAMBDA[1,1] <- LAMBDA[1,1] + wo[i,1]*CR[i,1]
    LAMBDA[1,2] <- LAMBDA[1,2] + wo[i,2]*CR[i,2]
}

#Habitat contributions
cat("Habitat Contributions (CO): \nSummer/Fall  Winter/Spring \n")
print(CR)
cat("\n\n")

##Growth Rate
cat("Growth Rate (Lambda) for each aniversary date:\nSummer/Fall:\n")
print(LAMBDA[1,1])
cat("Winter/Spring:\n")
print(LAMBDA[1,2])
cat("\n\n")

#Transition contributions for o to i 
cat("Pathway Contrubutions (COI):\nSummer/Fall\n")
print(COI[,,1])
cat("\n\n")


#Transition contributions for i to d 
cat("Pathway Contrubutions (COI):\nWinter/Spring\n")
print(COI[,,2])
cat("\n\n")


## Calculate pathway transitions 
##Breeding season anniversary date
cat("Migratory Path Contributions (COID) Summer/Fall:\n")
cat("C121\n")
print(Coid(c(1,2,1),1, breedSeason, sA_node,sJ_node, R0, sA_edge, sJ_edge))  
cat("C131\n")
print(Coid(c(1,3,1),1, breedSeason, sA_node,sJ_node, R0, sA_edge, sJ_edge))  
cat("C133\n")
print(Coid(c(1,3,3),1, breedSeason, sA_node,sJ_node, R0, sA_edge, sJ_edge)) 
cat("C331\n")
print(Coid(c(3,3,1),1,  breedSeason, sA_node,sJ_node, R0, sA_edge, sJ_edge))  
cat("C333\n")
print(Coid(c(3,3,3),1,  breedSeason, sA_node,sJ_node,R0, sA_edge, sJ_edge))  

cat("\n\n")
##Nonbreeding season anniversary date
cat("Migratory Path Contributions (COID) Winter/Spring:\n")
cat("C313\n")
print(Coid(c(3,1,3),2,  breedSeason, sA_node,sJ_node, R0, sA_edge, sJ_edge))  
cat("C312\n")
print(Coid(c(3,1,2),2,  breedSeason, sA_node,sJ_node, R0, sA_edge, sJ_edge))  
cat("C212\n")
print(Coid(c(2,1,2),2,  breedSeason, sA_node,sJ_node, R0, sA_edge, sJ_edge))  
cat("C213\n")
print(Coid(c(2,1,3),2,  breedSeason, sA_node,sJ_node, R0, sA_edge, sJ_edge))  
cat("C333\n")
print(Coid(c(3,3,3),2,  breedSeason, sA_node,sJ_node, R0, sA_edge, sJ_edge))  
 
