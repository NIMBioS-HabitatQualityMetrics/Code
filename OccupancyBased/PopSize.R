# This code calculates Metapopulation Size for migrating Metapopulation Model. Followes the work of:
# C.M. Taylor and R.J. Hall 2012 Metapopulation models for seasonally migratory animals
# O. Ivaskainen and I. Hanski 2001 Spatially Structured Metapopulation Models: Global and Local Assessment of Metapopulation Capacity

# REQUIRED INPUTS
# si a vector that contains weight parameter - here assume si=Ai 
# gstar a vector that contains the equlibrium solutions

si <- 3*AN
SIZEB <- 0
SIZEW <- 0
for (i in 1:NUM){
  if (i <= NB){ SIZEB <- SIZEB + si[i]*gstar[i] }
  if (i > NB){ SIZEW <- SIZEW + si[i]*gstar[i] }
}