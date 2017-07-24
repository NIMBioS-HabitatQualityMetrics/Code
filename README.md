This code is provided to supplement the paper:

A guide to calculating habitat-quality metrics to inform conservation of highly mobile species. Bieri et al [Submitted to Natural Resource Modeling - 2017]

This code was writen by Joanna A. Bieri and Christine Sample.

https://zenodo.org/badge/latestdoi/98235744


** Occupancy Based:

This code contains 5 R files, all five are needed in the working directory. To run the code you need to source ElkOccupancyCODE.R, this will calculate the occupancy metrics reported in the paper, results are outputted to the console. See comments within the code for more information. It contains:
	* Contribution.R
	* ElkOccupancyCODE.R
	* LandscapeMatric.R
	* OccupancyProbabilities.R
	* PopSize.R



** Demographic Instantansous:

This code contains 2 R files, both are needed in the working directory. To run the code you need to source ELK_example.R, this will calculate C^o, C^oi, and C^oid values along with outputing lambda, results are outputted to the console. All parameter values are set at the top of the ELK_example.R file. See comments within the code for more information. It contains:
	* CR_ONEEQ.R
	* ELK_example.R



** Demographic Perturbation Sensitivities:

This code contains 1 R file, it must be in your working directory. To run the code you need to source VecPermutationElk.R, this will calculate the sensitivity of the growth rate, lambda, to the summer/fall and winter/spring demographic parameters, results are outputted to the console. It contains:
	* VecPermutationElk.R


** Demographic Perturbation Approximation:

NOTE: The file structure is important. Any changes to the structure will require changing directory information in the code.

This code contains two files:

1. ElkSpecies - This is where parameters and species functions specific to the elk are saved. To run the code you need to source ElkSimulation.R, which is kept in this folder. This will output final population values along with graphs for the time dependent iterations and final season population values. It contains:
	* Baseline1 Folder - contains two .xlsx files with baseline parameter values for adults and juvelniles.
	* ElkSimulation.R - The main code that runs the simulation.
	* SpeciesFunctions.R - This contains the density dependent functions: f_i(t), p_ij(t) and s_ij(t)


2. NetworkCode1.1 - This contains the general network code. Users should not need to interact with this code. See comments within each file for more information about what each script does. It contains:
	* NetworkOutputs.R
	* NetworkSetup.R
	* NetworkSimulation.R
