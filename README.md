# diagonal

Supporting material for the manuscript "Self-regulation and the stability of large ecological networks". This repository contains:

* The folder "EmpiricalMatrices" which contains the adjacency matrices and parameterized versions of the empirical networks in an R-compatible format. The file naming convention follows
 
  `[name of food web]_[number of species]_Distribution-Gamma-[value of g]_[value of f]`

The parameters *g* and *f* are defined in the Methods section of the manuscript, and Section 8 of the Supplementary Information. To obtain the data, simply run "load([filename])" from R; this will load into the workspace the matrices "Adj", "M1", "M2", and "M3". "Adj" is the adjacency matrix; "M1" is the parameterized network without indirect effects; "M2" is the matrix of positive indirect effects; and "M3" is the matrix of indirect negative effects.
* The folder "GeneticAlgorithm", with source code and makefiles for hill climbing (HC) and genetic algorithms (GA), plus an R wrapper to run them ("launcher.R"). To compile, simply execute "make" in the GA and HC folders.
* The folder "Videos" containing animations of the analytical approximation to network stability, applied to all parameterizations of the empirical networks used in the manuscript. The naming convention follows

  `[name of food web]_[number of species]_Distribution-Gamma-[value of *g*]_[value of *f*]_type-[type number]`

The parameters *g* and *f* are defined in the Methods section of the manuscript, and Section 8 of the Supplementary Information. "Type" refers to whether and what kinds of indirect effects are included; type=1: no indirect effects; type=2: indirect negative effects; type=3: indirect positive effects; type=4: indirect positive and negative effects (see Supplementary Information, Section 8 for details).
