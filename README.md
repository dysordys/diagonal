# diagonal

Supporting material for the manuscript "Self-regulation and the stability of large ecological networks". This repository contains:

* The folder "Distributions" which contains three R-compatible files: "distr-g-55.RData", "distr-g-75.RData", and "distr-g-95.RData". They can be loaded in R using `load("distr-g-55.RData")` etc. They contain the empirical bivariate distributions from which the matrix entries are drawn. The numbers in the file names refer to the value of the parameter *g* (See Methods in the manuscript, and Section 8 of the Supplementary Information): -0.55, -0.75, or -0.95.
* The folder "EmpiricalWebs" containing the adjacency matrices of the empirical networks.
* The folder "GeneticAlgorithm", with source code and makefiles for hill climbing (HC) and genetic algorithms (GA), plus an R wrapper to run them ("launcher.R"). To compile, simply run "make" in the GA and HC folders.
* The folder "Videos" containing animations of the analytical approximation to network stability, applied to all parameterizations of the empirical networks used in the manuscript. The naming convention follows

  `[name of food web]_[number of species]_Distribution-Gamma-[value of g]_[value of f]_type-[type number]`

  The parameters *g* and *f* are defined in the Methods section of the manuscript, and Section 8 of the Supplementary Information. "Type" refers to whether and what kinds of indirect effects are included; type=1: no indirect effects; type=2: indirect negative effects; type=3: indirect positive effects; type=4: indirect positive and negative effects (see Supplementary Information, Section 8 for details). The videos are essentially animated versions of Figures 15, 18, and 32 in the Supplementary Information.
* The file "prob_stability.R". In this file one can
