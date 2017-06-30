
# This R script reads in a file containing the adjacency matrix of a food web,
# parameterizes the matrix, assigns a given fraction P of all diagonal
# entries some nonzero value, determines whether the resulting community is
# stable, and repeats this assessment a given number of times. The final
# output is the fraction of cases that ended up stable.
# Usage, from command line:
#
#   Rscript prob_stability.R [foodweb] [distr] [P] [q] [f] [type] [NREPS]
#
# where the command-line arguments are described below.

arguments <- commandArgs(trailingOnly=TRUE) # command-line arguments
foodweb <- arguments[1] # file name (with path) of food web's adjacency matrix
distrfile <- arguments[2] # empirical distribution for parameterization;
                          # with path, file name, and extension ".RData"
P <- as.numeric(arguments[3]) # fraction of self-regulating species
q <- as.numeric(arguments[4]) # strength of self-regulation; will be set
                              # to -2^q times the leading eigenvalue
f <- as.numeric(arguments[5]) # strength of indirect effects relative to direct
                              # effects; set to either 2, 5, or 10
type <- as.numeric(arguments[6]) # what kinds of indirect effects to include:
          # type = 1: direct effects only
          # type = 2: direct effects + indirect negative effects
          # type = 3: direct effects + indirect positive effects
          # type = 4: all direct and indirect effects
NREPS <- as.numeric(arguments[7]) # number of replicates


STABVALUE <- -10^-6 # only matrices with a leading eigenvalue smaller than
                    # STABVALUE are considered stable, to prevent misclassifying
                    # matrices as stable due to numerical rounding errors


# This function creates parameterized matrices out of the adjacency matrix.
# Input:
# - foodweb: path & file name of text file containing the adjacency matrix
# - distrfile: an .RData file containing the bivariate distribution of
#              nonzero pairs of entries
# - f: the factor by which indirect effects are weaker than direct effects
# Output:
# - a list with the following elements:
#   - Adj: the original adjacency matrix
#   - M1: the parameterized matrix of direct effects
#   - M2: the parameterized matrix of indirect negative effects
#   - M3: the parameterized matrix of indirect positive effects
build_matrices <- function(foodweb, distrfile, f){
    load(distrfile) # load the distribution
    Adj <- as.matrix(read.table(foodweb)) # load adjacency matrix of food web
    diag(Adj) <- 0 # remove diagonal entries
    Adj <- Adj * ((Adj + t(Adj)) < 2) # remove bidirectional edges
    hasinter <- (rowSums(Adj) + colSums(Adj)) > 0 # remove noninteracting species
    Adj <- Adj[hasinter, hasinter]
    n <- nrow(Adj) # number of nodes
    l <- sum(Adj) # number of links
    M1 <- matrix(0, n, n) # Initialize matrix of direct interactions
    M2 <- matrix(0, n, n) # Initialize matrix of indirect negative interactions
    M3 <- matrix(0, n, n) # Initialize matrix of indirect positive interactions
    # Matrix M1 -- predator-prey interactions
    pairs <- sample.int(nrow(DistributionPairs), sum(Adj), replace=TRUE)
    M1[Adj == 1] <- DistributionPairs[pairs, 1]
    M1[t(Adj) == 1] <- DistributionPairs[pairs, 2]
    # Matrix M2 -- mutual interference of predators sharing prey
    maxneg <- mean(M1[M1 < 0]) / f
    exploitative_competition <- ((t(Adj) %*% Adj) > 0) & !diag(rep.int(1, n))
    M2[exploitative_competition] <- runif(sum(exploitative_competition), maxneg, 0)
    # Matrix M3 -- mutual benefit of prey from sharing predators
    maxpos <- mean(M1[M1 > 0]) / f
    apparent_competition <- ((Adj %*% t(Adj)) > 0) & !diag(rep.int(1, n))
    M3[apparent_competition] <- runif(sum(apparent_competition), 0, maxpos)
    return(list(Adj=Adj, M1=M1, M2=M2, M3=M3))
}

# This function repeatedly assigns the nonzero diagonal entries randomly and
# evaluates the stability of the resulting matrices, and then returns the
# number of cases that ended up stable.
# Input:
# - M: parameterized matrix with zero diagonal
# - num_nonzero: how many diagonal entries should be nonzero
# - val_nonzero: the value to which all nonzero diagonal entries should be set
# - nreps: number of times to repeat the assignment and evaluation of stability
# Output:
# - the number of cases (out of all the nreps ones) which ended up stable
p_stab <- function(M, num_nonzero, val_nonzero, nreps){
    toshuffle <- c(rep.int(val_nonzero, num_nonzero),
                   rep.int(0, nrow(M) - num_nonzero))
    stable <- unlist(lapply(1:nreps, function(xx) {
        max(Re(eigen(M + diag(sample(toshuffle)),
                     only.values=TRUE)$values)) < STABVALUE
    }))
    return(sum(stable))
}

Alist <- build_matrices(foodweb, distrfile, f) # generate parameterized matrices
S <- nrow(Alist$Adj) # number of species
if (type==1) A <- Alist$M1 # choose appropriate parameterization based on type
if (type==2) A <- Alist$M1 + Alist$M2
if (type==3) A <- Alist$M1 + Alist$M3
if (type==4) A <- Alist$M1 + Alist$M2 + Alist$M3

leading_ew <- max(Re(eigen(A, only.values=TRUE)$values)) # leading eigenvalue
                                               # without any diagonal entries
d <- (-2^q)*leading_ew # strength of self-regulation (for nonzero diagonal entries)


# Write result to stdout. Here round(P*S) is the actual number
# (as opposed to the fraction) of self-regulating species.
write(paste0("Fraction of stable communities: ",
             p_stab(A, round(P*S), d, NREPS), "/", NREPS), stdout())

