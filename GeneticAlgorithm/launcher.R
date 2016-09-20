library(igraph)
library(tools)
library(dplyr)

#### COMMAND LINE ARGUMENTS
args <- commandArgs(trailingOnly=TRUE)
filename <- args[1] # the file with the matrices in it
ratio <- args[2] # on the diagonal -ratio * Re(lambda_1)
link_type <- args[3] # specify link-types to be included in matrix

#### Parameters for HC (Hill Climber) and GA (Genetic Algorithm) --- Adjust as needed
REPEATHC <- 50 # Number of times the HC should be repeated before trying the GA
MAXSTEPSSHC <- "500" # How many steps with no improvement before quitting the HC (keep in string)

REPEATGA <- 3 # Number of times the GA should be repeated before giving up
POPGA <- 500 # Pop size for GA
MAXSTEPSSGA <- "500" # How many steps with no improvement before quitting the GA (keep in string)

STARTSEED <- 650

f_base <- basename(file_path_sans_ext(filename))

#### Compute the structural rank. If the matrix M has size S, then we'll need at least
# S - compute_structural_rank(M)
# nonzero elements on the diagonal
# This is used to see whether we've reached the theoretical minimum
compute_structural_rank <- function(M){
    FW <- (M != 0) * 1
    diag(FW) <- 0
    colnames(FW) <- paste0("C", 1:ncol(FW))
    rownames(FW) <- paste0("R", 1:nrow(FW))
    B <- graph_from_incidence_matrix(FW)
    R <- maximum.bipartite.matching(B)$matching_size
    return(R)
}

## get the location of the temporary file to feed to C
tmpmatfile <- paste0(f_base, "_", link_type, "__TMP")
outfilename <- paste0("output/", f_base, "_", ratio, "_find_mink.csv")
## and read it in to get some information before processing
M <- read.table(tmpmatfile)
S <- nrow(M)
## find the lowest number of diagonal elements reached in past attempts
already_run <- read.csv(outfilename) %>% filter(linktype == link_type) %>% select(current_k)
## if no past attempts, use S
maxk <- ifelse(nrow(already_run) > 0, min(already_run), S)
currentk <- maxk
## get the theoretical minimum number of diagonal elements required (without constraining value)
structural_rank <- compute_structural_rank(M)
mink <- S - structural_rank
## initialize a template for each row
status <- data.frame(filename=filename, linktype=link_type, ratio=ratio,
                     current_k=-2, min_k=mink, type="XX", replicate=0)
## search for the minimum number given constraints on values
found_solution <- TRUE
## search until no new solutions can be found
while(found_solution == TRUE){
    found_solution <- FALSE
    currentk <- currentk - 1
    if (currentk <= mink) { # No point computing, it is the theoretical minimum
        write.table(status %>% mutate(current_k=currentk), file=outfilename,
                    col.names=FALSE, row.names=FALSE, append=TRUE, sep=",")
        break
    }
    # Try to HC several times
    for (i in STARTSEED:(STARTSEED+REPEATHC-1)){
        # Launch the hill climber to see whether the food web can be stabilized
        launchstring <- paste("./HC/FindArrangHC", S, tmpmatfile, MAXSTEPSSHC, i, currentk, ratio)
        print(launchstring)
        output <- system(launchstring, ignore.stderr=TRUE, intern=TRUE)
        print(output)
        if (length(output) != 0) {
            if (!is.na(output) && output == "SUCCESS") {
                found_solution <- TRUE
                ## update the template and insert it as a new row in the results file
                write.table(status %>% mutate(current_k=currentk, replicate=i, type="HC"),
                            file=outfilename, col.names=FALSE, row.names=FALSE, append=TRUE, sep=",")
                ## skip remaining HC's and move on to next value of k (one lower)
                break
            }
        }
    }
    if (found_solution == FALSE){
        # Bring out the big guns!
        # Run a GA for REPEATGA times
        for (i in STARTSEED:(STARTSEED+REPEATGA-1)){
            # Launch the GA to see whether the food web can be stabilized
            launchstring <- paste("./GA/FindArrangGA", S, tmpmatfile,
                                  MAXSTEPSSGA, POPGA, i, currentk, ratio)
            print(launchstring)
            output <- system(launchstring, ignore.stderr=TRUE, intern=TRUE)
            print(output)
            if (length(output) != 0) {
                if (!is.na(output) && output == "SUCCESS") {
                    found_solution <- TRUE
                    ## update the template and insert it as a new row in the results file
                    write.table(status %>% mutate(current_k=currentk, replicate=i, type="GA"),
                                file=outfilename, col.names=FALSE, row.names=FALSE, append=TRUE, sep=",")
                    ## skip remaining GA's and move on to next value of k (one lower)
                    break
                }
            }
        }
    }
}
