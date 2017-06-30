
arguments <- commandArgs(trailingOnly=TRUE)
path <- arguments[1] # path to data files
base_selfreg <- as.numeric(arguments[2]) # baseline self-regulation strength
h <- as.numeric(arguments[3]) # Hill exponent for the functional response
numdiaggens <- as.numeric(arguments[4]) # no. times to resample diagonal entries

Jac <- function(B) { # return Jacobian of following system:
    ## dB/dt = as.numeric(r*B + B*X*((e*W)%*%(B^h)) - (B^h)*(t(W)%*%(B*X)) - d*B^2)
    ## where X is defined as below
    X <- 1 / as.numeric(B0^h + W %*% (B^h))
    J <- diag(r - 2*d*B + X*as.numeric((e*W) %*% (B^h)) -
              h * (B^(h-1)) * as.numeric(t(W) %*% (B*X))) +
         h * diag(B*X) %*% (e*W) %*% diag(B^(h-1)) -
         h * diag(B*X^2 * as.numeric((e*W) %*% (B^h))) %*% W %*% diag(B^(h-1)) -
         diag(B^h) %*% t(W) %*% diag(X) +
         h * ((B^h) %o% (B^(h-1))) * ((t(W) %*% diag(B*X^2)) %*% W)
    return(as.matrix(J))
}

getr <- function(B) { # choose r to make B the equilibrium biomasses
    X <- 1 / as.numeric(B0^h + W %*% (B^h))
    return((B^(h-1)*as.numeric(t(W)%*%(B*X)) + d*B - X*as.numeric((e*W)%*%(B^h))))
}

webname <- tail(unlist(strsplit(path, "/")), n=1) # rip name from path string
w <- read.table(paste0(path, "data.txt")) # load adjacency matrix
colnames(w) <- NULL # get rid of column names (for converting to numeric matrix)
w <- data.matrix(w) # convert to matrix

S <- nrow(w) # number of species
omega <- numeric(0) # initialize the matrix omega
rsum <- rowSums(w) # omega[i,j] is the proportion of i's maximum consumption rate,
for (i in 1:S) omega <- rbind(omega, rsum) # targeted at consuming j
omega <- t(omega)
omega[omega!=0] <- 1/omega[omega!=0]
W <- w*omega # only the product is needed in the calculations
A <- w # we use A to calculate trophic levels via the Leontief inverse
for (i in 1:S) if (rsum[i] > 0) A[i,] <- w[i,] / rsum[i]
tl <- as.numeric(solve(diag(rep(1, S)) - A) %*% rep(1, S)) # get trophic levels
tl <- round(tl, 6) # round to 6-digit precision (discard numerical errors)
e <- matrix(0.2, S, S) # species j's assimilation efficiency of species i:
for (i in 1:S) if (tl[i] == 2) e[i,] <- 0.1 # 0.1 for herbivores, otherwise 0.2
B0 <- 1 # half-saturation constant for functional responses
bodymass <- read.table(paste0(path, "bodymass.txt")) # species' body masses [kg]
bodymass <- as.numeric(unlist(bodymass)) # convert to vector
B <- bodymass^(-3/4) # equilibrium biomasses, based on Damuth's law

qs <- seq(0, 7, by=1) + 0.5 # increasing self-regulation strengths
nonzeros <- seq(floor(0.79*S), S, by=1) # how many diagonal entries are nonzero
write(paste0("web\th\tk\tP\tq\tpstab"), stdout()) # header of output

for (q in qs) { # main loop
    selfregstr <- base_selfreg * 2^q # determine strength of self-regulation
    for (k in nonzeros) {
        toshuffle <- c(rep(selfregstr, k), rep(0, S-k)) # set of diagonal coeffs
        stable <- 0 # initialize no. of stable matrices out of numdiaggens ones
        for (count in 1:numdiaggens) { # test numdiaggens matrices for stability
            d <- sample(toshuffle) # randomize diagonal coeffs
            r <- getr(B) # choose r to yield B as the equilibrium solution
            J <- Jac(B) # calculate the Jacobian around that equilibrium
            eJ <- eigen(J, only.values=TRUE)$values # obtain its eigenvalues
            lew <- round(max(Re(eJ)), 6) # get leading eigenvalue, rounded
            if (lew < 0) stable <- stable + 1 # the matrix is stable if lew < 0
        }
        write(paste0(webname, "\t", h, "\t", k, "\t", k/S, "\t",
                     q, "\t", stable/numdiaggens), stdout())
    } ## output: name of web, h (type of funcional response), k, P (equal to k/S),
}     ## q, and the fraction of stable matrices out of numdiaggens ones
