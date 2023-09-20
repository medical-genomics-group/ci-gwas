library("pcalg")
library("fastmatch")
library("Matrix")

myargs = commandArgs(trailingOnly=TRUE)

# stem of input fileset
stem=myargs[1]
# number of individuals
n=as.numeric(myargs[2])
alpha=as.numeric(myargs[3])
srfci_mode = myargs[4] # one of [mpu, mpd, std]

if (is.na(stem) ||  is.na(n) || is.na(alpha)) {
    stop("Received NA input argument, or argument is missing")
}

print(sprintf("got args: stem: %s, n: %s, alpha: %s", stem, n, alpha))

if (srfci_mode == "mpu") {
    unsh_triple_pheno_only = TRUE
    force_marker_to_trait_before_R1 = FALSE
    force_marker_to_trait_in_the_end = TRUE
    external_ambiguous_triples = FALSE
} else if (srfci_mode == "mpd") {
    unsh_triple_pheno_only = TRUE
    force_marker_to_trait_before_R1 = TRUE
    force_marker_to_trait_in_the_end = FALSE
    external_ambiguous_triples = FALSE
} else if (srfci_mode == "std") {
    unsh_triple_pheno_only = FALSE
    force_marker_to_trait_before_R1 = FALSE
    force_marker_to_trait_in_the_end = FALSE
    external_ambiguous_triples = FALSE
} else if (srfci_mode == "cusk2") {
    unsh_triple_pheno_only = FALSE
    force_marker_to_trait_before_R1 = FALSE
    force_marker_to_trait_in_the_end = TRUE
    external_ambiguous_triples = TRUE
} else {
    print("mode has to be one of [mpu, mpd, std, cusk2]")
    exit()
}

mdim_path <- paste0(stem, ".mdim")
corr_path <- paste0(stem, ".corr")
adj_path <- paste0(stem, ".adj")
sep_path <- paste0(stem, ".sep")
atr_path <- paste0(stem, ".atr")

x <- read.csv(mdim_path, sep="\t", header=FALSE)
num_var <- as.numeric(x[1])
num_phen <- as.numeric(x[2])
max_level <- as.numeric(x[3])
num_atr <- as.numeric(x[4])

form_sepset <- function(seps, nvar) {
    res <- rep(list(rep(list(NULL), times=nvar)), nvar)
    for (i in 1:nvar) {
        for (j in 1:nvar) {
            for (k in 1:max_level) {
                val <- as.numeric(seps[ (((i - 1) * nvar + j) - 1) * max_level + k ])
                    if (val != -1) {
                        if (is.null(res[[i]][[j]])) {
                            res[[i]][[j]] <- c(val + 1)
                        } else {
                            res[[i]][[j]] <- c(res[[i]][[j]], val + 1)
                        }
                    }
                }
            }
        }
    return(res)
}

corr_file <- file(corr_path, "rb")
corrs <- readBin(corr_file, n=num_var * num_var, double(), size=4)
corrmat <- matrix(corrs, nrow=num_var, byrow=TRUE)

adj_file <- file(adj_path, "rb")
adjs <- readBin(adj_file, n=num_var * num_var, integer(), size=4)
adjmat <- matrix(adjs, nrow=num_var, byrow=TRUE)

sep_file <- file(sep_path, "rb")
seps <- readBin(sep_file, n=num_var * num_var * max_level, integer(), size=4)
sepset <- form_sepset(seps, num_var)

atr_file <- file(atr_path, "rb")
# there is a chance that num_var*num_var is too small
atr <- readBin(atr_file, integer(), size=4, n=num_atr*3)
atrmat <- matrix(atr, ncol=3, byrow=TRUE)

conservative = FALSE 
maj.rule = FALSE
num_individuals = n
suffStat <- list(C=corrmat, n=num_individuals)
indepTest = gaussCItest

adjmat[adjmat==2]<-1
adjmat[adjmat==3]<-1

num_marker <- num_var - num_phen

print("Searching for unshielded triples")
if (unsh_triple_pheno_only) {
    u.t <- find.unsh.triple(adjmat[(num_marker + 1):num_var, (num_marker + 1):num_var], check = FALSE)
    unshTripl <- u.t$unshTripl + num_marker
    m <- ncol(unshTripl)
    p <- num_phen
    unshVect <- vapply(seq_len(m), function(k) triple2numb(p,
        unshTripl[1, k], unshTripl[2, k], unshTripl[3, k]), numeric(1))
    u.t <- list(unshTripl = unshTripl, unshVect = unshVect)
} else {
    u.t <- find.unsh.triple(adjmat, check = FALSE)
}

r.v. <- rfci.vStruc(suffStat = suffStat,
                    indepTest = gaussCItest,
                    alpha,
                    sepset,
                    adjmat, 
                    unshTripl = u.t$unshTripl,
                    unshVect = u.t$unshVect,
                    conservative = (conservative || maj.rule),
                    version.unf = c(1, 1),
                    maj.rule = maj.rule, 
                    verbose = FALSE)

A <- r.v.$amat
sepset <- r.v.$sepset

MMfile_pag_vStruc <- paste0(stem, sprintf("_estimated_pag_%s_after_vStruc.mtx", srfci_mode))
writeMM(as(A, "sparseMatrix"), file=MMfile_pag_vStruc)

if (force_marker_to_trait_before_R1) {
    A[1:num_marker, (num_marker + 1):num_var][A[1:num_marker, (num_marker + 1):num_var]!=0]<-2
    A[(num_marker + 1):num_var, 1:num_marker][A[(num_marker + 1):num_var, 1:num_marker]!=0]<-3
}

# optionally mark unshielded triples a-b-c as unfaithful / ambiguous,
# if b is in maximal sepset, but not in minimal sepset.
if (external_ambiguous_triples) {
    p <- num_var
    for (i in 1:nrow(atrmat)) {
        # increment because of one based indexing
        x = atrmat[i, 1] + 1
        y = atrmat[i, 2] + 1
        z = atrmat[i, 3] + 1
        r.v.$unfTripl <- c(r.v.$unfTripl, triple2numb(p, x, y, z))
    }
}

res <- udag2apag(A, suffStat = suffStat, indepTest = gaussCItest, alpha, sepset, rules = rep(TRUE, 10), 
                 unfVect = r.v.$unfTripl, verbose = FALSE)
Amat <- res$graph

if (force_marker_to_trait_in_the_end) {
    Amat[1:num_marker, (num_marker + 1):num_var][Amat[1:num_marker, (num_marker + 1):num_var] !=0 ] <- 2
    Amat[(num_marker + 1):num_var, 1:num_marker][Amat[(num_marker + 1):num_var, 1:num_marker] !=0 ] <- 3
}

pag_mtx=as(Amat, "sparseMatrix")
writeMM(pag_mtx, file=paste0(stem, sprintf("_estimated_pag_%s.mtx", srfci_mode))
) 
