library("Matrix")

gen_rand_dag <- function(
    n,
    SNP,
    Tr,
    nL,
    deg,
    prob_pleio,
    lo_mp,
    hi_mp,
    lo_pp,
    hi_pp
) {
    pq <- SNP + Tr + nL
    prob1 <- deg / (SNP)
    prob2 <- deg / (Tr)
    if(prob2 > 1) {
        prob2 = 1
    }
    G <- sparseMatrix(1:pq, 1:pq, x=0)

    for (i in 1:SNP) {
        my_rbin <- rbinom(pq-i, 1, prob1)
        which_element <- which(my_rbin != 0)
        if (length(which_element) != 0) {
            G[i, i + which_element] <- 1
        }
    }

    for (i in 1:SNP) {
        iv <- which(G[i,(SNP + nL + 1):pq] == 1)
        if(length(iv) == 1) {
            my_rbin <- rbinom(pq - (SNP + nL + 1), 1, prob_pleio)
            which_element <- which(my_rbin != 0)
            if (length(which_element) != 0) {
                G[i, (SNP + nL + 1):pq ][-iv] <- my_rbin
            }
        }
    }

    for (j in (SNP+1):pq) {
        my_rbin <- rbinom(pq-j, 1, prob2)
        which_element <- which(my_rbin != 0)
        if (length(which_element) != 0) {
            G[j, j + which_element] <- 1
        }
    }

    AA <- G
    x <- data.frame(matrix(vector(), n, pq))

    # marker effects
    for (i in 1:SNP) {
        # marker descendants
        descendants = which(G[i, 1:SNP] == 1)
        nb <- length(descendants)
        ub <- runif(nb)
        b <- runif(nb, lo_pp, hi_pp) * sign(rnorm(nb))
        AA[i, descendants] <- b

        # trait descendants
        descendants = which(G[i, (SNP + 1):pq] == 1)
        nb <- length(descendants)
        ub <- runif(nb)
        b <- runif(nb, lo_mp, hi_mp) * sign(rnorm(nb))
        AA[i, (descendants + SNP)] <- b
    }

    # trait effects
    for (i in (SNP+1):pq) {
        descendants = which(G[i, ] == 1)
        nb <- length(descendants)
        ub <- runif(nb)
        b <- runif(nb, lo_pp, hi_pp) * sign(rnorm(nb))
        AA[i, descendants] <- b
    }

    for (i in 1:pq) {
        if (sum(G[, i] != 0) == 0) {
            x[, i] <- rnorm(n)
        } else {
            anc <- which(G[, i] == 1)
            wa <- x[, anc, drop = FALSE]
            wa <- model.matrix(as.formula(paste("~0+", paste(names(wa), collapse="+"))), data.frame(wa))
            b <- AA[anc, i]
            g <- (wa %*% b)
            x[, i] <- g + rnorm(n, 0, sqrt(1 - var(g)))
        }
    }

    V1 = paste("X", 1:SNP, sep = "")
    V2 = paste("L", 1:nL, sep = "")
    V3 = paste("Y", 1:Tr, sep = "")
    V  = c(V1, V2, V3)
    colnames(G) <- rownames(G) <- colnames(AA) <- rownames(AA) <- colnames(x) <- V
    return(list(G=G, x=x, A=AA))
}

args = commandArgs(trailingOnly=TRUE)
id = as.numeric(args[1])
outdir = args[2]

set.seed(id)

n = 16000
SNP = 1600
Tr = 10
nL = 2
deg = 3
prob_pleio = 0.2
lo_mp = 0.001
hi_mp = 0.05
lo_pp = 0.01
hi_pp = 0.2

dag <- gen_rand_dag(n, SNP, Tr, nL, deg, prob_pleio, lo_mp, hi_mp, lo_pp, hi_pp)
dag_data <- dag$x[,-c((SNP + 1):(SNP + nL))]
dag_data_mat <- data.matrix(dag_data)
writeMM(dag$A, file=file.path(outdir, paste("true_adj_mat_n",  toString(n), "_SNP_", toString(SNP),"_it_",toString(id),".mtx", sep = "")))

corrs <- cor(dag_data_mat)
writeMM(as(corrs, "sparseMatrix"), file=file.path(outdir, paste("corr_n", toString(n), "_SNP_", toString(SNP),"_it_",toString(id) ,".mtx", sep = "" )))

direct_effects = Diagonal(pq) - t(dag$A)
all_effects <- Matrix::tcrossprod(Matrix::solve(direct_effects, sparse=TRUE, tol=.Machine$double.eps))
writeMM(all_effects, file=file.path(outdir, paste("true_causaleffects_n", toString(n), "_SNP_", toString(SNP),"_it_",toString(id),".mtx", sep = "" )))

trait_effects = all_effects[(SNP + nL + 1):pq, (SNP + nL + 1):pq]
trait_effects[lower.tri(trait_effects, diag = TRUE)] <- 0
writeMM(trait_effects, file=file.path(outdir, paste("true_trait_causaleffects_n", toString(n), "_SNP_", toString(SNP),"_it_",toString(id),".mtx", sep = "" )))

