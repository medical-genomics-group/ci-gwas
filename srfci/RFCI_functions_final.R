library(Matrix)
#######################################
#######################################

# functions from pcalg package for running RFCI v_ structure part and orientations with sparce matrix

########################################
########################################
find.unsh.triple <- function(g, check = TRUE) {
    require(Matrix)
    if (check) {
        stopifnot(all(g == t(g)))
    }
    m <- 0L
    unshTripl <- matrix(integer(), 3, m)
    if (any(g != 0)) {
        p <- nrow(g)
        indS <- which(g == 1, arr.ind = TRUE)
        for (i in seq_len(nrow(indS))) {
            xy <- indS[i, ]
            x <- xy[1]
            y <- xy[2]
            allZ <- setdiff(which(g[y, ] == 1), x)
            for (z in allZ) {
                if (g[x, z] == 0 && g[z, x] == 0) {
                    unshTripl <- cbind(unshTripl, c(xy, z))
                }
            }
        }
        if ((m <- ncol(unshTripl)) > 0) {
            deleteDupl <- logical(m)
            for (i in seq_len(m)) {
                if (unshTripl[1, i] > unshTripl[
                    3,
                    i
                ]) {
                    deleteDupl[i] <- TRUE
                }
            }
            if (any(deleteDupl)) {
                m <- ncol(unshTripl <- unshTripl[, !deleteDupl,
                    drop = FALSE
                ])
            }
        }
    }
    unshVect <- vapply(seq_len(m), function(k) {
        triple2numb(
            p,
            unshTripl[1, k], unshTripl[2, k], unshTripl[3, k]
        )
    }, numeric(1))
    list(unshTripl = unshTripl, unshVect = unshVect)
}

rule1_order_indp <- function(apag, unfVect = NULL) {
    p <- ncol(apag)
    search_apag <- apag
    ind <- which((apag != 0 & t(apag) == 1), arr.ind = TRUE)
    for (i in seq_len(nrow(ind))) {
        b <- ind[i, 1]
        c <- ind[i, 2]
        indA <- which(search_apag[b, ] != 0 & search_apag[, b] == 2 & search_apag[c, ] == 0 & search_apag[, c] == 0)
        indA <- setdiff(indA, c)
        for (a in indA) {
            if (any(unfVect == triple2numb(p, a, b, c), na.rm = TRUE) ||
                any(unfVect == triple2numb(p, c, b, a), na.rm = TRUE)) {
                # skip if the triple is ambiguous, we are not using those for orientations
                next
            }
            if (apag[c, b] == 1 & apag[b, c] != 0) { # not oriented yet
                apag[b, c] <- 2
                apag[c, b] <- 3
            } else if (apag[c, b] == 2) { # has been orientated in opposing direction before
                apag[b, c] <- 2
            }
        }
    }
    apag
}


rule2_order_indp <- function(apag, unfVect = NULL) {
    search_apag <- apag
    ind <- which((apag == 1 & t(apag) != 0), arr.ind = TRUE)
    for (i in seq_len(nrow(ind))) {
        a <- ind[i, 1]
        c <- ind[i, 2]
        indB <- which(
            ((search_apag[a, ] == 2 & search_apag[, a] == 3) & (search_apag[c, ] != 0 & search_apag[, c] == 2)) |
                ((search_apag[a, ] == 2 & search_apag[, a] != 0) & (search_apag[c, ] == 3 & search_apag[, c] == 2))
        )
        if (length(indB) > 0) {
            # it doesn't matter if it has been orientated in the other direction before or not,
            # we only modify the edge mark in one direction, this may result in a <->
            apag[a, c] <- 2
        }
    }
    apag
}

rule3_order_indp <- function(apag, unfVect = NULL) {
    search_apag <- apag
    p <- ncol(apag)
    ind <- which(apag != 0 & t(apag) == 1, arr.ind = TRUE)
    for (i in seq_len(nrow(ind))) {
        b <- ind[i, 1]
        d <- ind[i, 2]
        indAC <- which(search_apag[b, ] != 0 & search_apag[, b] == 2 & search_apag[, d] == 1 & search_apag[d, ] != 0)
        if (length(indAC) >= 2) {
            comb.indAC <- combn(indAC, 2)
            for (j in seq_len(ncol(comb.indAC))) {
                a <- comb.indAC[1, j]
                c <- comb.indAC[2, j]
                if (apag[a, c] == 0 && apag[c, a] == 0 && c != a) {
                    if (any(unfVect == triple2numb(p, a, d, c), na.rm = TRUE) ||
                        any(unfVect == triple2numb(p, c, d, a), na.rm = TRUE)) {
                        apag[d, b] <- 2
                    }
                }
            }
        }
    }
    apag
}

rule4_order_indp <- function(apag, unfVect = NULL) {
    search_apag <- apag
    ind <- which((apag != 0 & t(apag) == 1), arr.ind = TRUE)
    for (i in seq_len(nrow(ind))) {
        b <- ind[i, 1]
        c <- ind[i, 2]
        indA <- which((search_apag[b, ] == 2 & search_apag[, b] != 0) & (search_apag[c, ] == 3 & search_apag[, c] == 2))
        for (a in indA) {
            if (apag[a, b] == 2 && apag[b, c] == 2 && apag[c, b] == 2) {
                break
            }
            md.path <- minDiscrPath(apag, a, b, c)
            N.md <- length(md.path)
            if (N.md > 1) {
                if (b %in% sepset[[md.path[1]]][[md.path[N.md]]] ||
                    b %in% sepset[[md.path[N.md]]][[md.path[1]]]) {
                    apag[b, c] <- 2
                    if (apag[c, b] != 2) { # don't modify <-> caused by conflicting information
                        apag[c, b] <- 3
                    }
                } else {
                    apag[a, b] <- apag[b, c] <- apag[c, b] <- 2
                }
            }
        }
    }
    apag
}

rule5_order_indp <- function(apag, unfVect = NULL) {
    search_apag <- apag
    p <- ncol(apag)
    ind <- which((apag == 1 & t(apag) == 1), arr.ind = TRUE)
    for (i in seq_len(nrow(ind))) {
        a <- ind[i, 1]
        b <- ind[i, 2]
        indC <- which((search_apag[a, ] == 1 & search_apag[, a] == 1) & (search_apag[b, ] == 0 & search_apag[, b] == 0))
        indC <- setdiff(indC, b)
        indD <- which((search_apag[b, ] == 1 & search_apag[, b] == 1) & (search_apag[a, ] == 0 & search_apag[, a] == 0))
        indD <- setdiff(indD, a)
        for (c in indC) {
            for (d in indD) {
                if (search_apag[c, d] == 1 && search_apag[d, c] == 1) {
                    if (faith.check(path2check, unfVect, p)) {
                        apag[a, b] <- apag[b, a] <- 3
                        apag[a, c] <- apag[c, a] <- 3
                        apag[c, d] <- apag[d, c] <- 3
                        apag[d, b] <- apag[b, d] <- 3
                    }
                } else {
                    ucp <- minUncovCircPath(
                        p,
                        pag = search_apag,
                        path = c(a, c, d, b),
                        unfVect = unfVect
                    )
                    if (length(ucp) > 1) {
                        n <- length(ucp)
                        apag[ucp[1], ucp[n]] <- apag[ucp[n], ucp[1]] <- 3
                        for (j in seq_len(length(ucp) - 1)) {
                            apag[ucp[j], ucp[j + 1]] <- apag[ucp[j + 1], ucp[j]] <- 3
                        }
                    }
                }
            }
        }
    }
    apag
}

rule6_order_indp <- function(apag, unfVect = NULL) {
    search_apag <- apag
    ind <- which((apag != 0 & t(apag) == 1), arr.ind = TRUE)
    for (i in seq_len(nrow(ind))) {
        b <- ind[i, 1]
        c <- ind[i, 2]
        indA <- which(search_apag[b, ] == 3 & search_apag[, b] == 3)
        if (length(indA) > 0) {
            apag[c, b] <- 3
        }
    }
    apag
}

rule7_order_indp <- function(apag, unfVect = NULL) {
    search_apag <- apag
    p <- ncol(apag)
    ind <- which((apag != 0 & t(apag) == 1), arr.ind = TRUE)
    for (i in seq_len(nrow(ind))) {
        b <- ind[i, 1]
        c <- ind[i, 2]
        indA <- which((search_apag[b, ] == 3 & search_apag[, b] == 1) & (search_apag[c, ] == 0 & search_apag[, c] == 0))
        indA <- setdiff(indA, c)
        for (a in indA) {
            if (apag[c, b] == 3) {
                break
            }
            if (any(unfVect == triple2numb(p, a, b, c), na.rm = TRUE) ||
                any(unfVect == triple2numb(p, c, b, a), na.rm = TRUE)) {
                next
            }
            apag[c, b] <- 3
        }
    }
    apag
}

rule8_order_indp <- function(apag, unfVect = NULL) {
    search_apag <- apag
    p <- ncol(apag)
    ind <- which((apag == 2 & t(apag) == 1), arr.ind = TRUE)
    for (i in seq_len(nrow(ind))) {
        a <- ind[i, 1]
        c <- ind[i, 2]
        indB <- which(
            ((search_apag[a, ] == 2 & search_apag[, a] == 3) | (search_apag[a, ] == 1 & search_apag[, a] == 3)) &
                (search_apag[c, ] == 3 & search_apag[, c] == 2)
        )
        if (length(indB) > 0) {
            apag[c, a] <- 3
        }
    }
    apag
}

rule9_order_indp <- function(apag, unfVect = NULL) {
    search_apag <- apag
    p <- ncol(apag)
    ind <- which((apag == 2 & t(apag) == 1), arr.ind = TRUE)
    for (i in seq_len(nrow(ind))) {
        a <- ind[i, 1]
        c <- ind[i, 2]
        indB <- which(
            (search_apag[a, ] == 2 | search_apag[a, ] == 1) &
                (search_apag[, a] == 1 | search_apag[, a] == 3) &
                (search_apag[c, ] == 0 & search_apag[, c] == 0)
        )
        indB <- setdiff(indB, c)
        for (b in indB) {
            if (apag[c, a] == 3) { # nothing to do here anymore
                break
            }
            upd <- minUncovPdPath(p, search_apag, a, b, c, unfVect = unfVect)
            if (length(upd) > 1) {
                apag[c, a] <- 3
            }
        }
    }
    apag
}

rule10_order_indp <- function(apag, unfVect = NULL) {
    search_apag <- apag
    p <- ncol(apag)
    ind <- which((apag == 2 & t(apag) == 1), arr.ind = TRUE)
    for (i in seq_len(nrow(ind))) {
        a <- ind[i, 1]
        c <- ind[i, 2]
        indB <- which((search_apag[c, ] == 3 & search_apag[, c] == 2))
        for (b in indB) {
            if (apag[c, a] != 1) {
                break
            }
            for (d in indB) {
                if (b == d || apag[c, a] != 1) {
                    next
                }
                if (
                    (search_apag[a, b] == 1 || search_apag[a, b] == 2) &&
                        (search_apag[b, a] == 1 || search_apag[b, a] == 3) &&
                        (search_apag[a, d] == 1 || search_apag[a, d] == 2) &&
                        (search_apag[d, a] == 1 || search_apag[d, a] == 3) &&
                        (search_apag[d, b] == 0 && search_apag[b, d] == 0)) {
                    if (any(unfVect == triple2numb(p, b, a, d), na.rm = TRUE) ||
                        any(unfVect == triple2numb(p, d, a, b), na.rm = TRUE)) {
                        next
                    }
                    apag[c, a] <- 3
                } else {
                    indX <- which(
                        (seach_apag[a, ] == 1 | seach_apag[a, ] == 2) &
                            (seach_apag[, a] == 1 | seach_apag[, a] == 3),
                        arr.ind = TRUE
                    )
                    indX <- setdiff(indX, c)
                    for (pos.1 in indX) {
                        if (apag[c, a] != 1) {
                            break
                        }
                        for (pos.2 in indX) {
                            if (pos.1 == pos.2 || apag[c, a] != 1) {
                                next
                            }
                            tmp1 <- minUncovPdPath(p, search_apag, a, pos.1, b, unfVect = unfVect)
                            tmp2 <- minUncovPdPath(p, search_apag, a, pos.2, d, unfVect = unfVect)
                            if (length(tmp1) > 1 && length(tmp2) > 1 && pos.1 != pos.2 && apag[pos.1, pos.2] == 0) {
                                if (!any(unfVect == triple2numb(p, pos.1, a, pos.2), na.rm = TRUE) &&
                                    !any(unfVect == triple2numb(p, pos.2, a, pos.1), na.rm = TRUE)) {
                                    apag[c, a] <- 3
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    apag
}

udag2apag_ci_gwas <- function(apag, sepset, rules = rep(TRUE, 10), unfVect = NULL) {
    p <- ncol(apag)
    old_apag <- matrix(0, nrow = p, ncol = p)
    while (any(old_apag != apag)) {
        old_apag <- apag
        if (rules[1]) {
            cat("Applying rule 1 \n")
            apag <- rule1_order_indp(apag, unfVect = unfVect)
        }
        if (rules[2]) {
            cat("Applying rule 2 \n")
            apag <- rule2_order_indp(apag, unfVect = unfVect)
        }
        if (rules[3]) {
            cat("Applying rule 3 \n")
            apag <- rule3_order_indp(apag, unfVect = unfVect)
        }
        if (rules[4]) {
            cat("Applying rule 4 \n")
            apag <- rule4_order_indp(apag, unfVect = unfVect)
        }
        if (rules[5]) {
            cat("Applying rule 5 \n")
            apag <- rule5_order_indp(apag, unfVect = unfVect)
        }
        if (rules[6]) {
            cat("Applying rule 6 \n")
            apag <- rule6_order_indp(apag, unfVect = unfVect)
        }
        if (rules[7]) {
            cat("Applying rule 7 \n")
            apag <- rule7_order_indp(apag, unfVect = unfVect)
        }
        if (rules[8]) {
            cat("Applying rule 8 \n")
            apag <- rule8_order_indp(apag, unfVect = unfVect)
        }
        if (rules[9]) {
            cat("Applying rule 9 \n")
            apag <- rule9_order_indp(apag, unfVect = unfVect)
        }
        if (rules[10]) {
            cat("Applying rule 10 \n")
            apag <- rule10_order_indp(apag, unfVect = unfVect)
        }
    }
    list(graph = apag, sepset = sepset)
}

udag2apag <- function(
    apag, suffStat, indepTest, alpha, sepset, rules = rep(TRUE, 10), unfVect = NULL, verbose = FALSE) {
    require(Matrix)
    if (any(apag != 0)) {
        p <- ncol(apag)
        old_apag1 <- matrix(0, nrow = p, ncol = p)
        while (any(old_apag1 != apag)) {
            old_apag1 <- apag
            if (rules[1]) {
                cat("Applying rule 1 \n")
                ind <- which((apag == 2 & t(apag) != 0), arr.ind = TRUE)
                for (i in seq_len(nrow(ind))) {
                    a <- ind[i, 1]
                    b <- ind[i, 2]
                    indC <- which(apag[b, ] != 0 & apag[, b] ==
                        1 & apag[a, ] == 0 & apag[, a] == 0)
                    indC <- setdiff(indC, a)
                    if (length(indC) > 0) {
                        if (length(unfVect) == 0) {
                            apag[b, indC] <- 2
                            apag[indC, b] <- 3
                            if (verbose) {
                                cat("\nRule 1", "\n")
                                cat(
                                    "Orient:", a, "*->", b, "o-*", indC,
                                    "as:", b, "->", indC, "\n"
                                )
                            }
                        } else {
                            for (j in seq_along(indC)) {
                                c <- indC[j]
                                if (!any(unfVect == triple2numb(
                                    p, a, b,
                                    c
                                ), na.rm = TRUE) && !any(unfVect ==
                                    triple2numb(p, c, b, a), na.rm = TRUE)) {
                                    apag[b, c] <- 2
                                    apag[c, b] <- 3
                                    if (verbose) {
                                        cat("\nRule 1", "\n")
                                        cat(
                                            "Conservatively orient:", a, "*->",
                                            b, "o-*", c, "as:", b, "->", c, "\n"
                                        )
                                    }
                                }
                            }
                        }
                    }
                }
            }
            if (rules[2]) {
                cat("Applying rule 2 \n")
                ind <- which((apag == 1 & t(apag) != 0), arr.ind = TRUE)
                for (i in seq_len(nrow(ind))) {
                    a <- ind[i, 1]
                    c <- ind[i, 2]
                    indB <- which(((apag[a, ] == 2 & apag[, a] ==
                        3) & (apag[c, ] != 0 & apag[, c] == 2)) |
                        ((apag[a, ] == 2 & apag[, a] != 0) & (apag[c, ] == 3 & apag[, c] == 2)))
                    if (length(indB) > 0) {
                        apag[a, c] <- 2
                        if (verbose) {
                            cat("\nRule 2", "\n")
                            cat(
                                "Orient:", a, "->", indB, "*->", c,
                                "or", a, "*->", indB, "->", c, "with",
                                a, "*-o", c, "as:", a, "*->", c, "\n"
                            )
                        }
                    }
                }
            }
            if (rules[3]) {
                cat("Applying rule 3 \n")
                ind <- which(apag != 0 & t(apag) == 1, arr.ind = TRUE)
                for (i in seq_len(nrow(ind))) {
                    b <- ind[i, 1]
                    d <- ind[i, 2]
                    indAC <- which(apag[b, ] != 0 & apag[, b] ==
                        2 & apag[, d] == 1 & apag[d, ] != 0)
                    if (length(indAC) >= 2) {
                        if (length(unfVect) == 0) {
                            counter <- 0
                            while ((counter < (length(indAC) - 1)) &&
                                (apag[d, b] != 2)) {
                                counter <- counter + 1
                                ii <- counter
                                while ((ii < length(indAC)) && (apag[
                                    d,
                                    b
                                ] != 2)) {
                                    ii <- ii + 1
                                    if (apag[indAC[counter], indAC[ii]] ==
                                        0 && apag[indAC[ii], indAC[counter]] ==
                                        0) {
                                        apag[d, b] <- 2
                                        if (verbose) {
                                            cat(
                                                "\nRule 3", "\nOrient:", d,
                                                "*->", b, "\n"
                                            )
                                        }
                                    }
                                }
                            }
                        } else {
                            comb.indAC <- combn(indAC, 2)
                            for (j in seq_len(ncol(comb.indAC))) {
                                a <- comb.indAC[1, j]
                                c <- comb.indAC[2, j]
                                if (apag[a, c] == 0 && apag[c, a] ==
                                    0 && c != a) {
                                    if (!any(unfVect == triple2numb(
                                        p,
                                        a, d, c
                                    ), na.rm = TRUE) && !any(unfVect ==
                                        triple2numb(p, c, d, a), na.rm = TRUE)) {
                                        apag[d, b] <- 2
                                        if (verbose) {
                                            cat(
                                                "\nRule 3", "\nConservatively orient:",
                                                d, "*->", b, "\n"
                                            )
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            if (rules[4]) {
                cat("Applying rule 4 \n")
                ind <- which((apag != 0 & t(apag) == 1), arr.ind = TRUE)
                while (length(ind) > 0) {
                    b <- ind[1, 1]
                    c <- ind[1, 2]
                    ind <- ind[-1, , drop = FALSE]
                    indA <- which((apag[b, ] == 2 & apag[, b] !=
                        0) & (apag[c, ] == 3 & apag[, c] == 2))
                    while (length(indA) > 0 && apag[c, b] == 1) {
                        a <- indA[1]
                        indA <- indA[-1]
                        Done <- FALSE
                        while (!Done && apag[a, b] != 0 && apag[a, c] != 0 && apag[b, c] != 0) {
                            cat("")
                            md.path <- minDiscrPath(apag, a, b, c,
                                verbose = verbose
                            )
                            N.md <- length(md.path)
                            if (N.md == 1) {
                                Done <- TRUE
                            } else {
                                if (b %in% sepset[[md.path[1]]][[md.path[N.md]]] ||
                                    b %in% sepset[[md.path[N.md]]][[md.path[1]]]) {
                                    if (verbose) {
                                        cat("\nRule 4", "\n")
                                        cat(
                                            "There is a discriminating path between",
                                            md.path[1], "and", c, "for",
                                            b, ",and", b, "is in Sepset of",
                                            c, "and", md.path[1], ". Orient:",
                                            b, "->", c, "\n"
                                        )
                                    }
                                    apag[b, c] <- 2
                                    apag[c, b] <- 3
                                } else {
                                    if (verbose) {
                                        cat("\nRule 4", "\n")
                                        cat(
                                            "There is a discriminating path between:",
                                            md.path[1], "and", c, "for",
                                            b, ",and", b, "is not in Sepset of",
                                            c, "and", md.path[1], ". Orient",
                                            a, "<->", b, "<->", c, "\n"
                                        )
                                    }
                                    apag[a, b] <- apag[b, c] <- apag[c, b] <- 2
                                }
                                Done <- TRUE
                                # chkE <- checkEdges(suffStat, indepTest,
                                #                    alpha = alpha, apag = apag, sepset = sepset,
                                #                    path = md.path, unfVect = unfVect,
                                #                    verbose = verbose)
                                # sepset <- chkE$sepset
                                # apag <- chkE$apag
                                # unfVect <- c(unfVect, chkE$unfTripl)
                                # if (!chkE$deleted) {
                                #   if (b %in% sepset[[md.path[1]]][[md.path[N.md]]] ||
                                #       b %in% sepset[[md.path[N.md]]][[md.path[1]]]) {
                                #     if (verbose) {
                                #       cat("\nRule 4", "\n")
                                #       cat("There is a discriminating path between",
                                #           md.path[1], "and", c, "for",
                                #           b, ",and", b, "is in Sepset of",
                                #           c, "and", md.path[1], ". Orient:",
                                #           b, "->", c, "\n")
                                #     }
                                #     apag[b, c] <- 2
                                #     apag[c, b] <- 3
                                #   }
                                #   else {
                                #     if (verbose) {
                                #       cat("\nRule 4", "\n")
                                #       cat("There is a discriminating path between:",
                                #           md.path[1], "and", c, "for",
                                #           b, ",and", b, "is not in Sepset of",
                                #           c, "and", md.path[1], ". Orient",
                                #           a, "<->", b, "<->", c, "\n")
                                #     }
                                #     apag[a, b] <- apag[b, c] <- apag[c,
                                #                                      b] <- 2
                                #   }
                                #   Done <- TRUE
                                # }
                            }
                        }
                    }
                }
            }
            if (rules[5]) {
                cat("Applying rule 5 \n")
                ind <- which((apag == 1 & t(apag) == 1), arr.ind = TRUE)
                while (length(ind) > 0) {
                    a <- ind[1, 1]
                    b <- ind[1, 2]
                    ind <- ind[-1, , drop = FALSE]
                    indC <- which((apag[a, ] == 1 & apag[, a] ==
                        1) & (apag[b, ] == 0 & apag[, b] == 0))
                    indC <- setdiff(indC, b)
                    indD <- which((apag[b, ] == 1 & apag[, b] ==
                        1) & (apag[a, ] == 0 & apag[, a] == 0))
                    indD <- setdiff(indD, a)
                    if (length(indC) > 0 && length(indD) > 0) {
                        counterC <- 0
                        while ((counterC < length(indC)) && apag[
                            a,
                            b
                        ] == 1) {
                            counterC <- counterC + 1
                            c <- indC[counterC]
                            counterD <- 0
                            while ((counterD < length(indD)) && apag[
                                a,
                                b
                            ] == 1) {
                                counterD <- counterD + 1
                                d <- indD[counterD]
                                if (apag[c, d] == 1 && apag[d, c] ==
                                    1) {
                                    if (length(unfVect) == 0) {
                                        apag[a, b] <- apag[b, a] <- 3
                                        apag[a, c] <- apag[c, a] <- 3
                                        apag[c, d] <- apag[d, c] <- 3
                                        apag[d, b] <- apag[b, d] <- 3
                                        if (verbose) {
                                            cat("\nRule 5", "\n")
                                            cat(
                                                "There exists an uncovered circle path between",
                                                a, "and", b, ". Orient", a, "-",
                                                b, "and", a, "-", c, "-", d,
                                                "-", b, "\n"
                                            )
                                        }
                                    } else {
                                        path2check <- c(a, c, d, b)
                                        if (faith.check(
                                            path2check, unfVect,
                                            p
                                        )) {
                                            apag[a, b] <- apag[b, a] <- 3
                                            apag[a, c] <- apag[c, a] <- 3
                                            apag[c, d] <- apag[d, c] <- 3
                                            apag[d, b] <- apag[b, d] <- 3
                                            if (verbose) {
                                                cat("\nRule 5", "\n")
                                                cat(
                                                    "There exists a faithful uncovered circle path between",
                                                    a, "and", b, ". Conservatively orient:",
                                                    a, "-", b, "and", a, "-", c,
                                                    "-", d, "-", b, "\n"
                                                )
                                            }
                                        }
                                    }
                                } else {
                                    ucp <- minUncovCircPath(p,
                                        pag = apag,
                                        path = c(a, c, d, b), unfVect = unfVect,
                                        verbose = verbose
                                    )
                                    if (length(ucp) > 1) {
                                        n <- length(ucp)
                                        apag[ucp[1], ucp[n]] <- apag[
                                            ucp[n],
                                            ucp[1]
                                        ] <- 3
                                        for (j in seq_len(length(ucp) - 1)) {
                                            apag[ucp[j], ucp[j + 1]] <- apag[ucp[j +
                                                1], ucp[j]] <- 3
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            if (rules[6]) {
                cat("Applying rule 6 \n")
                ind <- which((apag != 0 & t(apag) == 1), arr.ind = TRUE)
                for (i in seq_len(nrow(ind))) {
                    b <- ind[i, 1]
                    c <- ind[i, 2]
                    indA <- which(apag[b, ] == 3 & apag[, b] ==
                        3)
                    if (length(indA) > 0) {
                        apag[c, b] <- 3
                        if (verbose) {
                            cat(
                                "\nRule 6", "\nOrient:", b, "o-*",
                                c, "as", b, "-*", c, "\n"
                            )
                        }
                    }
                }
            }
            if (rules[7]) {
                cat("Applying rule 7 \n")
                ind <- which((apag != 0 & t(apag) == 1), arr.ind = TRUE)
                for (i in seq_len(nrow(ind))) {
                    b <- ind[i, 1]
                    c <- ind[i, 2]
                    indA <- which((apag[b, ] == 3 & apag[, b] ==
                        1) & (apag[c, ] == 0 & apag[, c] == 0))
                    indA <- setdiff(indA, c)
                    if (length(indA) > 0) {
                        if (length(unfVect) == 0) {
                            apag[c, b] <- 3
                            if (verbose) {
                                cat("\nRule 7", "\n")
                                cat(
                                    "Orient", indA, "-o", b, "o-*", c,
                                    "as", b, "-*", c, "\n"
                                )
                            }
                        } else {
                            for (j in seq_along(indA)) {
                                a <- indA[j]
                                if (!any(unfVect == triple2numb(
                                    p, a,
                                    b, c
                                ), na.rm = TRUE) && !any(unfVect ==
                                    triple2numb(p, c, b, a), na.rm = TRUE)) {
                                    apag[c, b] <- 3
                                    if (verbose) {
                                        cat("\nRule 7", "\n")
                                        cat(
                                            "Conservatively orient:", a,
                                            "-o", b, "o-*", c, "as", b, "-*",
                                            c, "\n"
                                        )
                                    }
                                }
                            }
                        }
                    }
                }
            }
            if (rules[8]) {
                cat("Applying rule 8 \n")
                ind <- which((apag == 2 & t(apag) == 1), arr.ind = TRUE)
                for (i in seq_len(nrow(ind))) {
                    a <- ind[i, 1]
                    c <- ind[i, 2]
                    indB <- which(((apag[a, ] == 2 & apag[, a] ==
                        3) | (apag[a, ] == 1 & apag[, a] == 3)) &
                        (apag[c, ] == 3 & apag[, c] == 2))
                    if (length(indB) > 0) {
                        apag[c, a] <- 3
                        if (verbose) {
                            cat("\nRule 8", "\n")
                            cat(
                                "Orient:", a, "->", indB, "->", c,
                                "or", a, "-o", indB, "->", c, "with",
                                a, "o->", c, "as", a, "->", c, "\n"
                            )
                        }
                    }
                }
            }
            if (rules[9]) {
                cat("Applying rule 9 \n")
                ind <- which((apag == 2 & t(apag) == 1), arr.ind = TRUE)
                while (length(ind) > 0) {
                    a <- ind[1, 1]
                    c <- ind[1, 2]
                    ind <- ind[-1, , drop = FALSE]
                    indB <- which((apag[a, ] == 2 | apag[a, ] ==
                        1) & (apag[, a] == 1 | apag[, a] == 3) &
                        (apag[c, ] == 0 & apag[, c] == 0))
                    indB <- setdiff(indB, c)
                    while ((length(indB) > 0) && (apag[c, a] ==
                        1)) {
                        b <- indB[1]
                        indB <- indB[-1]
                        upd <- minUncovPdPath(p, apag, a, b, c,
                            unfVect = unfVect,
                            verbose = verbose
                        )
                        if (length(upd) > 1) {
                            apag[c, a] <- 3
                            if (verbose) {
                                cat("\nRule 9", "\n")
                                cat(
                                    "There exists an uncovered potentially directed path between",
                                    a, "and", c, ". Orient:", a, " ->",
                                    c, "\n"
                                )
                            }
                        }
                    }
                }
            }
            if (rules[10]) {
                cat("Applying rule 10 \n")
                ind <- which((apag == 2 & t(apag) == 1), arr.ind = TRUE)
                while (length(ind) > 0) {
                    a <- ind[1, 1]
                    c <- ind[1, 2]
                    ind <- ind[-1, , drop = FALSE]
                    indB <- which((apag[c, ] == 3 & apag[, c] ==
                        2))
                    if (length(indB) >= 2) {
                        counterB <- 0
                        while (counterB < length(indB) && (apag[
                            c,
                            a
                        ] == 1)) {
                            counterB <- counterB + 1
                            b <- indB[counterB]
                            indD <- setdiff(indB, b)
                            counterD <- 0
                            while ((counterD < length(indD)) && (apag[
                                c,
                                a
                            ] == 1)) {
                                counterD <- counterD + 1
                                d <- indD[counterD]
                                if ((apag[a, b] == 1 || apag[a, b] ==
                                    2) && (apag[b, a] == 1 || apag[b, a] ==
                                    3) && (apag[a, d] == 1 || apag[a, d] ==
                                    2) && (apag[d, a] == 1 || apag[d, a] ==
                                    3) && (apag[d, b] == 0 && apag[b, d] ==
                                    0)) {
                                    if (length(unfVect) == 0) {
                                        apag[c, a] <- 3
                                        if (verbose) {
                                            cat(
                                                "\nRule 10", "\nOrient:", a,
                                                "->", c, "\n"
                                            )
                                        }
                                    } else {
                                        if (!any(unfVect == triple2numb(
                                            p,
                                            b, a, d
                                        ), na.rm = TRUE) && !any(unfVect ==
                                            triple2numb(p, d, a, b), na.rm = TRUE)) {
                                            apag[c, a] <- 3
                                            if (verbose) {
                                                cat(
                                                    "\nRule 10", "\nConservatively orient:",
                                                    a, "->", c, "\n"
                                                )
                                            }
                                        }
                                    }
                                } else {
                                    indX <- which((apag[a, ] == 1 | apag[a, ] == 2) & (apag[, a] == 1 | apag[
                                        ,
                                        a
                                    ] == 3), arr.ind = TRUE)
                                    indX <- setdiff(indX, c)
                                    if (length(indX >= 2)) {
                                        i1 <- 0
                                        while (i1 < length(indX) && apag[
                                            c,
                                            a
                                        ] == 1) {
                                            i1 <- i1 + 1
                                            pos.1 <- indX[i1]
                                            indX2 <- setdiff(indX, pos.1)
                                            i2 <- 0
                                            while (i2 < length(indX2) && apag[
                                                c,
                                                a
                                            ] == 1) {
                                                i2 <- i2 + 1
                                                pos.2 <- indX2[i2]
                                                tmp1 <- minUncovPdPath(p, apag,
                                                    a, pos.1, b,
                                                    unfVect = unfVect,
                                                    verbose = verbose
                                                )
                                                tmp2 <- minUncovPdPath(p, apag,
                                                    a, pos.2, d,
                                                    unfVect = unfVect,
                                                    verbose = verbose
                                                )
                                                if (length(tmp1) > 1 && length(tmp2) >
                                                    1 && pos.1 != pos.2 && apag[
                                                    pos.1,
                                                    pos.2
                                                ] == 0) {
                                                    if (length(unfVect) == 0) {
                                                        apag[c, a] <- 3
                                                        if (verbose) {
                                                            cat(
                                                                "\nRule 10", "\nOrient:",
                                                                a, "->", c, "\n"
                                                            )
                                                        }
                                                    } else {
                                                        if (!any(unfVect == triple2numb(
                                                            p,
                                                            pos.1, a, pos.2
                                                        ), na.rm = TRUE) &&
                                                            !any(unfVect == triple2numb(
                                                                p,
                                                                pos.2, a, pos.1
                                                            ), na.rm = TRUE)) {
                                                            apag[c, a] <- 3
                                                            if (verbose) {
                                                                cat(
                                                                    "\nRule 10", "\nConservatively orient:",
                                                                    a, "->", c, "\n"
                                                                )
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    list(graph = apag, sepset = sepset)
}

#######################################

minUncovPdPath <- function(p, pag, a, b, c, unfVect, verbose = FALSE) {
    ## Purpose: find a minimal uncovered pd path for a,b,c saved in path.
    ## Check also for the conservative case that it is unambiguous
    ## If a path exists this is the output, otherwise NA
    ## ----------------------------------------------------------------------
    ## Arguments: - p: number of nodes in the graph
    ##            - pag: adjacency matrix
    ##            - a,b,c : nodes under interest
    ##            - unfVect: vector containing the ambiguous triples
    ## ----------------------------------------------------------------------
    ## Author: Diego Colombo, Date: 19 Oct 2011; small changes: Martin Maechler, Joris Mooij

    ## first check whether a,b,c is already a upd path
    stopifnot((pag[a, b] == 1 | pag[a, b] == 2) &
        (pag[b, a] == 1 | pag[b, a] == 3))
    min.upd.path <- NA
    done <- FALSE
    if ((pag[b, c] == 1 | pag[b, c] == 2) &
        (pag[c, b] == 1 | pag[c, b] == 3) &
        (pag[c, a] == 0)) {
        mpath <- c(a, b, c)
        if (length(unfVect) == 0 || ## <<- normal version: just save
            ## conservative version, check the path to be faithful:
            faith.check(mpath, unfVect, p)) {
            ## save the path to be returned
            min.upd.path <- mpath
            if (verbose) {
                cat("    minUncovPdPath: path found: ", mpath, ", uncovered: ", TRUE, "\n")
            }
            done <- TRUE
        }
    }

    ## now check paths of 4 or more nodes of the form <a,b,...,c>
    if (!done) {
        visited <- rep(FALSE, p)
        visited[c(a, b, c)] <- TRUE
        min.upd.path <- NA
        ## find all neighbours of b not visited yet
        indD <- which((pag[b, ] == 1 | pag[b, ] == 2) &
            (pag[, b] == 1 | pag[, b] == 3) &
            (pag[, a] == 0) & !visited)
        if (length(indD) > 0) {
            path.list <- updateList(b, indD, NULL)
            done <- FALSE
            while ((length(path.list) > 0) && (!done)) {
                ## next element in the queue
                mpath <- path.list[[1]]
                m <- length(mpath)
                d <- mpath[m]
                path.list[[1]] <- NULL
                visited[d] <- TRUE
                if (any(pag[d, c] == 1:2) && any(pag[c, d] == c(1, 3))) {
                    ## pd path found
                    mpath <- c(a, mpath, c)
                    n <- length(mpath)
                    ## check the path to be uncovered
                    uncov <- TRUE
                    for (l in seq_len(n - 2)) {
                        if (!(pag[mpath[l], mpath[l + 2]] == 0 &&
                            pag[mpath[l + 2], mpath[l]] == 0)) {
                            uncov <- FALSE
                            break ## speed up!
                        }
                    }
                    if (verbose) {
                        cat("    minUncovPdPath: path found: ", mpath, ", uncovered: ", uncov, "\n")
                    }
                    ## if it is uncovered
                    if (uncov) {
                        if (length(unfVect) == 0 || ## <<- normal version: just save
                            ## conservative version, check the path to be faithful:
                            faith.check(mpath, unfVect, p)) {
                            ## save the path to be returned
                            min.upd.path <- mpath
                            done <- TRUE
                        }
                    }
                } else {
                    ## d and c are either not connected or connected with a "wrong" edge -----> search iteratively
                    ## find all neighbours of d not visited yet
                    indR <- which((pag[d, ] == 1 | pag[d, ] == 2) &
                        (pag[, d] == 1 | pag[, d] == 3) & !visited)
                    if (length(indR) > 0) {
                        ## update the queues
                        path.list <- updateList(mpath, indR, path.list)
                    }
                }
            } ## {while}
        }
    }
    min.upd.path
} ## {minUncovPdPath}

## R5
minUncovCircPath <- function(p, pag, path, unfVect, verbose = FALSE) {
    ## Purpose: find a minimal uncovered circle path for a,c,d,b saved in path.
    ## Check also for the conservative case that it is unambiguous
    ## If a path exists this is the output, otherwise NA
    ## ----------------------------------------------------------------------
    ## Arguments: - p: number of nodes in the graph
    ##            - pag: adjacency matrix
    ##            - path: a,c,d,b under interest
    ##            - unfVect: vector containing the unfaithful triples
    ## ----------------------------------------------------------------------
    ## Author: Diego Colombo, Date: 19 Oct 2011, 13:11

    visited <- rep(FALSE, p)
    visited[path] <- TRUE # (a,b,c,d) all 'visited'
    a <- path[1]
    c <- path[2]
    d <- path[3]
    b <- path[4]
    min.ucp.path <- NA
    ## find all neighbours of c not visited yet
    indX <- which(pag[c, ] == 1 & pag[, c] == 1 & !visited) ## c o-o x
    if (length(indX) > 0) {
        path.list <- updateList(c, indX, NULL)
        done <- FALSE
        while (!done && length(path.list) > 0) {
            ## next element in the queue
            mpath <- path.list[[1]]
            x <- mpath[length(mpath)]
            path.list[[1]] <- NULL
            visited[x] <- TRUE
            if (pag[x, d] == 1 && pag[d, x] == 1) {
                ## circle path found
                mpath <- c(a, mpath, d, b)
                n <- length(mpath)
                ## check the path to be uncovered
                uncov <- TRUE
                for (l in seq_len(n - 2)) {
                    if (!(pag[mpath[l], mpath[l + 2]] == 0 &&
                        pag[mpath[l + 2], mpath[l]] == 0)) {
                        uncov <- FALSE
                        break ## speed up!
                    }
                }
                ## if it is uncovered
                if (uncov) {
                    if (length(unfVect) == 0 || ## <<- normal version: just save
                        ## conservative version, check the path to be faithful:
                        faith.check(mpath, unfVect, p)) {
                        ## save the path to be returned
                        min.ucp.path <- mpath
                        done <- TRUE
                    }
                }
            } else {
                ## x and d are either not connected or connected with an edge which is not o--o -----> search iteratively
                ## find all neighbours of x not visited yet
                indR <- which(pag[x, ] == 1 & pag[, x] == 1 & !visited) ## x o--o r
                if (length(indR) > 0) {
                    ## update the queues
                    path.list <- updateList(mpath, indR, path.list)
                }
            }
        } ## {while}
    }
    return(min.ucp.path)
}

## R4
minDiscrPath <- function(pag, a, b, c, verbose = FALSE) {
    ## Purpose: find a minimal discriminating path for a,b,c.
    ## If a path exists this is the output, otherwise NA
    ## ----------------------------------------------------------------------
    ## Arguments: - pag: adjacency matrix
    ##            - a,b,c: node positions under interest
    ## ----------------------------------------------------------------------
    ## Author: Diego Colombo, Date: 25 Jan 2011; speedup: Martin Maechler

    p <- as.numeric(dim(pag)[1])
    visited <- rep(FALSE, p)
    visited[c(a, b, c)] <- TRUE # {a,b,c} "visited"
    ## find all neighbours of a  not visited yet
    indD <- which(pag[a, ] != 0 & pag[, a] == 2 & !visited) ## d *-> a
    if (length(indD) > 0) {
        path.list <- updateList(a, indD, NULL)
        while (length(path.list) > 0) {
            ## next element in the queue
            mpath <- path.list[[1]]
            m <- length(mpath)
            d <- mpath[m]
            if (pag[c, d] == 0 & pag[d, c] == 0) {
                ## minimal discriminating path found :
                return(c(rev(mpath), b, c))
            }

            ## else :
            pred <- mpath[m - 1]
            path.list[[1]] <- NULL


            ## d is connected to c -----> search iteratively
            if (pag[d, c] == 2 && pag[c, d] == 3 && pag[pred, d] == 2) {
                visited[d] <- TRUE
                ## find all neighbours of d not visited yet
                indR <- which(pag[d, ] != 0 & pag[, d] == 2 & !visited) ## r *-> d
                if (length(indR) > 0) {
                    ## update the queues
                    path.list <- updateList(mpath[-1], indR, path.list)
                }
            }
        } ## {while}
    }
    ## nothing found:  return
    NA
} ## {minDiscrPath}

###################

updateList <- function(path, set, old.list) {
    ## Purpose: update the list of all paths in the iterative functions
    ## minDiscrPath, minUncovCircPath and minUncovPdPath
    ## ----------------------------------------------------------------------
    ## Arguments: - path: the path under investigation
    ##            - set: (integer) index set of variables to be added to path
    ##            - old.list: the list to update
    ## ----------------------------------------------------------------------
    ## Author: Diego Colombo, Date: 21 Oct 2011; Without for() by Martin Maechler
    c(old.list, lapply(set, function(s) c(path, s)))
}

## R9-R10
minUncovPdPath <- function(p, pag, a, b, c, unfVect, verbose = FALSE) {
    ## Purpose: find a minimal uncovered pd path for a,b,c saved in path.
    ## Check also for the conservative case that it is unambiguous
    ## If a path exists this is the output, otherwise NA
    ## ----------------------------------------------------------------------
    ## Arguments: - p: number of nodes in the graph
    ##            - pag: adjacency matrix
    ##            - a,b,c : nodes under interest
    ##            - unfVect: vector containing the ambiguous triples
    ## ----------------------------------------------------------------------
    ## Author: Diego Colombo, Date: 19 Oct 2011; small changes: Martin Maechler, Joris Mooij

    ## first check whether a,b,c is already a upd path
    stopifnot((pag[a, b] == 1 | pag[a, b] == 2) &
        (pag[b, a] == 1 | pag[b, a] == 3))
    min.upd.path <- NA
    done <- FALSE
    if ((pag[b, c] == 1 | pag[b, c] == 2) &
        (pag[c, b] == 1 | pag[c, b] == 3) &
        (pag[c, a] == 0)) {
        mpath <- c(a, b, c)
        if (length(unfVect) == 0 || ## <<- normal version: just save
            ## conservative version, check the path to be faithful:
            faith.check(mpath, unfVect, p)) {
            ## save the path to be returned
            min.upd.path <- mpath
            if (verbose) {
                cat("    minUncovPdPath: path found: ", mpath, ", uncovered: ", TRUE, "\n")
            }
            done <- TRUE
        }
    }

    ## now check paths of 4 or more nodes of the form <a,b,...,c>
    if (!done) {
        visited <- rep(FALSE, p)
        visited[c(a, b, c)] <- TRUE
        min.upd.path <- NA
        ## find all neighbours of b not visited yet
        indD <- which((pag[b, ] == 1 | pag[b, ] == 2) &
            (pag[, b] == 1 | pag[, b] == 3) &
            (pag[, a] == 0) & !visited)
        if (length(indD) > 0) {
            path.list <- updateList(b, indD, NULL)
            done <- FALSE
            while ((length(path.list) > 0) && (!done)) {
                ## next element in the queue
                mpath <- path.list[[1]]
                m <- length(mpath)
                d <- mpath[m]
                path.list[[1]] <- NULL
                visited[d] <- TRUE
                if (any(pag[d, c] == 1:2) && any(pag[c, d] == c(1, 3))) {
                    ## pd path found
                    mpath <- c(a, mpath, c)
                    n <- length(mpath)
                    ## check the path to be uncovered
                    uncov <- TRUE
                    for (l in seq_len(n - 2)) {
                        if (!(pag[mpath[l], mpath[l + 2]] == 0 &&
                            pag[mpath[l + 2], mpath[l]] == 0)) {
                            uncov <- FALSE
                            break ## speed up!
                        }
                    }
                    if (verbose) {
                        cat("    minUncovPdPath: path found: ", mpath, ", uncovered: ", uncov, "\n")
                    }
                    ## if it is uncovered
                    if (uncov) {
                        if (length(unfVect) == 0 || ## <<- normal version: just save
                            ## conservative version, check the path to be faithful:
                            faith.check(mpath, unfVect, p)) {
                            ## save the path to be returned
                            min.upd.path <- mpath
                            done <- TRUE
                        }
                    }
                } else {
                    ## d and c are either not connected or connected with a "wrong" edge -----> search iteratively
                    ## find all neighbours of d not visited yet
                    indR <- which((pag[d, ] == 1 | pag[d, ] == 2) &
                        (pag[, d] == 1 | pag[, d] == 3) & !visited)
                    if (length(indR) > 0) {
                        ## update the queues
                        path.list <- updateList(mpath, indR, path.list)
                    }
                }
            } ## {while}
        }
    }
    min.upd.path
} ## {minUncovPdPath}


#################
checkEdges <- function(suffStat, indepTest, alpha, apag, sepset, path,
                       unfVect = NULL, verbose = FALSE) {
    ## Purpose: check if every edge on the path should exist in R4
    ## ----------------------------------------------------------------------
    ## Values: - updated sepset and apag
    ##         - deleted==FALSE no edge has been deleted on the path
    ##                  ==TRUE the discriminating path doesn't exist anymore
    ## ----------------------------------------------------------------------
    ## Author: Diego Colombo, Date: 17 Aug 2010, 16:21

    stopifnot((n.path <- length(path)) >= 2)
    ## did we delete an edge?
    found <- FALSE
    ## conservative v-structures or not?
    conservative <- (length(unfVect) > 0)

    ## set that at the beginning there are no new unfaithful v-structures
    unfTripl <- NULL
    ## define the sepset
    SepSet.tot <- unique(c(
        sepset[[path[1]]][[path[n.path]]],
        sepset[[path[n.path]]][[path[1]]]
    ))
    if (length(SepSet.tot) != 0) {
        if (verbose) {
            cat(
                "\nCheck discriminating path:", path,
                "for dependence of any edge given sepset", SepSet.tot, "\n"
            )
        }
        p <- nrow(apag)
        ## check every edge on the path for independence given every possible subset of SepSet
        for (i in seq_len(n.path - 1)) {
            x <- path[i]
            y <- path[i + 1]
            SepSet <- setdiff(SepSet.tot, c(x, y))
            x. <- min(x, y)
            y. <- max(x, y)
            if (verbose >= 2) {
                cat("Edge: ", x, "*-*", y, "; Sepset=", SepSet, "; |S|=", length(SepSet), "\n")
            }
            if (length(SepSet) != 0) {
                j <- 0
                while (!found && j < length(SepSet)) {
                    j <- j + 1
                    ## all combinations of SepSet of size j
                    S.j <- if (j == 1 && length(SepSet) == 1) matrix(SepSet, 1, 1) else combn(SepSet, j)
                    ii <- 0
                    while (!found && ii < ncol(S.j)) {
                        ii <- ii + 1
                        pval <- indepTest(x, y, S.j[, ii], suffStat)
                        if (verbose) {
                            cat("x=", x, " y=", y, " S=", S.j[, ii], ": pval =", pval, "\n")
                        }
                        if (pval >= alpha) {
                            if (verbose) cat("Independence found: delete edge between", x, "and", y, "\n")
                            found <- TRUE
                            ## delete edge and save set in sepset
                            apag[x, y] <- apag[y, x] <- 0
                            sepset[[x]][[y]] <- sepset[[y]][[x]] <- S.j[, ii]
                            ## before we had a triangle and now it is an unshielded triple x-m-y
                            indM <- setdiff(
                                which(apag[x, ] != 0 & apag[, x] != 0 &
                                    apag[y, ] != 0 & apag[, y] != 0),
                                c(x, y)
                            ) ## just to be sure
                            ## create the list with all the new unshielded triples to be tested
                            if ((nI <- length(indM)) > 0) {
                                triplM <- matrix(integer(), 3, nI)
                                newVect <- numeric(nI)
                                for (jj in seq_len(nI)) {
                                    m <- indM[jj]
                                    triplM[, jj] <- c(x., m, y.)
                                    newVect[jj] <- triple2numb(p, x., m, y.)
                                }
                                ## new unshielded triple to be tested
                                r.v <- rfci.vStruc(suffStat, indepTest, alpha, sepset, apag,
                                    unshTripl = triplM, unshVect = newVect,
                                    conservative = conservative, verbose = verbose
                                )
                                ## save the modified graph g in apag
                                apag <- r.v$amat
                                ## save the new sepset
                                sepset <- r.v$sepset
                                ## save the new unfTripl, since we tested new v-structures and some can be unfaithful
                                ## for the conservative version
                                unfTripl <- r.v$unfTripl
                            }
                        }
                    } ## while(!found && i < *)
                } ## while(!found && j < *)
            }
        } ## for(i ....)
    }
    ## if SepSet is the empty set do nothing because surely the vertices are dependent
    list(deleted = found, apag = apag, sepset = sepset, unfTripl = unfTripl)
}
