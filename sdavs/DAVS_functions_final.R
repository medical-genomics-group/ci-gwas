####################################################################
####################################################################
#Functions for run DAVS_cor
####################################################################
####################################################################



Davs.con.causaleffect_cor_new<-function(datcor,w.pos,y.pos ,Q,pag,n,alpha=0.05,models="match",
                                    method = "subclass", subclass=6, type="ps",possdsep="small",
                                    max_k=NULL, max_depth=NULL, force_directed=TRUE){
  
  ###  dat:observational data; PCv: PCwy and candidate s
  ###  type      "iv" (for inverse variance) or "PS" (for propensity score)
  # require(pcalg)
  
  starttime<-proc.time()
  pag <- matrix(pag,nrow(pag),ncol(pag),dimnames = NULL)
  if(force_directed) {
    if(pag[w.pos,y.pos]==2 & pag[y.pos,w.pos]==2){
      pag[w.pos,y.pos]<-2
      pag[y.pos,w.pos]<-3
    }
  }
  
  # w.pos <- ncol(expdat) - 1 
  # y.pos <- ncol(expdat)
  pdes <- possibleDe(pag,w.pos) ## trying this out 11/27/2016
  #print(pdes)
  if (!(y.pos %in% pdes)) return(0)
  # possible d-sep
  if(possdsep=="small"){
    pdsep <- union(pdsepset.reach(w.pos, y.pos, -1, pag, max_depth=max_depth), y.pos) ##added y.pos to this set 8/29
  }
  if(possdsep=="big"){
    pdsep <- union(reach(w.pos,y.pos,-1,pag),y.pos) ##added y.pos to this set 8/29
  }
  pdsep <- sort(setdiff(pdsep,c(w.pos,y.pos)))
  
  pdsepset <- as.vector(setdiff(pdsep,c(Q,pdes)))
  
  valid_Z<-list()
  DAVS.ACE<-list()
  DAVS.sd<-list()
  
  
  colnames(datcor)<-paste("X", 1:dim(datcor)[1])
  rownames(datcor)<-paste("X", 1:dim(datcor)[1])
  #tempdat=as(tempdat,"matrix")
  #colnames(tempdat)[,2:(length(pdsepset)+1)]<-paste("X", 1:(length(pdsepset)), sep = "")
  
  # tempdata<-expdat[pdsepset]
  # VarQ <- expdat[Q]
  # # Tr<-expdat$Tr;Y<-expdat$Y
  # Tr<-expdat[,w.pos];Y<-expdat[,y.pos]
  # tempdat<-cbind(VarQ,tempdata, Tr,Y)
  tempdat<-datcor[c(Q,pdsepset, w.pos,y.pos),c(Q,pdsepset, w.pos,y.pos)]
  #print(colnames(tempdat))
  
  colnames(tempdat)[1]<-c("Q")
  rownames(tempdat)[1]<-c("Q")
  colnames(tempdat)[dim(tempdat)[1]]<-c("Y")
  colnames(tempdat)[dim(tempdat)[1]-1]<-c("Tr")
  rownames(tempdat)[dim(tempdat)[1]]<-c("Y")
  rownames(tempdat)[dim(tempdat)[1]-1]<-c("Tr")
  suffStat<-list(C=tempdat,n=n)
  ordnum <- 1; Q <- 1;nvar<-2:(length(pdsepset)+1)
  
  
  ###    Rule 2  #### 
  Pvaluewy<-gaussCItest(ncol(tempdat)-1,ncol(tempdat),c(),suffStat)
  Pvaluesw<-gaussCItest(Q,ncol(tempdat)-1,c(),suffStat)
  Pvaluesy<-gaussCItest(Q,ncol(tempdat),c(),suffStat)
  if(Pvaluewy > alpha || (Pvaluesw < alpha && Pvaluesy > alpha)){
    DAVS.ACE[[ordnum]]<-NULL;#DAVS.sd[[ordnum]]<-NULL
    ordnum <- ordnum+1
    # cat("The causal effect of W on Y does not be estimated from this dataset.","\n")
  }else{
    ##Rule 1  ####
    # create the size is equal to 1 Called C1(candidate set)  --- L1
    L<-combn(nvar,1)
    #print(L)
    L.tmp <- list() 
    for (ii in 1:ncol(L)) {
      Z<-L[,ii]
      #Print(Z)
      Pvalue01<-gaussCItest(Q,ncol(tempdat),Z,suffStat)
      # cat("P-VALUE01:",Pvalue01,"\n")
      if(Pvalue01<alpha){
        M<-c(ncol(tempdat)-1,Z) # W \cup Z
        Pvalue02<-gaussCItest(Q,ncol(tempdat),M,suffStat)
        if(Pvalue02>alpha){
          if(!(list(pdsepset[Z-1]) %fin% valid_Z)){
            cat("Q:",Q,"\n");cat("Z:",pdsepset[Z-1],"\n")
            valid_Z[[ordnum]]<-pdsepset[Z-1]
            L.tmp[[ordnum]]<-L[,ii]  
            if(models == "match"){
              rest<-est_match_con(Z,tempdat,type) 
            } 
            if(models == "logitreg"){
              rest<-est_reg_bin(Z,tempdat) 
            }
            if(models == "linearreg"){
              #print(Z)
              rest<-est_reg_con_cor(Z,tempdat) 
            }
            if(models == "matchit" && method == "subclass"){
              rest<-est_matchit_con(Z, tempdat, method = "subclass", subclass=subclass)
            }
            if(models == "matchit" && method == "nearest"){
              rest<-est_matchit_con(Z, tempdat, method = "nearest")
            }
            if(models == "matchit" && method == "cem"){
              rest<-est_matchit_con(Z, tempdat, method = "cem", subclass=subclass)
            }
            DAVS.ACE[[ordnum]]<-rest[1]
            DAVS.sd[[ordnum]]<-rest[2]
            ordnum<-ordnum+1 
          }
        }
      }
    }
    Fk<-setdiff(L, L.tmp)
    if(length(Fk)==2){
      Pvalue01<-gaussCItest(Q,ncol(tempdat),Fk,suffStat)
      if(Pvalue01<alpha){
        M<-c(ncol(tempdat)-1,Fk) # W \cup Z
        Pvalue02<-gaussCItest(Q,ncol(tempdat),Fk,suffStat) 
        if(Pvalue02>alpha){
          if(!(list(pdsepset[Z]) %fin% valid_Z)){
            valid_Z[[ordnum]]<-pdsepset[Fk] 
            if(models == "match"){
              rest<-est_match_con(Z,tempdat,type) 
            } 
            if(models == "logitreg"){
              rest<-est_reg_bin(Z,tempdat) 
            }
            if(models == "linearreg"){
              #print(Z)
              rest<-est_reg_con_cor(Z,tempdat) 
            }
            if(models == "matchit" && method == "subclass"){
              rest<-est_matchit_con(Z,tempdat, method="subclass",subclass=subclass)
            }
            if(models == "matchit" && method == "nearest"){
              rest<-est_matchit_con(Z, tempdat, method = "nearest")
            }
            if(models == "matchit" && method == "cem"){
              rest<-est_matchit_con(Z,tempdat, method="cem",subclass=subclass)
            }
            DAVS.ACE[[ordnum]]<-rest[1]
            DAVS.sd[[ordnum]]<-rest[2]
            ordnum<-ordnum+1 
          }
        }
      }
    }
    if(length(Fk)>2){
      for (k in 2:(length(Fk))) {
        if((!is.null(max_k)) && (k > max_k)) {
            break
        }
        Ck<-create_Ck(Fk,k-1)   # Candidate generation functiOn
        if(length(Ck)==0){
          break
        }
        Ltmp <- list()
        for (jj in 1:length(Ck)) {
          Z<-Ck[[jj]]
          Pvalue01<-gaussCItest(Q,ncol(tempdat),Z,suffStat)
          if(Pvalue01<alpha){
            M<-c(ncol(tempdat)-1,Z) # W \cup Z
            Pvalue02<-gaussCItest(Q,ncol(tempdat),M,suffStat) 
            if(Pvalue02>alpha){
              if(!(list(pdsepset[Z]) %fin% valid_Z)){
                valid_Z[[ordnum]]<-pdsepset[Z] 
                if(models == "match"){
                  rest<-est_match_con(Z,tempdat,type) 
                } 
                if(models == "logitreg"){
                  rest<-est_reg_bin(Z,tempdat) 
                }
                if(models == "linearreg"){
                  #print(Z)
                  rest<-est_reg_con_cor(Z,tempdat) 
                }
                if(models == "matchit" && method == "subclass"){
                  rest<-est_matchit_con(Z,tempdat,method="subclass",subclass=subclass)
                }
                if(models == "matchit" && method == "nearest"){
                  rest<-est_matchit_con(Z,tempdat,method="nearest")
                }
                if(models == "matchit" && method == "cem"){
                  rest<-est_matchit_con(Z,tempdat,method="cem",subclass=subclass)
                }
                DAVS.ACE[[ordnum]]<-rest[1]
                DAVS.sd[[ordnum]]<-rest[2]
                ordnum<-ordnum+1 
              }
            }else{
              Ltmp[[jj]]<-Ck[[jj]]
            }
          }
        }
        if(length(Ltmp)==0){
          break
        }
        Ltmp1 = Ltmp[-which(sapply(Ltmp, is.null))]
        if(length(Ltmp1)==0){
          Fk<-Ck 
        }else{
          Fk<-Ltmp1  
        }
        if(length(Fk)==1){
          break
        }
      }
    }
  }
  
  if(length(DAVS.ACE)==0 ){
    DAVS_ACE<-list(); ACE_DCE<-vector();SD_DCE<-vector() 
  }else{
    DAVS_ACE<-unlist(DAVS.ACE)
    ACE_DCE<-mean(DAVS_ACE)
    DAVS_sd<-unlist(DAVS.sd)
    SD_DCE<-mean(DAVS_sd)
  }
  runtime <- proc.time()-starttime
  Runtime<-runtime[1]
  
  #retu<-list(DAVS_ACE,ACE_DCE,SD_DCE,valid_Z,Runtime)
  #names(retu) <- c("DAVS_ACE","ACE_DCE","SD_DCE","valid_Z","Runtime")
  if(is.logical(ACE_DCE) ) ACE_DCE=0
  #print(ACE_DCE)
  retu<-(ACE_DCE)
  #names(retu) <- c("ACE_DCE")
  return(retu)
}



########################
#functions for running DAVS for  max_k, max_depth 
#######################

pdsepset.reach <- function(a, b, c, adjacency, max_depth=NULL)
{
  ## reachable      set of vertices;
  ## edgeslist      array[1..maxvertex] of list of edges
  ## numvertex      integer
  ## labeled        array (by depth) of list of edges that have been labeled
  
  makeedge <- function(x,y) list(list(x,y))
  
  legal.dsep <- function(r,s) {
    ## Modifying global 'edgeslist'
    if (((adjacency[r[[1]],r[[2]]] == 2 &&
          adjacency[s,     r[[2]]] == 2 && r[[1]] != s) || ((adjacency[r[[1]],s] != 0 && r[[1]] != s))) &&  (is.poss.ancestor(s,a,adjacency) || is.poss.ancestor(s,b,adjacency))    && (is.poss.ancestor(r[[2]],a,adjacency) || is.poss.ancestor(r[[2]],b,adjacency))) {
      edgeslist[[r[[2]]]] <<- setdiff(edgeslist[[r[[2]]]],s)
      makeedge(r[[2]],s)
    }
  }
  
  
  initialize.pdsep <- function(x,y) mapply(makeedge, x=x, y=y)
  
  labeled <- list()
  numvertex <- dim(adjacency)[1]
  edgeslist <- list()
  
  for (i in 1:numvertex)
    edgeslist <- c(edgeslist,list(which(adjacency[,i] != 0)))
  #print( edgeslist)
  labeled[[1]] <- initialize.pdsep(a, edgeslist[[a]])
  #print(labeled)
  edgeslist[[a]] <- list()
  depth <- 2
  repeat {
    labeled[[depth]] <- list()
    for (i in seq_along(labeled[[depth-1]])) {
      lab.i <- labeled[[depth-1]][[i]]
      edgestemp <- edgeslist[[lab.i[[2]]]]
      if (length(edgestemp) == 0) break
      for (j in seq_along(edgestemp))
        labeled[[depth]] <- union(legal.dsep(lab.i, edgestemp[[j]]),
                                  labeled[[depth]])
    }
    if (length(labeled[[depth]]) == 0)
      break
    ## else :
    depth <- depth  + 1
    if ((!is.null(max_depth)) && (depth >= max_depth))
      break
  }
  dsep <- unique(unlist(labeled))
  dsep
} # end function






#########################
pdsepset.reach_orginal <- function(a,b,c,adjacency)
{
  ## reachable      set of vertices;
  ## edgeslist      array[1..maxvertex] of list of edges
  ## numvertex      integer
  ## labeled        array (by depth) of list of edges that have been labeled
  
  makeedge <- function(x,y) list(list(x,y))
  
  legal.dsep <- function(r,s) {
    ## Modifying global 'edgeslist'
    if (((adjacency[r[[1]],r[[2]]] == 2 &&
          adjacency[s,     r[[2]]] == 2 && r[[1]] != s) || ((adjacency[r[[1]],s] != 0 && r[[1]] != s))) &&  (is.poss.ancestor(s,a,adjacency) || is.poss.ancestor(s,b,adjacency))    && (is.poss.ancestor(r[[2]],a,adjacency) || is.poss.ancestor(r[[2]],b,adjacency))) {
      edgeslist[[r[[2]]]] <<- setdiff(edgeslist[[r[[2]]]],s)
      makeedge(r[[2]],s)
    }
  }
  
  
  initialize.pdsep <- function(x,y) mapply(makeedge, x=x, y=y)
  
  labeled <- list()
  numvertex <- dim(adjacency)[1]
  edgeslist <- list()
  
  for (i in 1:numvertex)
    edgeslist <- c(edgeslist,list(which(adjacency[,i] != 0)))
  #print( edgeslist)
  labeled[[1]] <- initialize.pdsep(a, edgeslist[[a]])
  #print(labeled)
  edgeslist[[a]] <- list()
  depth <- 2
  repeat {
    labeled[[depth]] <- list()
    for (i in seq_along(labeled[[depth-1]])) {
      lab.i <- labeled[[depth-1]][[i]]
      edgestemp <- edgeslist[[lab.i[[2]]]]
      if (length(edgestemp) == 0) break
      for (j in seq_along(edgestemp))
        labeled[[depth]] <- union(legal.dsep(lab.i, edgestemp[[j]]),
                                  labeled[[depth]])
    }
    if (length(labeled[[depth]]) == 0)
      break
    ## else :
    depth <- depth  + 1
   # if (depth >= max_depth) 
   #   break
  }
  dsep <- unique(unlist(labeled))
  dsep
} # end function





reach <- function(a,b,c,adjacency)
{
  ## reachable      set of vertices;
  ## edgeslist      array[1..maxvertex] of list of edges
  ## numvertex      integer
  ## labeled        array (by depth) of list of edges that have been labeled
  
  makeedge <- function(x,y) list(list(x,y))
  
  legal.pdsep <- function(r,s) {
    ## Modifying global 'edgeslist'
    if ((adjacency[r[[1]],r[[2]]] == 2 &&
         adjacency[s,     r[[2]]] == 2 && r[[1]] != s) ||
        (adjacency[r[[1]],s] != 0 && r[[1]] != s)) {
      edgeslist[[r[[2]]]] <<- setdiff(edgeslist[[r[[2]]]],s)
      makeedge(r[[2]],s)
    }
  }
  
  initialize.pdsep <- function(x,y) mapply(makeedge, x=x, y=y)
  
  labeled <- list()
  numvertex <- dim(adjacency)[1]
  edgeslist <- list()
  for (i in 1:numvertex)
    edgeslist <- c(edgeslist,list(which(adjacency[,i] != 0)))
  labeled[[1]] <- initialize.pdsep(a, edgeslist[[a]])
  edgeslist[[a]] <- list()
  depth <- 2
  repeat {
    labeled[[depth]] <- list()
    for (i in seq_along(labeled[[depth-1]])) {
      lab.i <- labeled[[depth-1]][[i]]
      edgestemp <- edgeslist[[lab.i[[2]]]]
      if (length(edgestemp) == 0) break
      for (j in seq_along(edgestemp))
        labeled[[depth]] <- union(legal.pdsep(lab.i, edgestemp[[j]]),
                                  labeled[[depth]])
    }
    if (length(labeled[[depth]]) == 0)
      break
    ## else :
    depth <- depth  + 1
  }
  unique(unlist(labeled))
}

possibleDe <- function(amat,x)
{
  ## Purpose: in a DAG, CPDAG, MAG, or PAG determine which nodes are
  ##          possible descendants of x on definite status paths
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - amat: matrix corresponding to the DAG, CPDAG, MAG, or PAG
  ## - x: node of interest
  ## ----------------------------------------------------------------------
  ## Value:
  ## - de.list: array containing the possible descendants of x
  ## ----------------------------------------------------------------------
  ## Author: Diego Colombo, Date: 26 Apr 2012, 16:58
  
  stopifnot(is.matrix(amat))
  p <- nrow(amat)
  is.de <- rep.int(FALSE, p) ##
  ## 1. case: x is a possible child of itself
  is.de[x] <- TRUE
  ## 2. case: find all the possible children of x
  indD <- which(amat[x,] != 0  & amat[,x] != 2 & !is.de) ## x (o,-)-* d
  i.pr <- rep(x,length(indD))
  while (length(indD) > 0) {
    ##next element in the queue
    d <- indD[1]
    indD <- indD[-1]
    pred <- i.pr[1]
    i.pr <- i.pr[-1]
    is.de[d] <- TRUE
    a.d <- amat[,d]
    a.d.p <- a.d[pred]
    ## find all possible children of d not visited yet
    indR <- which(amat[d,] != 0 & a.d != 2 & !is.de) ## d (o,-)-* r
    for(j in seq_along(indR)) {
      ## check that the triple <pred,d,r> is of a definite status
      ## 1. d is a collider on this subpath; this is impossible
      ##    because the edge between d and r cannot be into d
      ## 2. d is a definite non-collider
      r <- indR[j]
      if (a.d.p == 3 || a.d[r] == 3 ||
          (a.d.p == 1 && a.d[r] == 1 && amat[pred,r] == 0)) {
        ## update the queues
        indD <- c(indD, r)
        i.pr <- c(i.pr, d)
      }
    }
  }
  ## return 'de.list' :
  which(is.de)
  
} ## {possibleDe}

is.poss.ancestor <- function(a, b, g,visited=NULL){
  if(a==b) return(TRUE)
  foundpath <- c()
  ind1 <- which(g==3, arr.ind=TRUE, useNames=FALSE) #tails
  ind11 <- which(g==1, arr.ind=TRUE, useNames=FALSE) #circles
  ind1 <- rbind(ind1,ind11) ## tails and circles
  ind1 <- subset(ind1,ind1[,2]==a) # array of indices for tails and circles out of A
  if(nrow(ind1)==0) return(FALSE) # addition: if there are no tails or circles at A
  for(x in 1:nrow(ind1)){ # loop through tails and circles out of A
    if(ind1[x,1] %in% visited) next
    if(g[ind1[x,2],ind1[x,1]]==2 || g[ind1[x,2],ind1[x,1]]==1){ # if there is an arrowhead or circle at the other end of the x-th tail (call this C)
      if(ind1[x,1]==b){
        foundpath <- append(foundpath,TRUE)
        break
      }
      if(any(g[,ind1[x,1]]==3 | g[,ind1[x,1]]==1)){ # if there are any tails or circles out of C
        a_old <- a
        a2 <- ind1[x,1]
        if(a2==a_old) next
        foundpath <- append(foundpath,is.poss.ancestor(a2,b,g,visited=c(visited,a_old)))
        if(any(foundpath)==TRUE) break
      }
    } # if there isn't an arrowhead at C - !(A-->C) - don't return anything
  } # for x in 1:nrow(ind1)
  if(any(foundpath)==TRUE) return(TRUE)
  else return(FALSE)
} # end function





standardize_data <- function(dat)
{
  N <- nrow(dat)
  sd_dat <- apply(dat, 2, sd)
  m_dat <- colMeans(dat)
  dat <- (dat - t(matrix(rep(m_dat, N), nrow = length(m_dat)))) / t(matrix(rep(sd_dat, N), nrow = length(sd_dat)))
  return(dat)
}



################
lm.cov<- function (C, y, x) {
 # print(x)
 solve(C[x, x], C[x, y, drop = FALSE])[1, ]
}
################
est_reg_con_cor<- function(adjset, cor.data) {
  # adjset    selected adjustment set (logical vector)
  # obs       correlation matrix
  
  covs <- colnames(cor.data)[1:(ncol(cor.data)-2)]
  
  adjcov <- covs[adjset]
  
  
  x <- c( "Tr", adjcov)
  y <- "Y"
  # C <- obs
  fit <- lm.cov(cor.data, y, x)
  #print(fit)
  est_ATE <- fit[1]
  print(est_ATE)
  #est_sd <- sqrt((Cyy - t(Cyx) %*% est_ATE) / Cxx)
  est_sd <- (fit[1])
  #}
  return(c(est_ATE, est_sd))
}


create_Ck<-function(Lksub1,k){
  countn<-1
  Ck<-list()
  len_Lksub1<-length(Lksub1)
  for (i in 1:(len_Lksub1-1)) {
    for (j in (i+1):len_Lksub1) {
      L1<-Lksub1[[i]]
      L2<-Lksub1[[j]]
      if(all((L1==L2)[-k])==TRUE & (L1==L2)[k]==FALSE){
        Ck_item<-union(L1,L2)
        temps <- combn(Ck_item, k,simplify = FALSE)
        require("fastmatch")
        if(all(temps %fin% Lksub1) == TRUE){
          Ck[[countn]]<- Ck_item
          countn<-countn+1
        }
      }
    }
  }
  return(Ck)
}



########################################################################
#####Cheng's estimators for estimating causal effects ##################
####################################  ###############################

est_reg_con<-function(adjset,obs){
  # adjset    selected adjustment set (logical vector)
  # obs       data frame with observation data
  covs <- names(obs)[1:(ncol(obs)-2)]
  if (anyNA(adjset)) {
    est_ATE <- NA
    est_sd <- NA
  } else { 
    adjcov <- covs[adjset]
    if (length(adjcov)==0) {form <- formula("Y~Tr")} else {
      form <- formula(paste("Y~Tr", paste(adjcov, collapse="+"), sep="+"))}
    fit <- lm(form, data=obs)
    est_ATE <- summary(fit)$coefficients["Tr", "Estimate"]
    est_sd <- summary(fit)$coefficients["Tr", "Std. Error"]
  }
  return(c(est_ATE, est_sd))
}
est_reg_bin <- function(adjset, obs){
  # adjset    selected adjustment set (logical vector)
  # obs       data frame with observation data
  
  covs <- names(obs)[1:(ncol(obs)-2)]
  if (anyNA(adjset)) {
    return(NA)
  } else { 
    adjcov <- covs[adjset]
    if (length(adjcov)==0) {form <- formula("Y~Tr")} else {
      form <- formula(paste("Y~Tr", paste(adjcov, collapse="+"), sep="+"))}
    fit <- glm(form, data=obs, family="binomial")
    stdreg <- tryCatch(stdGlm(fit, obs, X="Tr"), error=function(x){NA})
    if (anyNA(stdreg)) {return(NA)} else {
      p_0 <- stdreg$est[1]
      p_1 <- stdreg$est[2]
      est_logMCOR <- log(p_1/(1-p_1) / (p_0/(1-p_0)))
    }
    return(est_logMCOR)
  }
}


est_matchit_con<-function(adjset,obs, method = "subclasss", subclass=6){
  covs <- names(obs)[1:(ncol(obs)-2)]
  if (anyNA(adjset)) {
    est_ATE <- NA
    est_sd <- NA
  } else if (sum(adjset)==0) { 
    est_ATE <- mean(obs$Y[obs$Tr==1]) - mean(obs$Y[obs$Tr==0])
    n_1 <- sum(obs$Tr==1)
    n_0 <- sum(obs$Tr==0)
    s2_pooled <- (var(obs$Y[obs$Tr==1])*(n_1-1) +
                    var(obs$Y[obs$Tr==0])*(n_1-0)) / (n_1+n_0-2)
    est_sd <- sqrt(( 1/n_1 + 1/n_0 ) * s2_pooled)
  } else {
    require(Zelig)
    if(method == "subclass"){
      fmla1<-as.formula(paste("Tr~", paste(covs[adjset],collapse = "+")))
      txnam<-c("Tr", covs[adjset])
      fmla2<-as.formula(paste("Y~", paste(txnam, collapse = "+")))
      m.out <- matchit(fmla1, data = obs, method = "subclass", subclass = subclass)
      # learest square regression for con outcome
      zmodel.out<-zelig(fmla2, data=match.data(m.out), model = "ls", cite = FALSE, by="subclass")
      control.out<-setx(zmodel.out,Tr=0) #control
      treat.out<-setx(zmodel.out,Tr=1) #Tread
      s.out<-sim(zmodel.out, x = control.out, x1 = treat.out)          
      est_ATE <- mean(unlist(s.out$sim.out$x1$fd))
      est_sd <- sd(unlist(s.out$sim.out$x1$fd))  
    }
    if(method == "nearest"){
      fmla1<-as.formula(paste("Tr~", paste(covs[adjset],collapse = "+")))
      txnam<-c("Tr", covs[adjset])
      fmla2<-as.formula(paste("Y~", paste(txnam, collapse = "+")))
      m.out <- matchit(fmla1, data = obs, method = "nearest")
      z.out <- zelig(fmla2, data =match.data(m.out), model = "ls",cite = FALSE)
      x.out <- setx(z.out, Tr = 0)
      x1.out <- setx(z.out, Tr = 1)
      s.out <- sim(z.out, x = x.out, x1 = x1.out)
      
      est_ATE <- mean(unlist(s.out$sim.out$x1$fd))
      est_sd <- sd(unlist(s.out$sim.out$x1$fd)) 
    }
    if(method == "cem"){
      require(cem)
      fmla1<-as.formula(paste("Tr~", paste(covs[adjset], collapse = "+")))
      txnam<-c("Tr", covs[adjset])
      fmla2<-as.formula(paste("Y~", paste(txnam, collapse = "+")))
      glm1<-glm(fmla1,family = binomial,  data = obs)
      m.out<-matchit(glm1, method = "cem",data = obs)
      z.out <- zelig(fmla2, data =match.data(m.out), model = "ls",cite = FALSE)
      x.out <- setx(z.out, Tr = 0)
      x1.out <- setx(z.out, Tr = 1)
      s.out <- sim(z.out, x = x.out, x1 = x1.out)
      
      est_ATE <- mean(unlist(s.out$sim.out$x1$fd))
      est_sd <- sd(unlist(s.out$sim.out$x1$fd))  
    }
  }
  
  return(c(est_ATE, est_sd))
}
est_matchit_bin<-function(adjset,obs, method = "subclasss", subclass=6){
  covs <- names(obs)[1:(ncol(obs)-2)]
  if (anyNA(adjset)) {
    est_ATE <- NA
    est_sd <- NA
  } else if (sum(adjset)==0) { 
    est_ATE <- mean(obs$Y[obs$Tr==1]) - mean(obs$Y[obs$Tr==0])
    n_1 <- sum(obs$Tr==1)
    n_0 <- sum(obs$Tr==0)
    s2_pooled <- (var(obs$Y[obs$Tr==1])*(n_1-1) +
                    var(obs$Y[obs$Tr==0])*(n_1-0)) / (n_1+n_0-2)
    est_sd <- sqrt(( 1/n_1 + 1/n_0 ) * s2_pooled)
  } else {
    require(Zelig)
    if(method == "subclass"){
      fmla1<-as.formula(paste("Tr~", paste(covs[adjset],collapse = "+")))
      txnam<-c("Tr", covs[adjset])
      fmla2<-as.formula(paste("Y~", paste(txnam, collapse = "+")))
      m.out <- matchit(fmla1, data = obs, method = "subclass", subclass = subclass)
      # learest square regression for con outcome
      zmodel.out<-zelig(fmla2, data=match.data(m.out), model = "logit", cite = FALSE, by="subclass")
      control.out<-setx(zmodel.out,Tr=0) #control
      treat.out<-setx(zmodel.out,Tr=1) #Tread
      s.out<-sim(zmodel.out, x = control.out, x1 = treat.out)          
      est_ATE <- mean(unlist(s.out$sim.out$x1$fd))
      est_sd <- sd(unlist(s.out$sim.out$x1$fd))  
    }
    if(method == "nearest"){
      fmla1<-as.formula(paste("Tr~", paste(covs[adjset],collapse = "+")))
      txnam<-c("Tr", covs[adjset])
      fmla2<-as.formula(paste("Y~", paste(txnam, collapse = "+")))
      m.out <- matchit(fmla1, data = obs, method = "nearest")
      z.out <- zelig(fmla2, data =match.data(m.out), model = "logit",cite = FALSE)
      x.out <- setx(z.out, Tr = 0)
      x1.out <- setx(z.out, Tr = 1)
      s.out <- sim(z.out, x = x.out, x1 = x1.out)
      
      est_ATE <- mean(unlist(s.out$sim.out$x1$fd))
      est_sd <- sd(unlist(s.out$sim.out$x1$fd)) 
    }
    if(method == "cem"){
      require(cem)
      fmla1<-as.formula(paste("Tr~", paste(covs[adjset], collapse = "+")))
      txnam<-c("Tr", covs[adjset])
      fmla2<-as.formula(paste("Y~", paste(txnam, collapse = "+")))
      glm1<-glm(fmla1,family = binomial,  data = obs)
      m.out<-matchit(glm1, method = "cem",data = obs)
      z.out <- zelig(fmla2, data =match.data(m.out), model = "logit",cite = FALSE)
      x.out <- setx(z.out, Tr = 0)
      x1.out <- setx(z.out, Tr = 1)
      s.out <- sim(z.out, x = x.out, x1 = x1.out)
      
      est_ATE <- mean(unlist(s.out$sim.out$x1$fd))
      est_sd <- sd(unlist(s.out$sim.out$x1$fd))  
    }
  }
  
  return(c(est_ATE, est_sd))
}



est_match_con <- function(adjset, obs, type){
  # adjset    selected adjustment set (logical vector)
  # obs       data frame with observation data
  # type      "iv" (for inverse variance) or "PS" (for propensity score)
  
  # if(!(exists("run"))){run <<- 0}
  # run <<- run+1
  # cat(run, "\n")
  
  covs <- names(obs)[1:(ncol(obs)-2)]
  if (anyNA(adjset)) {
    est_ATE <- NA
    est_sd <- NA
  } else if (sum(adjset)==0) { 
    est_ATE <- mean(obs$Y[obs$Tr==1]) - mean(obs$Y[obs$Tr==0])
    n_1 <- sum(obs$Tr==1)
    n_0 <- sum(obs$Tr==0)
    s2_pooled <- (var(obs$Y[obs$Tr==1])*(n_1-1) +
                    var(obs$Y[obs$Tr==0])*(n_1-0)) / (n_1+n_0-2)
    est_sd <- sqrt(( 1/n_1 + 1/n_0 ) * s2_pooled)
  } else {
    if (type=="PS") {
      if (sum(adjset)==0) {PSform <- formula("Tr~1")} else {
        PSform <- formula(paste("Tr~", paste(covs[adjset], collapse="+"),
                                sep="")) }
      PSmodel <- glm(PSform, family="binomial", data=obs)
      converged_PS <- PSmodel$converged
      PS <- predict(PSmodel, type="response")  
      covadj <- PS
    } else {
      all <- obs[ ,1:(ncol(obs)-2)]
      covadj <- all[ ,adjset]
    }
    Mres <- Match(Tr=obs$Tr, Y=obs$Y, X=covadj, estimand="ATE", M=1,replace=FALSE)
    est_ATE <- Mres$est
    est_sd <- Mres$se
  }
  return(c(est_ATE, est_sd))
}


est_match_bin <- function(adjset, obs, type){
  # adjset    selected adjustment set (logical vector)
  # obs       data frame with observation data
  # type      "iv" (for inverse variance) or "PS" (for propensity score)
  
  # if(!(exists("run"))){run <<- 0}
  # run <<- run+1
  # cat(run, "\n")
  
  covs <- names(obs)[1:(ncol(obs)-2)]
  if (anyNA(adjset)) {
    est_lOR <- NA
  } else if (sum(adjset)==0) {
    mu_0 <- mean(obs$Y[obs$Tr==0])
    mu_1 <- mean(obs$Y[obs$Tr==1])
    est_lOR <- log( mu_1/(1-mu_1) / (mu_0/(1-mu_0)) )
  } else {
    if (type=="PS") {
      PSform <- formula(paste("Tr~", paste(covs[adjset], collapse="+"),
                              sep=""))
      PSmodel <- glm(PSform, family="binomial", data=obs)
      converged_PS <- PSmodel$converged
      PS <- predict(PSmodel, type="response")
      names(PS) <- NULL
      covadj <- as.matrix(PS)
    } else {
      all <- obs[ ,1:(ncol(obs)-2)]
      covadj <- all[ ,adjset]
    }
    
    Mres <- tryCatch({Match(Tr=obs$Tr, Y=obs$Y, X=covadj, estimand="ATE",
                            Weight=1)},
                     error=function(e){return(NA)}
    )
    if (anyNA(Mres)) {est_lOR <- NA} else {
      wei <- Mres$weights
      m_Y <- Mres$mdata$Y*wei
      m_Tr <- Mres$mdata$Tr
      m_n <- Mres$orig.nobs
      mu_0 <- sum(m_Y[m_Tr==0])/m_n
      mu_1 <- sum(m_Y[m_Tr==1])/m_n
      est_lOR <- log( mu_1/(1-mu_1) / (mu_0/(1-mu_0)) )
    }
  }
  return(est_lOR)
}

est_dr_con <- function(adjset, obs){
  # adjset    selected adjustment set (logical vector)
  # obs       data frame with observation data
  # 
  #   if(!(exists("run"))){run <<- 0}
  #   run <<- run+1
  #   cat(run, "\n")
  
  covs <- names(obs)[1:(ncol(obs)-2)]
  if (anyNA(adjset)) {
    est_ATE <- NA
    est_sd <- NA
  } else { 
    adjcov <- covs[adjset]
    if (length(adjcov)==0) {
      form_Tr <- formula("Tr~1")
      form_Y <- formula("Y~1")
    } else {
      form_Tr <- formula(paste("Tr~", paste(adjcov, collapse="+"), sep="+"))
      form_Y <- formula(paste("Y~", paste(adjcov, collapse="+"), sep="+"))
    }
    ppi.glm <- glm(form_Tr, data=obs, family=binomial)
    X <- model.matrix(ppi.glm)
    ppi.hat <- ppi.glm$fitted
    
    eta1.glm <- glm(form_Y, subset=Tr==1, data=obs)
    eta1.hat <- predict.glm(eta1.glm, type="response", newdata=obs)
    eta0.glm <- glm(form_Y, subset=Tr==0, data=obs)
    eta0.hat <- predict.glm(eta0.glm, type="response", newdata=obs)
    
    #ppi.hat treated as known
    out.lik <- tryCatch({ate.clik(obs$Y, obs$Tr, ppi.hat, g0=cbind(1,eta0.hat),
                                  g1=cbind(1,eta1.hat))}, error=function(e){return(NA)}
    )
    if (anyNA(out.lik)) {est_ATE <- est_sd <- NA} else {
      est_ATE <- out.lik$diff
      est_sd <- sqrt(out.lik$v.diff)
    }
  }
  return(c(est_ATE, est_sd))
}


est_dr_bin <- function(adjset, obs){
  # adjset    selected adjustment set (logical vector)
  # obs       data frame with observation data
  
  # if(!(exists("run"))){run <<- 0}
  # run <<- run+1
  # cat(run, "\n")
  
  covs <- names(obs)[1:(ncol(obs)-2)]
  if (anyNA(adjset)) {
    est_logMCOR <- NA
  } else { 
    adjcov <- covs[adjset]
    if (length(adjcov)==0) {
      form_Tr <- formula("Tr~1")
      form_Y <- formula("Y~1")
    } else {
      form_Tr <- formula(paste("Tr~", paste(adjcov, collapse="+"), sep="+"))
      form_Y <- formula(paste("Y~", paste(adjcov, collapse="+"), sep="+"))
    }
    ppi.glm <- glm(form_Tr, data=obs, family=binomial)
    X <- model.matrix(ppi.glm)
    ppi.hat <- ppi.glm$fitted
    
    eta1.glm <- glm(form_Y, subset=Tr==1, data=obs, family=binomial)
    eta1.hat <- predict.glm(eta1.glm, type="response", newdata=obs)
    eta0.glm <- glm(form_Y, subset=Tr==0, data=obs, family=binomial)
    eta0.hat <- predict.glm(eta0.glm, type="response", newdata=obs)
    
    #ppi.hat treated as known
    out.lik <- tryCatch({ate.clik(obs$Y, obs$Tr, ppi.hat, g0=cbind(1,eta0.hat),
                                  g1=cbind(1,eta1.hat))}, error=function(e){return(NA)}
    )
    if (anyNA(out.lik)) {est_logMCOR <- NA} else {
      p_0 <- out.lik$mu[2]
      p_1 <- out.lik$mu[1]
      est_logMCOR <- log(p_1/(1-p_1) / (p_0/(1-p_0)))
    }
  }
}
