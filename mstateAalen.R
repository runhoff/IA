require(timereg)
require(data.table)

################################################################################
##
##  Functions for additive hazards multi-state modeling in a similar manner to
##  mstate-package in R by Putter, de Wreede and Fiocco.
## 
################################################################################


##
##  mstateaalen:
##    Description: Fit additive hazards models for transitions in a multi-state model
##    Arguments: msmdata    - expanded dataset as in mstate with 1 row per possible
##                          transition, status 1 indicates the actual transition.
##                          Needs the columns Tstart (day number), Tstop, status
##                          trans (number of transition type according to tmat)
##                          and a unique unit identifier called id
##               msmformula - formula on the form Surv(Tstart,Tstop,status) ~ covariates
##               tmat       - Transition matrix describing the states and transitions
##                            in the multi-state model.
##               See aalen() in timereg-package for additional arguments
mstateAalen <- function(msmformula, msmdata, tmat, start.time = 0, 
                        max.time = NULL, robust = 0, id = NULL, clusters = NULL, 
                        residuals = 0, n.sim = 1000, weighted.test = 0, covariance = 0, 
                        resample.iid = 0, deltaweight = 1, silent = 1, weights = NULL, 
                        max.clust = 1000, gamma = NULL, offsets = 0, caseweight = NULL) {
  
  ntrans <- length(sort(unique(msmdata$trans)))
  aalen.list <- vector(mode = "list", length = ntrans)
  for(k in 1:ntrans) {
    
    aalen.list[[k]] <- aalen(formula = msmformula, data = msmdata[msmdata$trans == k,],
                             start.time, 
                             max.time, robust, id, clusters, 
                             residuals, n.sim, weighted.test, covariance, 
                             resample.iid, deltaweight, silent, weights, 
                             max.clust, gamma, offsets, caseweight)
    #Sys.sleep(1)
    #cat("\r",k, "/", ntrans, "models fitted")
  }
  attr(aalen.list, "tvec") <- sort(unique(msmdata$Tstop))
  attr(aalen.list, "tmat") <- tmat
  attr(aalen.list, "start.time") <- start.time
  class(aalen.list) <- "aalen.list"
  return(aalen.list)
}

## Weighted version for doing IPT weighted analysis. Requires an additional
## collumn in msdata called ws.
mstateAalenWs <- function(msmformula, msmdata, tmat, start.time = 0, 
                         max.time = NULL, robust = 0, id = NULL, clusters = NULL, 
                         residuals = 0, n.sim = 1000, weighted.test = 0, covariance = 0, 
                         resample.iid = 0, deltaweight = 1, silent = 1, weights = NULL, 
                         max.clust = 1000, gamma = NULL, offsets = 0, caseweight = NULL) {
  
  ntrans <- length(sort(unique(msmdata$trans)))
  aalen.list <- vector(mode = "list", length = ntrans)
  for(k in 1:ntrans) {
    
    aalen.list[[k]] <- aalen(formula = msmformula, data = msmdata[msmdata$trans == k,],
                             start.time, 
                             max.time, robust, id, clusters, 
                             residuals, n.sim, weighted.test, covariance, 
                             resample.iid, deltaweight, silent,
                             weights = msmdata[msmdata$trans == k,]$ws, 
                             max.clust, gamma, offsets, caseweight)
    #Sys.sleep(1)
    #cat("\r",k, "/", ntrans, "models fitted")
  }
  attr(aalen.list, "tvec") <- sort(unique(msmdata$Tstop))
  attr(aalen.list, "tmat") <- tmat
  attr(aalen.list, "start.time") <- start.time
  class(aalen.list) <- "aalen.list"
  return(aalen.list)
}

## For bootstrap purposes. Needs extra tinkering to work atm.
mstateAalenBOOT <- function(msmformula, msmdata, tmat, start.time = 0, 
                          max.time = NULL, robust = 0, id = NULL, clusters = NULL, 
                          residuals = 0, n.sim = 1000, weighted.test = 0, covariance = 0, 
                          resample.iid = 0, deltaweight = 1, silent = 1, weights = NULL, 
                          max.clust = 1000, gamma = NULL, offsets = 0, caseweight = NULL) {
  
  ntrans <- length(sort(unique(msmdata$trans)))
  aalen.list <- vector(mode = "list", length = ntrans)
  for(k in 1:ntrans) {
    
    aalen.list[[k]] <- aalen(formula = msmformula, data = msmdata[msmdata$trans == k,],
                             start.time, 
                             max.time, robust, id, clusters, 
                             residuals, n.sim, weighted.test, covariance, 
                             resample.iid, deltaweight, silent,
                             weights = msmdata[msmdata$trans == k,]$bootW, 
                             max.clust, gamma, offsets, caseweight)
    #Sys.sleep(1)
    #cat("\r",k, "/", ntrans, "models fitted")
  }
  attr(aalen.list, "tvec") <- sort(unique(msmdata$Tstop))
  attr(aalen.list, "tmat") <- tmat
  attr(aalen.list, "start.time") <- start.time
  class(aalen.list) <- "aalen.list"
  return(aalen.list)
}

## Calculates cumulative transitions hazards for specified covariate vector x
## in the multi-state model stored in object 
msfitAalen <- function(object, x) {
  ntrans <- length(object)
  dt.list <-  vector(mode = "list", length = ntrans)
  n.jumps <- rep(0, ntrans)
  for(k in 1:ntrans) {
    dt.list[[k]] <- data.table(object[[k]]$cum)
    setkey(dt.list[[k]], "time")
  }
  A.list <- dt.list
  dt.time <- data.table(time = c(attr(object, "start.time"), attr(object,"tvec")), key = "time")
  for(k in 1:ntrans) {
    A.list[[k]] <- merge(dt.time, A.list[[k]], all.x = TRUE)
    A.list[[k]] <- na.locf(A.list[[k]], na.rm = FALSE)
    A.list[[k]][is.na(A.list[[k]])] <- 0
  }
  AA <- matrix(nrow = dim(dt.time)[1], ncol = ntrans + 1)
  AA[,1] <- dt.time$time
  for(k in 1:ntrans) {
    test <- as.matrix(A.list[[k]])[ , 2:(length(x) + 1)]
    AA[ , k + 1] <- (test %*% x)[ , 1]
  }
  AA <- data.table(AA)
  names(AA) <- c("time", paste("Trans ", 1:ntrans))
  attr(AA, "tmat") <- attr(object, "tmat")
  return(AA)
}

## Calculates cumulative transitions hazards for the full cohort, i.e no specific
## covariates
msfitAalenTot <- function(object, x = 1) {
  ntrans <- length(object)
  dt.list <-  vector(mode = "list", length = ntrans)
  n.jumps <- rep(0, ntrans)
  for(k in 1:ntrans) {
    dt.list[[k]] <- data.table(object[[k]]$cum)
    setkey(dt.list[[k]], "time")
  }
  A.list <- dt.list
  dt.time <- data.table(time = c(attr(object, "start.time"), attr(object,"tvec")), key = "time")
  for(k in 1:ntrans) {
    A.list[[k]] <- merge(dt.time, A.list[[k]], all.x = TRUE)
    A.list[[k]] <- na.locf(A.list[[k]], na.rm = FALSE)
    A.list[[k]][is.na(A.list[[k]])] <- 0
  }
  AA <- matrix(nrow = dim(dt.time)[1], ncol = ntrans + 1)
  AA[,1] <- dt.time$time
  for(k in 1:ntrans) {
    AA[ , k + 1] <- A.list[[k]]$`(Intercept)`
  }
  AA <- data.table(AA)
  names(AA) <- c("time", paste("Trans ", 1:ntrans))
  attr(AA, "tmat") <- attr(object, "tmat")
  return(AA)
}
### Calculates transition probabilities from the cumulative hazards
### object: as fitted with msfitAalen
probtransAalen <- function(object) {
  ntrans <- dim(object)[2]-1
  tmat <- attr(object, "tmat")
  dAmatrix <- apply(object[ , 2:(ntrans + 1)], 2, FUN = diff)
  S <- dim(tmat)[1]
  P <- diag(S)
  plist <- vector(mode = "list", length = S)
  for (s in 1:S) {
    plist[[s]] <- matrix(0, ncol = S, nrow = dim(dAmatrix)[1] + 1)
    plist[[s]][1, ] <- P[s, ]
  }
  
  for(k in 1:dim(dAmatrix)[1]) {
    dA <- rep(0, S*S)
    dA[match(1:ntrans, tmat)] <- c(dAmatrix[k, ])
    dA <- matrix(dA, ncol=S, nrow=S, byrow = F)
    diag(dA) <- -apply(dA, 1, FUN = sum)
    IplusdA <- diag(S) + dA
    P <- P %*% IplusdA
    for (s in 1:S) {
      plist[[s]][k + 1 ,] <- P[s, ]
    }
  }
  for (s in 1:S) {
    plist[[s]] <- data.table(plist[[s]])
    names(plist[[s]]) <- paste0("pstate", 1:S)
    plist[[s]] <- cbind(time = object$time, plist[[s]])
  }
  class(plist) <- "probtransAalen"
  return(plist)
}

## Plot all transition probabilities from a chosen state
plot.probtransAalen <- function(object, from = 1, ylim = c(0,1),
                                ylab = "Probability", xlab = "time",
                                legendxy = "topright") {
  transprob <- object[[from]]
  S <- dim(transprob)[2]
  tvec <- transprob$time
  plot(tvec, as.matrix(transprob)[,2], type = "l", ylim = ylim, ylab = ylab, xlab = xlab)
  if(S > 2) {
    for (s in 3:S) {
      lines(tvec, as.matrix(transprob)[,s], lty = s-1)
    }
  }
  leg <- paste0("P_", from, 1:(S-1), "(" ,min(tvec), ",t)")
  legend(legendxy, leg, lty = 1:(S-1), col = "black", cex = 1.5)
}
