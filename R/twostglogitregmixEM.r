#' @importFrom stats glm.fit dbinom coef
#' @title Performs an EM algorithm for mixtures of two-stage logistic regressions
#' @description A sub-function of \code{\link{StagedChoiceSplineMix}}. This function performs an EM algorithm for mixtures of two-stage logistic regressions.
#' @param yb The 1st stage binary variable (browsing). See \code{\link{StagedChoiceSplineMix}} for details.
#' @param idb Corresponding id for yb. See \code{\link{StagedChoiceSplineMix}} for details.
#' @param yw The 2nd stage binary variable (writing). See \code{\link{StagedChoiceSplineMix}} for details.
#' @param idw Corresponding id for yw. See \code{\link{StagedChoiceSplineMix}} for details.
#' @param xb Corresponding covariance matrix for yb. See \code{\link{StagedChoiceSplineMix}} for details.
#' @param xw Corresponding covariance matrix for yw. See \code{\link{StagedChoiceSplineMix}} for details.
#' @param xb.class Corresponding latent classes for yb.
#' @param xw.class Corresponding latent classes for yw.
#' @param sp.cols See \code{\link{StagedChoiceSplineMix}} for details.
#' @param num.knots See \code{\link{StagedChoiceSplineMix}} for details.
#' @param sp.knots See \code{\link{StagedChoiceSplineMix}} for details.
#' @param betab See \code{\link{StagedChoiceSplineMix}} for details.
#' @param betaw See \code{\link{StagedChoiceSplineMix}} for details.
#' @param lambda See \code{\link{StagedChoiceSplineMix}} for details.
#' @param k See \code{\link{StagedChoiceSplineMix}} for details.
#' @param epsilon See \code{\link{StagedChoiceSplineMix}} for details.
#' @param maxit See \code{\link{StagedChoiceSplineMix}} for details.
#' @param maxrestarts See \code{\link{StagedChoiceSplineMix}} for details.
#' @param maxer See \code{\link{StagedChoiceSplineMix}} for details.
#' @param verb See \code{\link{StagedChoiceSplineMix}} for details.
#' @seealso
#' \code{\link{StagedChoiceSplineMix}} \cr
#' \code{\link[mixtools]{regmixEM}}
#' @export
twostglogitregmixEM<-
  function (yb=NULL, idb=NULL, yw=NULL, idw=NULL, xb=NULL, xw=NULL, xb.class=NULL, xw.class=NULL, sp.cols=NULL, num.knots=NULL, sp.knots=NULL, betab = NULL, betaw = NULL, lambda = NULL, k = NULL, epsilon = 1e-06, maxit = 500, maxrestarts=100, maxer=20, verb = FALSE)
  {
    if (is.null(sp.cols)){
      num.knot<-NULL
      sp.knots<-NULL

      for (i in 1:k){
        xb.class[[i]]<-xb
        xw.class[[i]]<-xw
      }
    }

    logit <- function(x) 1/(1 + exp(-x))
    nb <- length(yb)
    pb <- ncol(xb.class[[1]])
    Nb <- rep(1, nb)
    nw <- length(yw)
    pw <- ncol(xw.class[[1]])
    Nw <- rep(1, nw)

    xbetab<-matrix(double(nb*k),ncol=k)
    xbetaw<-matrix(double(nw*k),ncol=k)

    for (i in 1:k){
      xbetab[,i]<-data.matrix(xb.class[[i]])%*%betab[,i]
      xbetaw[,i]<-data.matrix(xw.class[[i]])%*%betaw[,i]
    }

    n_id <-max(idb)
    z <- matrix(nrow = n_id, ncol = k)
    zb <- matrix(0, nb, k)
    zw <- matrix(0, nw, k)

    diff <- 1
    iter <- 0

    lgtb<- dbinom(yb, size = Nb, prob = logit(xbetab))
    compb <- t(t(log(1e-323+lgtb)))

    lgtw<- dbinom(yw, size = Nw, prob = logit(xbetaw))
    compw <- t(t(log(1e-323+lgtw)))

    compbaug<-data.frame(cbind(idb,compb))
    compwaug<-data.frame(cbind(idw,compw))

    compbidsum.temp<-ddply(compbaug,"idb",numcolwise(sum))[,-1]
    compwidsum.temp<-ddply(compwaug,"idw",numcolwise(sum))[,-1]

    compbidsum<-data.matrix(compbidsum.temp)
    compwidsum<-data.matrix(compwidsum.temp)

    compidsum<-compbidsum+compwidsum
    maxcomp<-apply(compidsum,1,max)
    max.col<-apply(compidsum,1,which.max)

    comp.diff<-compidsum-maxcomp
    lambda.diff<-matrix(rep(lambda,n_id),nrow=n_id,byrow=T)/lambda[max.col]
    exp.diff<-lambda.diff*exp(compidsum-maxcomp)

    logcompidsum<-log(lambda[max.col])+maxcomp+log(apply(exp.diff,1,sum))
    obsloglik<-sum(logcompidsum)

    ll <- obsloglik

    restarts <- 0
    rs<-0
    uf<-0

    betab.old<-betab
    betaw.old<-betaw
    lambda.old<-lambda

    while (abs(diff) > epsilon && iter < maxit|diff<0 ){
      if (diff<0) {iter<-0}

      for (i in 1:k){
        xbetab[,i]<-data.matrix(xb.class[[i]])%*%betab[,i]
        xbetaw[,i]<-data.matrix(xw.class[[i]])%*%betaw[,i]
      }

      log.fyb.temp<-matrix(double(nb*k),nb,k)
      for (i in 1:nb){
        for (j in 1:k){
          log.fyb.temp[i,j]<-log(logit(xbetab[i,j])^yb[i]*(1-logit(xbetab[i,j]))^(1-yb[i]))
          if (is.infinite(log.fyb.temp[i,j])) {log.fyb.temp[i,j]<--1e+308}
        }
      }

      log.fyw.temp<-matrix(double(nw*k),nw,k)
      for (i in 1:nw){
        for (j in 1:k){
          log.fyw.temp[i,j]<-log(logit(xbetaw[i,j])^yw[i]*(1-logit(xbetaw[i,j]))^(1-yw[i]))
          if (is.infinite(log.fyw.temp[i,j])) {log.fyw.temp[i,j]<--1e+308}
        }
      }

      log.fyb<-data.frame(log.fyb.temp,idb)
      log.prod_fyb_id_temp<-ddply(log.fyb,"idb",numcolwise(sum))[,-1]
      log.prod_fyb_id<-data.matrix(log.prod_fyb_id_temp)

      log.fyw<-data.frame(log.fyw.temp,idw)
      log.prod_fyw_id_temp<-ddply(log.fyw,"idw",numcolwise(sum))[,-1]
      log.prod_fyw_id<-data.matrix(log.prod_fyw_id_temp)

      for (i in 1:n_id) {
        for (j in 1:k) {
          z.denom<-0
          for (h in 1:k) {
            z.denom<-z.denom+exp(log(lambda[h])+log.prod_fyb_id[i,h]+log.prod_fyw_id[i,h]-log(lambda[j])-log.prod_fyb_id[i,j]-log.prod_fyw_id[i,j])
          }
          z[i, j] = 1/(z.denom)
        }
      }

      for (i in 1:k) {
        zb[,i]<-rep(z[,i],rle(idb)$length)
        zw[,i]<-rep(z[,i],rle(idw)$length)
      }

      if (sum(is.na(zb))>0 | sum(is.na(zw))>0|sum(apply(z, 2, mean)==0)>=1|is.infinite(sum(apply(z, 2, mean)))) {
        uf<-1
        stop("Underflow!!!")
      }
      else {

        lambda <- apply(z, 2, mean)

        lm.outb = lapply(1:k, function(j) try(glm.fit(xb.class[[j]], cbind(yb,Nb - yb), weights = zb[, j], family = binomial()),silent = TRUE))
        betab = sapply(lm.outb, coef)

        lm.outw = lapply(1:k, function(j) try(glm.fit(xw.class[[j]], cbind(yw,Nw - yw), weights = zw[, j], family = binomial()),silent = TRUE))
        betaw = sapply(lm.outw, coef)

        for (i in 1:k){
          xbetab[,i]<-data.matrix(xb.class[[i]])%*%betab[,i]
          xbetaw[,i]<-data.matrix(xw.class[[i]])%*%betaw[,i]
        }

        lgtb<- dbinom(yb, size = Nb, prob = logit(xbetab))
        compb <- t(t(log(1e-323+lgtb)))

        lgtw<- dbinom(yw, size = Nw, prob = logit(xbetaw))
        compw <- t(t(log(1e-323+lgtw)))

        compbaug<-data.frame(cbind(idb,compb))
        compwaug<-data.frame(cbind(idw,compw))

        compbidsum.temp<-ddply(compbaug,"idb",numcolwise(sum))[,-1]
        compwidsum.temp<-ddply(compwaug,"idw",numcolwise(sum))[,-1]

        compbidsum<-data.matrix(compbidsum.temp)
        compwidsum<-data.matrix(compwidsum.temp)

        compidsum<-compbidsum+compwidsum
        maxcomp<-apply(compidsum,1,max)
        max.col<-apply(compidsum,1,which.max)

        comp.diff<-compidsum-maxcomp
        lambda.diff<-matrix(rep(lambda,n_id),nrow=n_id,byrow=T)/lambda[max.col]
        exp.diff<-lambda.diff*exp(compidsum-maxcomp)

        logcompidsum<-log(lambda[max.col])+maxcomp+log(apply(exp.diff,1,sum))
        newobsloglik<-sum(logcompidsum)

        if (abs(newobsloglik) == Inf || is.na(newobsloglik)){
          cat("Singularity!!!","\n")

          restarts <- restarts + 1
          cat("restarts:",restarts,"\n")
          if (restarts > maxrestarts)
          {
            rs<-1
          }
          if (restarts > maxrestarts)
            stop("Too many tries!!!")

          betab <- betab.old+matrix(rnorm(dim(betab)[1]*dim(betab)[2],0,0.1*(1+restarts)),dim(betab)[1],dim(betab)[2])
          betaw <- betaw.old+matrix(rnorm(dim(betaw)[1]*dim(betaw)[2],0,0.1*(1+restarts)),dim(betaw)[1],dim(betaw)[2])
          lambda <-lambda.old

          lambda[which.max(lambda)]<-lambda[which.max(lambda)]-0.05
          lambda[which.min(lambda)]<-lambda[which.min(lambda)]+0.05

          xbetab<-matrix(double(length(yb)*k),ncol=k)
          xbetaw<-matrix(double(length(yw)*k),ncol=k)


          for (i in 1:k){
            xbetab[,i]<-data.matrix(xb.class[[i]])%*%betab[,i]
            xbetaw[,i]<-data.matrix(xw.class[[i]])%*%betaw[,i]
          }

          diff <- 1
          iter <- 0

          lgtb<- dbinom(yb, size = Nb, prob = logit(xbetab))
          compb <- t(t(log(1e-323+lgtb)))

          lgtw<- dbinom(yw, size = Nw, prob = logit(xbetaw))
          compw <- t(t(log(1e-323+lgtw)))

          compbaug<-data.frame(cbind(idb,compb))
          compwaug<-data.frame(cbind(idw,compw))

          compbidsum.temp<-ddply(compbaug,"idb",numcolwise(sum))[,-1]
          compwidsum.temp<-ddply(compwaug,"idw",numcolwise(sum))[,-1]

          compbidsum<-data.matrix(compbidsum.temp)
          compwidsum<-data.matrix(compwidsum.temp)

          compidsum<-compbidsum+compwidsum
          maxcomp<-apply(compidsum,1,max)
          max.col<-apply(compidsum,1,which.max)

          comp.diff<-compidsum-maxcomp
          lambda.diff<-matrix(rep(lambda,n_id),nrow=n_id,byrow=T)/lambda[max.col]
          exp.diff<-lambda.diff*exp(compidsum-maxcomp)

          logcompidsum<-log(lambda[max.col])+maxcomp+log(apply(exp.diff,1,sum))
          obsloglik<-sum(logcompidsum)

          ll <- obsloglik
        }
        else {
          diff <- newobsloglik - obsloglik
          obsloglik <- newobsloglik
          ll <- c(ll, obsloglik)
          iter <- iter + 1
          if (verb) {
            cat("EM iteration=", iter,",", "diff=", diff, ",","log-likelihood=", obsloglik, "\n")
          }
        }
      }
    }
    if (iter == maxit) {
      cat("WARNING! NOT CONVERGENT!", "\n")
    }
    betab <- matrix(betab, ncol = k)
    betaw <- matrix(betaw, ncol = k)

    if (is.null(sp.cols)){
      col.names<-c("intercept")
      for (i in 1:(dim(xb)[2]-1)){
        col.names<-c(col.names,paste("x",".",i,sep =""))
      }
    }
    else {
      col.names<-c("intercept")
      for (i in 1:(dim(xb)[2]-length(sp.cols)-1)){
        col.names<-c(col.names,paste("x",".",i,sep =""))
      }
      for (i in 1:(length(sp.cols))){
        col.names<-c(col.names,paste("sp",".",i,sep =""))
        col.names<-c(col.names,paste("sp",".",i," ","knot1",sep =""))
        col.names<-c(col.names,paste("sp",".",i," ","knot2",sep =""))
      }
    }
    rownames(betab) <- col.names
    rownames(betaw) <- col.names

    colnames(betab) <- c(paste("class", ".", 1:k, sep = ""))
    colnames(betaw) <- c(paste("class", ".", 1:k, sep = ""))

    names(lambda) <- c(paste("class", ".", 1:k, sep = ""))

    names(obsloglik) <-c("log-likelihood")

    colnames(zb) <- c(paste("class", ".", 1:k, sep = ""))
    colnames(zw) <- c(paste("class", ".", 1:k, sep = ""))
    colnames(z) <- c(paste("class", ".", 1:k, sep = ""))

    if (is.null(sp.cols)==F){
      for (i in 1:length(sp.cols)){
        colnames(sp.knots[[i]]) <- c("browsing knot1","browsing knot2", "writing knot1", "writing knot2")
        rownames(sp.knots[[i]]) <- c(paste("class", ".", 1:k, sep = ""))
      }
    }

    list(xb = xb, xw = xw, xb.class=xb.class, xw.class=xw.class
         ,k=k, yb = yb, yw = yw, idb=idb, idw=idw, lambda = lambda, betab = betab, betaw = betaw, loglik = obsloglik,
         posteriorzb = zb, posteriorzw = zw, posteriorz = z, all.loglik = ll,
         sp.cols=sp.cols, num.knots=num.knots, sp.knots=sp.knots)
  }
