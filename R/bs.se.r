#' @importFrom stats rbinom sd
#' @title Performs a parametric bootstrapping standard error approximation
#' @description The function performs a parametric bootstrapping standard error approximation for mixtures of two-stage logistic regressions.
#' @param output
#' output of \code{\link{StagedChoiceSplineMix}}
#' @param B
#' number of bootstrap samples (def:100)
#' @seealso
#' \code{\link{StagedChoiceSplineMix}} \cr
#' \code{\link[mixtools]{boot.se}}
#' @examples
#' ## parametric bootstrapping to calculate standard error:
#' ## use best output from StagedChoiceSplineMix exampl
#' #out.se<-bs.se(output=out$best, B=100)
#' #out.se$betab.se
#' #out.se$betaw.se
#' #out.se$lambda.se
#' @export
bs.se<- function (output, B = 100)
{
  yb <- output$yb
  yw <- output$yw
  xb.class<-output$xb.class
  xw.class<-output$xw.class
  xb<-output$xb
  xw<-output$xw
  k<-output$k

  betab<-output$betab
  betaw<-output$betaw
  lambda<-output$lambda
  idb<-output$idb
  idw<-output$idw
  sp.cols<-output$sp.cols
  num.knots<-output$num.knots
  sp.knots<-output$sp.knots

  n_id<-max(idb)
  nb<-length(yb)
  nw<-length(yw)

  xbetab<-matrix(double(length(yb)*k),ncol=k)
  xbetaw<-matrix(double(length(yw)*k),ncol=k)

  for (i in 1:k){
    xbetab[,i]<-data.matrix(xb.class[[i]])%*%betab[,i]
    xbetaw[,i]<-data.matrix(xw.class[[i]])%*%betaw[,i]
  }

  j <- 0
  betab.bs <- NULL
  betaw.bs <- NULL
  lambda.bs <- NULL

  num.idb<-as.vector(ddply(data.frame(idb),"idb",count)[,-1])
  num.idw<-as.vector(ddply(data.frame(idw),"idw",count)[,-1])

  while (j < B) {
    j <- j + 1

    cat("\n")
    cat("Bootstrap sample ",j,"\n")
    w<-sample(c(1:k), n_id, replace = T, prob = lambda)

    wb<-rep(w,num.idb)
    ww<-rep(w,num.idw)

    yb.sim.temp<-matrix(double(nb*k),ncol=k)
    yw.sim.temp<-matrix(double(nw*k),ncol=k)

    for (i in 1:k){
      yb.sim.temp[,i]<-rbinom((nb), size = 1, prob = (1/(1+exp(-xbetab[, i]))))
      yw.sim.temp[,i]<-rbinom((nw), size = 1, prob = (1/(1+exp(-xbetaw[, i]))))
    }

    yb.sim<-sapply(1:nb, function(i) yb.sim.temp[i,wb[i]])
    yw.sim<-sapply(1:nw, function(i) yw.sim.temp[i,ww[i]])

    bs.run <- try(twostglogitregmixEM(yb.sim, idb, yw.sim, idw, xb, xw, xb.class, xw.class, sp.cols, num.knots, sp.knots, betab, betaw, lambda, k, epsilon = 1e-06, maxit = 500, maxrestarts=100, maxer=20, verb = T), silent = T)


    if (class(bs.run) == "try-error"){
      j = j - 1
    }
    else {
      lambda.bs <-cbind(lambda.bs, bs.run$lambda)
      betab.bs <-cbind(betab.bs, as.vector(bs.run$betab))
      betaw.bs <- cbind(betaw.bs, as.vector(bs.run$betaw))
    }
  }

  lambda.se = apply(lambda.bs, 1, sd)
  betab.se = matrix(apply(betab.bs, 1, sd), ncol = k)
  betaw.se = matrix(apply(betaw.bs, 1, sd), ncol = k)

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
  rownames(betab.se) <- col.names
  rownames(betaw.se) <- col.names

  colnames(betab.se) <- c(paste("class", ".", 1:k, sep = ""))
  colnames(betaw.se) <- c(paste("class", ".", 1:k, sep = ""))

  names(lambda.se) <- c(paste("class", ".", 1:k, sep = ""))

  bs.list = list(lambda.bs = lambda.bs, lambda.se = lambda.se
                 ,betab.bs = betab.bs, betab.se = betab.se, betaw.bs = betaw.bs, betaw.se = betaw.se)
}

