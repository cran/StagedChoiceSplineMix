#' @import plyr
#' @importFrom stats quantile
#' @importFrom methods is
#' @title Performs iterations between an EM algorithm for a mixture of two-stage logistic regressions with fixed candidate knots and knot movements
#' @description The function performs iterations between an EM algorithm for a mixture of two-stage logistic regressions with fixed candidate knots and knot movements. The function generates candidate knots for each splined variable. Three sub-functions (\code{\link{gen.init}},\code{\link{twostglogitregmixEM}},\code{\link{move.knot}}) are used within the function.
#' @param data
#' Raw data for \emph{StagedChoiceSplineMix} \cr \cr
#' Format
#' \itemize{
#' \item 1st column: id
#' \item 2nd column: the 1st stage binary variable (browsing)
#' \item 3rd column: the 2nd stage binary variable (writing) conditional on the 1st stage binary variable. It should be left blank (or NA) if 1st stage variable is equal to 0
#' \item The rest of columns: covariates including splined variables
#' }
#' See \bold{Format}(below) for details \cr
#' @param M
#' number of iterations (def: 100)
#' @param sp.cols
#' vector of column numbers of splined variables in a data set (if sp.col is 0, \code{\link{twostglogitregmixEM}} function should be used)
#' @param num.knots
#' vector of numbers of knot candidates for splined variables. (def: a vector all of whose entries are "19")
#' @param sp.knots
#' list of knot configurations. For each splined variable, a knot configuration is a k by 4 matrix whose rows represent latent classes and columns represent knots [browsing knot 1, browsing knot 2, writing knot 1, writing knot 2]. (def: approximately 1/3 and 2/3 of knot candidates for knot 1 and knot 2 respectively)
#' @param betab
#' matrix of starting points for betab (browsing parameters). If not given, \code{\link{gen.init}} generates starting points.
#' @param betaw
#' matrix of starting points for betaw (writing parameters). If not given, \code{\link{gen.init}} generates starting points.
#' @param lambda
#' vector of starting points for lambda (membership proportion). If not given, \code{\link{gen.init}} generates starting points
#' @param k
#' number of latent classes (def: 2)
#' @param nst
#' number of random multiple starting points to try given a knot configuration. For each knot configuration, the output with the largest log-likelihood is stored among \emph{nst} trials. (def: 20)
#' @param epsilon
#' stopping tolerance for the EM algorithm. (def:1e-06)
#' @param maxit
#' maximum number of  the EM iterations allowed. If convergence is not declared before maxit, the EM algorithm stops with an error message and generates new starting points. (def: 500)
#' @param maxrestarts
#' maximum number of restarts (due to a singularity problem) allowed in the EM iterations. If convergence is not declared before maxrestarts, the algorithm stops with an error message and generates new starting points. (def: 100)
#' @param maxer
#' maximum number of errors allowed within a given knot configuration. If convergence is not declared before maxer, it tries a new knot configuration. (def: 20)
#' @return
#' \emph{StagedChoiceSplineMix} returns a list of the following items: \cr \cr
#' \bold{best}: best output, i.e., that giving the largest log-likelihood among M outputs.
#' \itemize{
#' \item{best$loglik: log-likelihood}
#' \item{best$betab: parameter estimates of betab}
#' \item{best$betaw: parameter estimates of betaw}
#' \item{best$lambda: parameter estimates of lambda}
#' \item{best$sp.knots: knot configuration}
#' }
#' \bold{loglike}: vector of  log-likelihoods for M outputs.
#' @references
#' Bruch, E., F. Feinberg, K. Lee (in press), "Detecting Cupid's Vector: Extracting Decision Rules from Online Dating Activity Data," \emph{Proceedings of the National Academy of Sciences}.
#' @seealso
#' \code{\link{gen.init}}
#' \code{\link{move.knot}}
#' \code{\link{twostglogitregmixEM}}\cr
#' "mixtools" package version 1.0.3
#' @examples
#' #######################################
#' ###### 1. Generate data (simdata) #####
#' set.seed(77)
#' k<-3
#'
#' betab<-matrix(c(1.5,2.0,0.3,-0.2,-0.1,0.6,-2.5,0.5,1,3,-1,1,2,-0.5,-1.5),5,k,byrow=TRUE)
#' betaw<-matrix(c(-2.0,-1,0.3,-0.3,-0.2,0.2,-3,-2,0,3,2,-1.5,3,2,-1),5,k,byrow=TRUE)
#' sp1.knots<-matrix(c(8,14,4,11,5,15,5,12,5,13,7,14),3,4,byrow=TRUE)
#' lambda<-c(0.3,0.3,0.4)
#'
#' # Large data set
#' #n_id<-700
#' #nobsb<-200
#'
#' # Small data set
#' n_id<-100
#' nobsb<-100
#'
#' nb<-n_id*nobsb
#'
#' idb <- rep(1:n_id, each=nobsb)
#'
#' xb.1<-rbinom(nb,size=1,prob=0.7)
#' xb.com<-cbind(1,xb.1)
#' xb.sp1<-runif((nb),-2,2)
#' xb<-cbind(xb.com,xb.sp1)
#'
#' sp1.mat.b<- matrix(double(20*nb),ncol=20)
#' sp1.mat.b[,1]<-xb.sp1
#' sp1.quan<-quantile(xb.sp1,c(0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,
#' 0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95))
#'
#' for (i in 1:19){
#'    sp1.mat.b[,(i+1)]<-xb.sp1-sp1.quan[i]
#' }
#' sp1.mat.b[sp1.mat.b<0]<-0
#' sp1.mat.b[,1]<-xb.sp1
#'
#' xb.class<-list()
#' for (i in 1:k){
#'   xb.class[[i]]<-cbind(xb.com,sp1.mat.b[,1],sp1.mat.b[,(sp1.knots[i,1:2])+1])
#' }
#'
#' xbetab<-matrix(double(nb*k),ncol=k)
#' for (i in 1:k){
#'   xbetab[,i]<-data.matrix(xb.class[[i]])%*%betab[,i]
#' }
#'
#' num.idb<-as.vector(ddply(data.frame(idb),"idb",count)[,-1])
#'
#' w<-sample(c(1:k), n_id, replace = TRUE, prob = lambda)
#' wb<-rep(w,num.idb)
#' yb.temp<-matrix(double(nb*k),ncol=k)
#'
#' for (i in 1:k){
#'   yb.temp[,i]<-rbinom((nb), size = 1, prob = (1/(1+exp(-xbetab[, i]))))
#' }
#' yb<-sapply(1:nb, function(i) yb.temp[i,wb[i]])
#'
#' idw<-idb[yb==1]
#' nw<-length(idw)
#' sp1.mat.w<-sp1.mat.b[yb==1,]
#' xw.com<-xb.com[yb==1,]
#' xw<-xb[yb==1,]
#'
#' xw.class<-list()
#' for (i in 1:k){
#'   xw.class[[i]]<-cbind(xw.com,sp1.mat.w[,1],sp1.mat.w[,(sp1.knots[i,3:4])+1])
#' }
#'
#' xbetaw<-matrix(double(nw*k),ncol=k)
#' for (i in 1:k){
#'   xbetaw[,i]<-data.matrix(xw.class[[i]])%*%betaw[,i]
#' }
#'
#' num.idw<-as.vector(ddply(data.frame(idw),"idw",count)[,-1])
#' ww<-wb[yb==1]
#'
#' yw.temp<-matrix(double(nw*k),ncol=k)
#' for (i in 1:k){
#'   yw.temp[,i]<-rbinom((nw), size = 1, prob = (1/(1+exp(-xbetaw[, i]))))
#' }
#' yw<-sapply(1:nw, function(i) yw.temp[i,ww[i]])
#'
#' yb.aug<-cbind(1:length(yb),yb)
#' yw.aug<-cbind(which(yb==1),yw)
#'
#' colnames(yb.aug)[1]<-"num"
#' colnames(yw.aug)[1]<-"num"
#'
#' ybyw<-merge(yb.aug,yw.aug, all = TRUE)[,-1]
#'
#' simdata<-cbind(idb,ybyw,xb.1,sp1.mat.b[,1])
#' simdata[is.na(simdata)]<-""
#' simdata[,3]<-as.integer(simdata[,3])
#' colnames(simdata)[1]<-"userid"
#' colnames(simdata)[2]<-"browsed"
#' colnames(simdata)[3]<-"wrote"
#' colnames(simdata)[4]<-"x1"
#' colnames(simdata)[5]<-"sp1"
#'
#' ########################################
#' ##### 2. Run StagedChoiceSplineMix #####
#' ## number of latent classes
#' set.seed(66)
#' k<-3
#'
#' ## starting points: true parameters used in the data generation (optional)
#' betab<-matrix(c(1.5,2.0,0.3,-0.2,-0.1,0.6,-2.5,0.5,1,3,-1,1,2,-0.5,-1.5),5,k,byrow=TRUE)
#' betaw<-matrix(c(-2.0,-1,0.3,-0.3,-0.2,0.2,-3,-2,0,3,2,-1.5,3,2,-1),5,k,byrow=TRUE)
#' lambda<-c(0.3,0.3,0.4)
#'
#' ## number of random multiple starting points to try given a knot configuration
#' nst<-1
#'
#' ## vector of the columns of spline variables in the data set (required)
#' sp.cols<-5
#'
#' # vector of the numbers of candidate knots for splined variables (optional)
#' num.knots<-19
#'
#' ## true knot configuration used in the data generation
#' sp1.knots<-matrix(c(8,14,4,11,5,15,5,12,5,13,7,14),3,4,byrow=TRUE)
#'
#' ## list of knot configuration of spline variables (optional)
#' sp.knots<-list(sp1.knots)
#'
#' ## Run "StagedChoiceSplineMix"
#' out<-StagedChoiceSplineMix(data=simdata, M=1, sp.cols=sp.cols, num.knots, sp.knots,
#' betab, betaw, lambda, k=k, nst=nst, epsilon=1e-06, maxit=500, maxrestarts=100, maxer=20)
#'
#' ## output
#' out$loglike  # vector of  M log-likelihoods
#' out$best$loglik # log-likelihood of the best output
#' out$best$betab # betab estimates of the best output
#' out$best$betaw # betaw estimates of the best output
#' out$best$lambda # lambda estimates of the best output
#' out$best$sp.knots # knot configuration of the best output
#'
#' @format
#' The simulated data(simdata) in the examples is generated using information below: \cr \cr
#' \bold{1. User identifier: }\emph{userid}
#' \itemize{
#' \item {Number of users: 700}
#' \item {Number of observations per user: 200}
#' }
#' \bold{2. Dependent variables:}
#' \itemize{
#' \item {\emph{browsed}: 1st stage binary dependent variable}
#' \item {\emph{wrote}: 2nd stage binary dependent variable}
#' }
#' \bold{3. Covariates:}
#' \itemize{
#' \item {\emph{x1 }(discrete variable): random draws from a binomial distribution with 0.7 success probability
#' }
#' \item {\emph{sp1 }(splined variable): random draws from a uniform distribution between -2 and 2}
#' }
#' \bold{4. Number of latent classes:} 3 \cr \cr
#' \bold{5. True parameters:} \cr \cr
#' \bold{i.Betab(browsing)}
#' \tabular{cccc}{
#' \tab \bold{Class1} \tab \bold{Class2} \tab \bold{Class3} \cr
#' \bold{intercept} \tab 1.5 \tab 2 \tab 0.3 \cr
#' \bold{x1} \tab -0.2 \tab -0.1 \tab 0.6 \cr
#' \bold{sp1} \tab -2.5 \tab 0.5 \tab 1 \cr
#' \bold{sp1 knot1} \tab 3 \tab -1 \tab 1 \cr
#' \bold{sp1 knot2} \tab 2 \tab -0.5 \tab -1.5 \cr
#' }
#' \bold{ii.Betaw(writing)}
#' \tabular{cccc}{
#' \tab \bold{Class1} \tab \bold{Class2} \tab \bold{Class3} \cr
#' \bold{intercept} \tab -2 \tab -1 \tab 0.3  \cr
#' \bold{x1} \tab -0.3 \tab -0.2 \tab 0.2 \cr
#' \bold{sp1} \tab -3 \tab -2 \tab 0 \cr
#' \bold{sp1 knot1} \tab 3 \tab 2 \tab -1.5 \cr
#' \bold{sp1 knot2} \tab 3 \tab 2 \tab -1 \cr
#' }
#' \bold{iii.Lambda(membership proportion)}
#' \tabular{ccc}{
#' \bold{Class1} \tab \bold{Class2} \tab \bold{Class3} \cr
#' 0.3 \tab 0.3 \tab 0.4 \cr
#' }
#' \bold{6. True knot configuration:} \cr
#' 19 candidate knots for sp1 (20-iles)
#' \tabular{ccccc}{
#' \tab \bold{knot1(browsing)} \tab \bold{knot2(browsing)} \tab \bold{knot1(writing)} \tab \bold{knot2(writing)} \cr
#' \bold{Class1} \tab 8 \tab 14 \tab 4 \tab 11 \cr
#' \bold{Class2} \tab 5 \tab 15 \tab 5 \tab 12 \cr
#' \bold{Class3} \tab 5 \tab 13 \tab 7 \tab 14 \cr
#' }
#' @export
StagedChoiceSplineMix<-
  function (data=NULL, M=100, sp.cols=NULL, num.knots=NULL, sp.knots=NULL, betab=NULL, betaw=NULL, lambda=NULL, k=2, nst=20, epsilon=1e-06, maxit=500, maxrestarts=100, maxer=20){

    idb<- c(as.factor(data[,1]))
    yb<-data[,2]
    xb<-cbind(1,data[,-c(1,2,3)])
    yw<-data[,3][yb==1]
    xw<-xb[yb==1,]
    idw<-idb[yb==1]
    num.sp<-length(sp.cols)

    if (is.null(num.knots)) {
      num.knots<-rep(19,num.sp)
    }

    if (is.null(sp.knots)) {
      sp.knots<-NULL
      for (i in num.sp){
        sp.knots[[i]]<-matrix(rep(c(ceiling(num.knots[i]/3),2*ceiling(num.knots[i]/3)),k*2),nrow=k,byrow=T)
      }
    }
    sp.mat.b<-NULL
    sp.mat.w<-NULL

    for (i in 1:num.sp){
      sp.mat.b[[i]]<-matrix(double((num.knots[i]+1)*length(yb)),ncol=(num.knots[i]+1))
      sp.mat.w[[i]]<-matrix(double((num.knots[i]+1)*length(yw)),ncol=(num.knots[i]+1))

      sp.mat.b[[i]][,1]<-xb[,(sp.cols[i]-2)]
      sp.mat.w[[i]][,1]<-xw[,(sp.cols[i]-2)]

      dis<-floor(1000*(1/(num.knots[i]+1)))/1000
      sp.quan<-quantile(xb[,(sp.cols[i]-2)],seq(from=dis,to=(1-dis),by=dis))

      for (j in 1:num.knots[i]){
        sp.mat.b[[i]][,(j+1)]<-xb[,(sp.cols[i]-2)]-sp.quan[j]
        sp.mat.w[[i]][,(j+1)]<-xw[,(sp.cols[i]-2)]-sp.quan[j]
      }
      sp.mat.b[[i]][sp.mat.b[[i]]<0]<-0
      sp.mat.w[[i]][sp.mat.w[[i]]<0]<-0

      sp.mat.b[[i]][,1]<-xb[,(sp.cols[i]-2)]
      sp.mat.w[[i]][,1]<-xw[,(sp.cols[i]-2)]
    }

    x.com.b<-xb[,-(sp.cols-2)]
    x.com.w<-xw[,-(sp.cols-2)]

    xb.class<-list()
    xw.class<-list()

    for (i in 1:k){
      xb.nc<-NULL
      xw.nc<-NULL
      for (j in 1:num.sp){
        xb.nc<-cbind(xb.nc,sp.mat.b[[j]][,1],sp.mat.b[[j]][,(sp.knots[[j]][i,1:2])+1])
        xw.nc<-cbind(xw.nc,sp.mat.w[[j]][,1],sp.mat.w[[j]][,(sp.knots[[j]][i,3:4])+1])
      }
      xb.class[[i]]<-cbind(x.com.b,xb.nc)
      xw.class[[i]]<-cbind(x.com.w,xw.nc)
    }

    nb <- length(yb)
    pb <- ncol(xb.class[[1]])
    Nb <- rep(1, nb)
    nw <- length(yw)
    pw <- ncol(xw.class[[1]])
    Nw <- rep(1, nw)

    er<-0

    if (is.null(betab) |is.null(betaw) |is.null(lambda)) {
      tt<-1
      while (tt==1) {
        test <- try(tmpb <- gen.init(yb, xb.class, k, er) , silent=F)
        if(is(test,"try-error")) {tt<-1}
        else {tt<-0}
      }

      tt<-1
      while (tt==1) {
        test <- try(tmpw <- gen.init(yw, xw.class, k, er) , silent=F)
        if(is(test,"try-error")) {tt<-1}
        else {tt<-0}
      }

      betab <- tmpb$beta
      betaw <- tmpw$beta

      while (sum(is.na(betab))>0 || sum(is.na(betaw))>0) {

        tt<-1
        while (tt==1) {
          test <- try(tmpb <- gen.init(yb, xb.class, k, er), silent=F)
          if(is(test,"try-error")) {tt<-1}
          else {tt<-0}
        }

        tt<-1
        while (tt==1) {
          test <- try(tmpw <- gen.init(yw, xw.class, k, er) , silent=F)
          if(is(test,"try-error")) {tt<-1}
          else {tt<-0}
        }

        betab <- tmpb$beta
        betaw <- tmpw$beta
      }
      lambda <- tmpb$lambda
    }

    best<-list()
    loglike<-double(M)
    ll.old<--1e+300
    j<-0

    while (j < M) {
      j = j + 1
      cat("\n")
      cat("*************************************************************","\n")
      cat("*** Knot configuration: ",j,"\n")
      cat("*************************************************************","\n")
      ex<-1
      er<-0
      st<-0

      llst.old<--1e+300
      bestst<-list()

      while (ex==1) {
        cat("Number of errors within a knot configuration: ",er,"\n")

        out1 <- try(twostglogitregmixEM(yb=yb, idb=idb, yw=yw, idw=idw, xb=xb, xw=xw, xb.class=xb.class, xw.class=xw.class, sp.cols=sp.cols, num.knots=num.knots, sp.knots=sp.knots, betab=betab, betaw=betaw, lambda=lambda, k=k, epsilon=epsilon, maxit=maxit, maxrestarts=maxrestarts, verb = TRUE), silent=F)
        if(is(out1,"try-error")){
          ex<-1
          er<-er+1
        }

        if(er>=maxer) {
          ex<-0
          j <- j - 1

          for (i in 1:num.sp){
            sp.knots[[i]]<-move.knot(num.knots,sp.knots[[i]],k)$newknot
          }

          xb.class<-list()
          xw.class<-list()

          for (i in 1:k){
            xb.nc<-NULL
            xw.nc<-NULL
            for (j in 1:num.sp){
              xb.nc<-cbind(xb.nc,sp.mat.b[[j]][,1],sp.mat.b[[j]][,(sp.knots[[j]][i,1:2])+1])
              xw.nc<-cbind(xw.nc,sp.mat.w[[j]][,1],sp.mat.w[[j]][,(sp.knots[[j]][i,3:4])+1])
            }
            xb.class[[i]]<-cbind(x.com.b,xb.nc)
            xw.class[[i]]<-cbind(x.com.w,xw.nc)
          }
        }

        if (is(out1,"try-error")==F && er<maxer ) {
          st<-st+1
          cat("* Starting point number:",st,"\n","\n")

          if (llst.old<=out1$loglik){
            bestst<-out1
            llst.old<-out1$loglik
          }

          if (st==nst){
            ex<-0
            cat("Best log-likelihood among multiple starting points:",llst.old,"\n","\n")

            for (i in 1:num.sp){
              sp.knots[[i]]<-move.knot(num.knots,sp.knots[[i]],k)$newknot
            }

            xb.class<-list()
            xw.class<-list()

            for (i in 1:k){
              xb.nc<-NULL
              xw.nc<-NULL
              for (m in 1:num.sp){
                xb.nc<-cbind(xb.nc,sp.mat.b[[m]][,1],sp.mat.b[[m]][,(sp.knots[[m]][i,1:2])+1])
                xw.nc<-cbind(xw.nc,sp.mat.w[[m]][,1],sp.mat.w[[m]][,(sp.knots[[m]][i,3:4])+1])
              }
              xb.class[[i]]<-cbind(x.com.b,xb.nc)
              xw.class[[i]]<-cbind(x.com.w,xw.nc)

              loglike[j]<-bestst$loglik
            }
          }
        }

        tt<-1
        while (tt==1) {
          test <- try(tmpb <- gen.init(yb, xb.class, k, er) , silent=F)
          if(is(test,"try-error")) {tt<-1}
          else {tt<-0}
        }

        tt<-1
        while (tt==1) {
          test <- try(tmpw <- gen.init(yw, xw.class, k, er) , silent=F)
          if(is(test,"try-error")) {tt<-1}
          else {tt<-0}
        }

        betab <- tmpb$beta
        betaw <- tmpw$beta

        while (sum(is.na(betab))>0 || sum(is.na(betaw))>0) {

          tt<-1
          while (tt==1) {
            test <- try(tmpb <- gen.init(yb, xb.class, k, er), silent=F)
            if(is(test,"try-error")) {tt<-1}
            else {tt<-0}
          }

          tt<-1
          while (tt==1) {
            test <- try(tmpw <- gen.init(yw, xw.class, k, er) , silent=F)
            if(is(test,"try-error")) {tt<-1}
            else {tt<-0}
          }

          betab <- tmpb$beta
          betaw <- tmpw$beta
        }
        lambda <- tmpb$lambda
      }

      if (ll.old<=bestst$loglik){
        best<-bestst
        ll.old<-bestst$loglik
      }

      betab<-bestst$betab
      betaw<-bestst$betaw
      lambda<-bestst$lambda

    }
    list(best = best, loglike = loglike)
  }
