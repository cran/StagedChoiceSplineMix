#' @title Generates a new set of knots for the following iteration
#' @description A sub-function of \code{\link{StagedChoiceSplineMix}}. This function generates a new set of knots for the following iteration. Please refer to Bruch et al. (in press) for the precise rule used.
#' @param num.knot See \code{\link{StagedChoiceSplineMix}} for details.
#' @param sp.knots See \code{\link{StagedChoiceSplineMix}} for details.
#' @param k See \code{\link{StagedChoiceSplineMix}} for details.
#' @references
#' Bruch, E., F. Feinberg, K. Lee (in press), "Detecting Cupid's Vector: Extracting Decision Rules from Online Dating Activity Data," \emph{Proceedings of the National Academy of Sciences}.
#' @seealso
#' \code{\link{StagedChoiceSplineMix}}
#' @export
move.knot<-function(num.knot,sp.knots,k)
{
  min.knot<-1
  max.knot<-num.knot
  km.prob<-c(0.2,0.6,0.2)

  newknot.b<-matrix(double(k*2),ncol=2)
  newknot.w<-matrix(double(k*2),ncol=2)

  for (i in 1:k){
    oldknot.b<-sp.knots[i,1:2]
    oldknot.w<-sp.knots[i,3:4]

    k1d<-sample(c(1:3), 1, replace = T, prob=km.prob)
    k1.new<-(k1d==1)*(oldknot.b[1]-1)+(k1d==2)*oldknot.b[1]+(k1d==3)*(oldknot.b[1]+1)

    k2d<-sample(c(1:3), 1, replace = T, prob=km.prob)
    k2.new<-(k2d==1)*(oldknot.b[2]-1)+(k2d==2)*oldknot.b[2]+(k2d==3)*(oldknot.b[2]+1)

    ind<-(k1.new==(min.knot-1)|k1.new==max.knot|k2.new<=k1.new|k2.new>max.knot)
    k1.new<-(ind)*oldknot.b[1]+(!ind)*k1.new
    k2.new<-(ind)*oldknot.b[2]+(!ind)*k2.new

    newknot.b[i,1]<-k1.new
    newknot.b[i,2]<-k2.new

    k1d<-sample(c(1:3), 1, replace = T, prob=km.prob)
    k1.new<-(k1d==1)*(oldknot.w[1]-1)+(k1d==2)*oldknot.w[1]+(k1d==3)*(oldknot.w[1]+1)

    k2d<-sample(c(1:3), 1, replace = T, prob=km.prob)
    k2.new<-(k2d==1)*(oldknot.w[2]-1)+(k2d==2)*oldknot.w[2]+(k2d==3)*(oldknot.w[2]+1)

    ind<-(k1.new==(min.knot-1)|k1.new==max.knot|k2.new<=k1.new|k2.new>max.knot)
    k1.new<-(ind)*oldknot.w[1]+(!ind)*k1.new
    k2.new<-(ind)*oldknot.w[2]+(!ind)*k2.new

    newknot.w[i,1]<-k1.new
    newknot.w[i,2]<-k2.new
  }
  newknot<-cbind(newknot.b,newknot.w)

  list(newknot=newknot)
}
