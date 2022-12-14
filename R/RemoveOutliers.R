#' RemoveOutliers
#'
#'
#' procedure to remove outliers from matrix of MC simulation using the z-score
#'
#'
#' @param tmp [matrix] (**required**) matrix of MCMC
#'
#' @details the function use the z-score.
#'
#' @return a tmp matrix modified after discard of the outliers
#'
#' @seealso outliers package especially grubbs.test (type=10) give the same result as z-score
#'
#' @export
#'
RemoveOutliers<-
  function(tmp){
      mu<-apply(tmp,2,mean)
      s<-apply(tmp,2,sd)
      z<-t((t(tmp)-mu)/s)

    if (max(z)>3) {
      tmp<-tmp[apply(z<3,1,all),]
      Recall(tmp=tmp)
    }
    else{
      coda::as.mcmc(tmp)
    }
  }
