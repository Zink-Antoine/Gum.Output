#'Output.Summary
#'
#' final quantity presentation (output) MC GUM suppl1
#'
#' @param output output quantity
#' @param p coverage probability
#'
#'
#' @return a list object with the following elements
#' @return $output.sort
#' @return $esp.sort
#' @return dev.sort
#' @return cov.star
#' @return low.cov
#' @return high.cov
#'
#'
#' @export
#'
Output.Summary <-
function (output,p){
output.sort<-sort(output,decreasing = FALSE)
# absolute growth test
m<-length(output.sort)
	while(length(m)>0){
	output.diff<-output.sort-data.table::shift(output.sort,n=1)
	m<-which(output.diff==0) #index where the deviation from the previous number is zero
	output.sort[m]<-output.sort[m]+1e-128 #includes minimal disturbance
}

#coverage calculate
M<-length(output)
q<-p*M

#symm?trique
#r<-(M-q)/2
#low<-r
#high<-r+q
#cov<-output.sort[high]-output.sort[low]

#shortest
cov.star<-M
for (r in 1:(M-q)){
low<-r
high<-r+q
cov<-output.sort[high]-output.sort[low]
if(cov<cov.star) {
cov.star<-cov
low.star<-low
high.star<-high
}

}

return(list(output.sort=output.sort,esp.sort=mean(output.sort),dev.sort=sd(output.sort),cov.star=cov.star,low.cov=output.sort[low.star],high.cov=output.sort[high.star]))
}

