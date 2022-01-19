#' MCMC.Summary
#'
#' final quantity presentation (output) MC GUM supplement 2
#'
#' @param output [matrix] (**required**) output quantity
#' @param p [integer] (**required**) coverage probability
#'
#' @return alist object with the following elements
#' @return $output.sort
#' @return $esp.sort
#' @return $dev.sort
#' @return $cov.map
#' @return $map.p
#' @return $low.cov high coverage value
#' @return $high.cov low coverage value
#'
#'
#' @export

MCMC.Summary <-
function (output,p){

output<-output[-seq(1,10),] #the initial values are eliminated, in case there is no burn.in.
m=ncol(output) #number of quantities
n.iter =length(output)/m #number of iteration

#average matrix
output.av =colMeans(output)
output.av.m<-t(as.matrix(output.av)) #other representation of average matrix

#coverage matrix
output.cov<-matrix(rep(0,m*m),m,m)
for (i in 1:n.iter){
output.cov<-output.cov+(output[i,]-output.av)%*%t(output[i,]-output.av)
}
output.cov<-output.cov/(n.iter-1)

#construct a rectangular region
L<-rbind(diag(output.cov))^0.5 #diagonal inverse
y.ro<-L[rep(1,n.iter),]^-1*(output-output.av.m[rep(1,n.iter),]) #transformed y points
output.sort<-0
for (i in 1:n.iter){
output.sort<-c(output.sort,max(y.ro[i,]))
}
output.sort<-output.sort[-1]
output.sort<-sort(output.sort,decreasing = FALSE)#rank the transformed points
q<-(1+p)/2
k.q<-output.sort[q*length(output.sort)]
Low<-output.av-k.q*L
High<-output.av+k.q*L

#subdivide this rectangular region into a mesh of small rectangle
mesh<-n.iter^(1/3)
Step<-(High-Low)/mesh

#assign each output quantity value to the small rectangle containing it
output<-cbind(output,floor((output-output.av.m[rep(1,n.iter),])/Step[rep(1,n.iter),]))
output<-cbind(output,abs(rowSums(output[,c(3,4)])))

x.range<-max(output[,3])-min(output[,3])+1
y.range<-max(output[,4])-min(output[,4])+1

output.map<-array(rep(0,mesh^2),dim=c(x.range,y.range))

xmin<-min(output[,3])-1
ymin<-min(output[,4])-1
for (i in 1:n.iter){
xi<-output[i,3]-xmin
yi<-output[i,4]-ymin
output.map[xi,yi]<-output.map[xi,yi]+1
}

#use the fraction of output quantity value as the approximate probability

#sort the rectangle in terms of decreasing probability
ii<-order(output.map[,],decreasing=TRUE)
output.ecdf<-output.map[ii]

# form the cumulative sum of probabilities for these listed rectangle, stopping when the sum is not smaller than p
n.cum<-1
M<-sum(output.ecdf)
cum<-0
while(cum<p*M){
cum<-sum(output.ecdf[seq(1,n.cum)])
n.cum<-n.cum+1
}


#take the corresponding set of rectangles as defining the smallest coverage region
output.map.rank<-array(rank(output.map,ties.method="min"),dim=c(x.range,y.range))
opp_n.cum<-x.range*y.range-n.cum
output.map.rank[output.map.rank<opp_n.cum]<-0
output.map.rank[output.map.rank>=opp_n.cum]<-1

Taxs<-seq(min(output[,2]),max(output[,2]),length.out=y.range)
Xaxs<-seq(min(output[,1]),max(output[,1]),length.out=x.range)

mfcol=c(2, 2)
#contour(x=Xaxs,y=Taxs,output.map)
#plot(output[,1:2])

filled.contour(x=Xaxs,y=Taxs,output.map,plot.axes={ axis(1); axis(2); points(output[,1:2])})
#filled.contour(x=Xaxs,y=Taxs,output.map.rank,plot.axes={axis(1); axis(2); lines(c(High[,1],High[,1],Low[,1],Low[,1],High[,1]),c(High[,2],Low[,2],Low[,2],High[,2],High[,2]))})

return(list(output.sort=output,esp.sort=output.av,dev.sort=output.cov,map=output.map,map.p=output.map.rank,low.cov=Low,high.cov=High))
}
