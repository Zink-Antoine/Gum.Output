#' Adaptive MC procedure
#'
#' adapted from GUM suppl.2
#' 7.8 Adaptive Monte Carlo procedure
#'
#' @param nOutput [integer] (**with default value**) number of model output(default nOutput=2)
#' @param ndig [integer] (**with default value**) number of significant decimal digits (default ndig=2)
#' @param p [integer] (**with default value**) coverage probability (default p= 0.95)
#' @param outliers [logical](**with default**) logical parameters indicating whether or not outliers are discarded. TRUE, the outliers are keeping, FALSE, the outliers are rejected ( default outliers=TRUE) 
#' @param FUN [function] (**required**) a simulation MC model
#' @param n.iter [integer] (**with default**) number of iteration for each simulation (default n.iter= 1e4)
#' @param ... further parameters passed to the FUN function
#'
#'
#' @return a list including
#' @return ndig [integer] number of significant decimal digits
#' @return M [integer] the number of draws M= h x n.iter (h number of simulation and n.iter number of iteration for each simulation)
#' @return y [vector] an estimate of measurand Y
#' @return uy [vector] the standard uncertainties associated with the estimates,
#' @return Ry [matrix] correlation coefficients rij = r(yi, yj) associated with pairs of the estimates,
#' @return kp [integer] a coverage factor defining a 100p % coverage region for Y,
#' 
#' @return as stored file
#' @return trial_h [mcmc] list of the simulations
#'
#' @export
#'
#' @references JCGM-WG1 (2011) data – Supplement 2 to the “Guide to the expression of uncertainty in measurement” – Extension to any number of output quantities. Guide JCGM 102:2011. Sèvres: BIPM, IEC, IFCC, ILAC, ISO, IUPAC, IUPAP and OIML.
#'
#' @examples
#'
#' 
#' ## example 9.2.2 GUM suppl.2 case1
#' Model<-
#' function(n.iter){
#' X1<-rnorm(n.iter,0,1)
#' X2<-rnorm(n.iter,0,1)
#' X3<-rnorm(n.iter,0,1)
#' Y1<-X1+X3
#' Y2<-X2+X3
#' cbind(Y1,Y2)
#' }
#'
#' AdaptivMC(nOutput=2,ndig=3,p=0.95,FUN=Model,n.iter=1000)
#' 
#' 
#' \dontrun{
#' AR(1) example
#' AR1<-
#' function(n.iter,rho=1,sig=1){
#'  X<-Xn<-rnorm(1,0,sig)
#'  Y<-0
#'  for (i in 2:n.iter){
#'    Yn<-rnorm(1,0,sig)
#'    Xn<-rho*Xn+Yn
#'    X<-rbind(X,Xn)
#'    Y<-rbind(Y,Yn)
#'  }
#'  out<-cbind(X,Y)
#' }
#' 
#' AdaptivMC(nOutput=2,ndig=2,p=0.95,FUN=AR1,n.iter=1000,rho=0.99)
#'
#' }
#'
#'
#' \dontrun{
#' require(TLpack)
#'
#' data(TLetru)
#' table<-Lum(TLetru,Doseb0=360,Dosea=0,alpha=FALSE,supra=FALSE)
#' B<-table$b
#' N<-table$n
#' table.norm<-B
#'  for (j in 1:3){
#'    for (i in 1:3)
#'      {
#'        table.norm[,i,j]<-(B[,i,j])/sum(N[seq(350,400),1,(j-1)*3+i])
#'        }
#'    }
#'  table.data<-cbind(table.norm[,,1],table.norm[,,2],table.norm[,,3])
#'  ii<-c(1,4,7,2,5,8,3,6,9)
#'  table.data<-table.data[seq(1,475),ii]
#'  Dose<-as.numeric(colnames(table.data))
#'  df.T<-matrix(rep(seq(26,500),9),475,9)
#'  df.y<-table.data[,1:9]
#' AdaptivMC(nOutput=5,ndig=2,outliers=FALSE,FUN=Slice5,Dose=Dose,df.T=df.T,df.y=df.y,n.thin=1)
#' }
#'


AdaptivMC<-
  function(nOutput=2,ndig=2,p=0.95,outliers=TRUE,FUN,n.iter=10^4,...){
    #7.8.3 Adaptive procedure
    # a) set ndig to an appropriate small positive integer (see 7.8.2);
    ndig<-ndig
    #b) set M = max(J, 104), where J is the smallest integer greater than or equal to 100/(1 − p);
    J<-ceiling(100/(1-p))
    M<-max(J,n.iter)
    #c) set h = 1, denoting the first application of MCM in the sequence;
    Y_h<-matrix(nrow=10,ncol=nOutput)
    U_h<-matrix(nrow=10,ncol=nOutput)
    R_h<-vector(mode="list",length=10)
    EigenMax_h<-array(dim=10)
    kp_h<-array(dim=10)
    trial_h<-vector(mode="list",length=10)

    h<-0
    #d) carry out M Monte Carlo trials, as in 7.3 and 7.4;
    #e) use the M vector output quantity values y1,..,yM so obtained to calculate y(h), u(y(h)), Ry(h) and k(h) p as an estimate of Y , the associated standard uncertainties, the associated correlation matrix and a coverage factor for a 100p % coverage region, respectively, i.e. for the hth member of the sequence;
    # f) if h ≤ 10, increase h by one and return to step d);

    AdaptivMC_Boucle(h,Y_h,U_h,R_h,EigenMax_h,kp_h,trial_h,ndig,nOutput,p,outliers,FUN=FUN,n.iter=M,...)
    }


AdaptivMC_Boucle<-
  function(h,Y_h,U_h,R_h,EigenMax_h,kp_h,trial_h,ndig,nOutput,p,outliers,FUN,n.iter,...){

    trial<-FUN(n.iter=n.iter,...) #MC model
    if (!outliers) trial<-RemoveOutliers(trial)

    h<-h+1
    #print(h)

    Y_h[h,]<-colMeans(trial)
    U_h[h,]<-apply(trial,2,sd)
    R_h[[h]]<-cor(trial)
    EigenMax_h[h]<-max(eigen(R_h[[h]],only.values = TRUE)$value)
    kp_h[h]<-CovFac(cov(trial),trial,p)
    trial_h[[h]]<-trial

    if (h<10) Recall(h=h,Y_h=Y_h,U_h=U_h,R_h=R_h,EigenMax_h=EigenMax_h,kp_h=kp_h,trial_h=trial_h,ndig=ndig,nOutput=nOutput,p=p,outliers=outliers,FUN=FUN,n.iter=n.iter,...=...)
    else
    {
      #g) for j = 1,...m, calculate the standard deviation syj associated with the average of the estimates yj,
      y<-colMeans(Y_h)
     # sy<-apply(Y_h,2,sd)/sqrt(h)
      sy<-sqrt(sapply(apply(Y_h,2,initseq),function(x){x$var.con})/h)

      #h) calculate the counterpart of this statistic for the components of u(y(h)) and for λmax and k(h) p;
      uy<-colMeans(U_h)
      suy<-apply(U_h,2,sd)/sqrt(h)

      EigenMax<-mean(EigenMax_h)
      sEigenMax<-sd(EigenMax_h)/sqrt(h)

      kp<-mean(kp_h)
      skp<-sd(kp_h)/sqrt(h)

      #i) use all h × M model values available so far to form values for u(y), Ry and kp;

      #for j = 1, . . m, calculate the numerical tolerances δj associated with u(yj) as in 7.8.2.1 and 7.8.2.2;
      delta<-lapply(uy,Tolerance,ndig)

      # k) calculate the numerical tolerance ρ associated with the matrix Ry of correlation coefficients
      rho<-Tolerance(EigenMax,ndig)

      # l) calculate the numerical tolerance κp associated with kp
      kap<-Tolerance(kp,ndig)

      print(c(2*sy>delta,2*suy>delta,2*sEigenMax>rho,2*skp>kap))

      #if for any j = 1, . . . , m, 2syj or 2su(yj) exceeds δj, or 2sλmax exceeds ρ, or 2skp exceeds κp, increase h by one and return to step d)
      if(any(2*sy>delta,2*suy>delta,2*sEigenMax>rho,2*skp>kap)){
        Y_h.new<-matrix(nrow=h+1,ncol=nOutput)
        U_h.new<-matrix(nrow=h+1,ncol=nOutput)

        Y_h.new[seq(1,h),]<-Y_h
        U_h.new[seq(1,h),]<-U_h

        Y_h<-Y_h.new
        U_h<-U_h.new

        save(trial_h,file=sprintf("chain%03d.rda",h))
        
        Recall(h=h,Y_h=Y_h,U_h=U_h,R_h=R_h,EigenMax_h=EigenMax_h,kp_h=kp_h,trial_h=trial_h,ndig=ndig,nOutput=nOutput,p=p,outliers=outliers,FUN=FUN,n.iter=n.iter,...=...)

        }
      else{
        Ry<-array(dim=c(nOutput,nOutput))
        Ry<-Reduce("+",R_h)/length(R_h)
        
        save(trial_h,file=sprintf("chain%03d.rda",h))

        return(list(ndig=ndig,M=h*n.iter,y=y,uy=uy,Ry=Ry,kp=kp))
      }
    }
  }





Tolerance<-
  function(z,ndig=2){
    # 7.8.2 Numerical tolerance associated with a numerical value
    l<-ceiling(log10(abs(z)))-ndig
    c<-trunc(ceiling(z/10^l)) #useful ?
    delta<-0.5*10^l
    return(delta)
  }

CovFac<-
  function(Uy,yr,p){
    #7.7.2 Hyper-ellipsoidal coverage region
    # a) Transform the points yr, denoting the transformed points by  ̊yr, according to  ̊yr = L−1(yr − y), r = 1, . . . , M, (21)
    # where L is given by the Cholesky decomposition of Uy;
    L<-chol(Uy)
    id<- diag(nrow=dim(L)[1],ncol=dim(L)[2])
    Linv<-solve(L,id)
    y0r<-t(t(yr)-colMeans(yr))%*%Linv
    n<-dim(y0r)[1]

    # b) Sort the transformed points̊ yr according to increasing value of dr, where d2 r =  ̊y> r  ̊yr = m ∑ j=1  ̊y2 j,r, r = 1, . . . , M ;
    d2r<-array(dim=n)
    for(i in 1:n){
      d2r[i]<-crossprod(y0r[i,])
    }
    ii<-order(d2r)
    y0r_sort<-cbind(d2r,y0r)[ii,]

    # c) Use the sorted  ̊yr to determine the coverage factor kp such that a fraction p of the  ̊yr satisfies dr < kp;
    i<-n*p
    kp<-sqrt(y0r_sort[i+1,1])[[1]]

    #d) Take the hyper-ellipsoid defined by equation (20) as the boundary of a 100p % coverage region for Y .

  }

####### modified from Geyer (2011)

LoadChain<-
  function(){
    files<-list.files(pattern="chain\\d+.rda")
    kcheck<-length(files)
    load(file=files[kcheck])
  }

