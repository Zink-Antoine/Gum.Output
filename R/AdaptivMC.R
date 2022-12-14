#' Adaptive MC procedure
#'
#' adapted from GUM suppl.2
#' 7.8 Adaptive Monte Carlo procedure
#'
#' @param nOutput [integer] (**with default value**) number of model output(default nOutput=2)
#' @param ndig [integer] (**with default value**) number of significant decimal digits (default ndig=2)
#' @param p [integer] (**with default value**) coverage probability (default p= 0.95)
#' @param FUN [function] (**required**) a simulation MC model
#' @param n.iter [integer] (**with default**) number of iteration for each simulation (default n.iter= 1e4)
#' @param ... further parameters passed to the FUN function
#'
#'
#' @return a list including
#' @return y [vector] an estimate of measurand Y
#' @return uy [vector] the standard uncertainties associated with the estimates,
#' @return Ry [matrix] correlation coefficients rij = r(yi, yj) associated with pairs of the estimates,
#' @return kp [integer] a coverage factor defining a 100p % coverage region for Y,
#'
#' @export
#'
#' @references JCGM-WG1 (2011) data – Supplement 2 to the “Guide to the expression of uncertainty in measurement” – Extension to any number of output quantities. Guide JCGM 102:2011. Sèvres: BIPM, IEC, IFCC, ILAC, ISO, IUPAC, IUPAP and OIML.
#'
#' @examples
#' \dontrun{
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
#' AdaptivMC(nOutput=5,ndig=2,FUN=Slice5,Dose=Dose,df.T=df.T,df.y=df.y,n.thin=1)
#' }
#'


AdaptivMC<-
  function(nOutput=2,ndig=2,p=0.95,FUN,n.iter=10^4,...){
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

    h<-1
    #d) carry out M Monte Carlo trials, as in 7.3 and 7.4;
    #e) use the M vector output quantity values y1,..,yM so obtained to calculate y(h), u(y(h)), Ry(h) and k(h) p as an estimate of Y , the associated standard uncertainties, the associated correlation matrix and a coverage factor for a 100p % coverage region, respectively, i.e. for the hth member of the sequence;
    # f) if h ≤ 10, increase h by one and return to step d);

    AdaptivMC_Boucle(M,h,Y_h,U_h,R_h,EigenMax_h,kp_h,ndig,nOutput,p,FUN,n.iter=M,...)
    }


AdaptivMC_Boucle<-
  function(M,h,Y_h,U_h,R_h,EigenMax_h,kp_h,ndig,nOutput,p,FUN,...){

    trial<-FUN(...) #MC model
    trial<-RemoveOutliers(trial)

    Y_h[h,]<-colMeans(trial)
    U_h[h,]<-apply(trial,2,sd)
    R_h[[h]]<-cor(trial)
    EigenMax_h[h]<-max(eigen(R_h[[h]],only.values = TRUE)$value)
    kp_h[h]<-CovFac(cov(trial),trial,p)

    h<-h+1
    #print(h)

    if (h<11) Recall(M=M,h=h,Y_h=Y_h,U_h=U_h,R_h=R_h,EigenMax_h=EigenMax_h,kp_h=kp_h,ndig=ndig,nOutput=nOutput,p=p,FUN=FUN,...=...)
    else
    {
      #g) for j = 1,...m, calculate the standard deviation syj associated with the average of the estimates yj,
      y<-colMeans(Y_h)
      sy<-apply(Y_h,2,sd)/h

      #h) calculate the counterpart of this statistic for the components of u(y(h)) and for λmax and k(h) p;
      uy<-colMeans(U_h)
      suy<-apply(U_h,2,sd)/h

      EigenMax<-mean(EigenMax_h)
      sEigenMax<-sd(EigenMax_h)/h

      kp<-mean(kp_h)
      skp<-sd(kp_h)/h

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
        Y_h.new<-matrix(nrow=h,ncol=nOutput)
        U_h.new<-matrix(nrow=h,ncol=nOutput)
        R_h.new<-vector(mode="list",length=h)
        EigenMax_h.new<-array(dim=h)
        kp_h.new<-array(dim=h)

        Y_h.new[seq(1,h-1),]<-Y_h
        U_h.new[seq(1,h-1),]<-U_h
        R_h.new[seq(1,h-1)]<-R_h
        EigenMax_h.new[seq(1,h-1)]<-EigenMax_h
        kp_h.new[seq(1,h-1)]<-kp_h

        Y_h<-Y_h.new
        U_h<-U_h.new
        R_h<-R_h.new
        EigenMax_h<-EigenMax_h.new
        kp_h<-kp_h.new

        Recall(M=M,h=h,Y_h=Y_h,U_h=U_h,R_h=R_h,EigenMax_h=EigenMax_h,kp_h=kp_h,ndig=ndig,nOutput=nOutput,p=p,FUN=FUN,... = ...)
      }
      else{
        Ry<-array(dim=c(nOutput,nOutput))
        Ry<-Reduce("+",R_h)/length(R_h)

        return(list(ndig=ndig,h=h,y=y,uy=uy,Ry=Ry,kp=kp))
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