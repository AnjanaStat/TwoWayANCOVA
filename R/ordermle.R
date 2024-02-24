#' This provides MLEs for the parameters of a two-way ANCOVA model when the treatment effects of row factor are in simple order restrictions under heterogeneous variances
#'
#' More detailed description
#'
#' @param Res a real data set
#' @param Cov a real data set
#' @param a a positive integer
#' @param b a positive integer
#' @param q a positive integer
#'
#' @return numeric vector
#'
#' @examples
#' a=3;b=2;k=a*b;q=3
#' N1=c(10,10,10,15,15,15,10,10,10,15,15,15,10,10,10,15,15,15)
#' S=c(1,1,2,1,1,2,3,1,2,4,6,3,2,4,2,3,2,1)
#' g=NULL
#' for(i in 1:(k*q))
#' {
#'  g[[i]]=rnorm(N1[i],0,sqrt(S[i]))
#' }
#' Cov=g
#' G2=NULL
#' N=c(10,15,10,15,10,15);m=c(1,2,3,4,2,1)
#' for(i in 1:k)
#' {
#' G2[[i]]=rnorm(N[i],m[i],sqrt(S[i]))
#' }
#' Res=G2
#' ordermle(Res,Cov,a,b,q)
#' @export
ordermle<-function(Res,Cov,a,b,q)
{
  Y=lapply(Res,function(col)col[!is.na(col)])
  X=lapply(Cov,function(col)col[!is.na(col)])
  N=unlist(rbind(lapply(Y,length)))
  yM=unlist(rbind(lapply(Y,mean)))
  xM=unlist(rbind(lapply(X,mean)))
  Yj_bar=NULL
  t1=NULL
  for(j in 1:b)
  {
    for(i in 1:a)
    {
      T1=yM[(i-1)*b+j]
      t1[i]=T1
    }
    Yj_bar[j]=sum(t1)/a
  }
  Yi_bar=NULL
  for(i in 1:a)
  {
    T2=sum(yM[(b*(i-1)+1):(b*(i-1)+b)])/b
    Yi_bar[i]=T2
  }
  Y_bar=sum(Yi_bar)/a
  xhj_bar=NULL
  for(h in 1:q)
  {

    for(j in 1:b)
    {
      t3=0
      for(i in 1:a)
      {
        T3=xM[(j-1)*q+h+(i-1)*(b*q)]
        t3=t3+T3
      }
      xhj_bar[(j-1)*q+h]=t3/a
    }
  }
  xhi_bar=NULL;x_bar=NULL;T5=NULL
  for(h in 1:q)
  {

    for(i in 1:a)
    {
      t3=0
      for(j in 1:b)
      {
        T3=xM[(i-1)*(b*q)+(j-1)*q+h]
        t3=t3+T3
      }
      xhi_bar[(i-1)*q+h]=t3/b
    }
  }
  x_bar=NULL
  for(h in 1:q)
  {
    t4=0
    for(i in 1:a)
    {
      T4=xhi_bar[(i-1)*q+h]
      t4=t4+T4
    }
    x_bar[h]=t4/a
  }
  t6=NULL
  for(i in 1:a)
  {
    for(j in 1:b)
    {
      for(h in 1:q)
      {
        T6=X[[(j-1)*q+(i-1)*(b*q)+h]]-xhj_bar[(j-1)*q+h]-xhi_bar[(i-1)*q+h]+x_bar[h]
        t6[[(j-1)*q+(i-1)*(b*q)+h]]=T6
      }
    }
  }
  t7=NULL;t8=NULL;t9=NULL
  for(h in 1:q)
  {
    for(l in 1:q)
    {
      for(j in 1:b)
      {
        for(i in 1:a)
        {
          T7=t6[[(i-1)*(b*q)+(j-1)*q+h]]*t6[[(i-1)*(b*q)+(j-1)*q+l]]
          t7[i]=sum(T7)
        }
        T8=sum(t7)
        t8[j]=T8
      }
      T9=sum(t8)
      t9[l+(h-1)*q]=T9
    }
  }
  Sxx=matrix(t9,nrow=q,ncol=q,byrow=TRUE)
  t10=NULL;t11=NULL;t12=NULL
  for(h in 1:q)
  {
    for(j in 1:b)
    {
      for(i in 1:a)
      {
        T10=t6[[(i-1)*(b*q)+(j-1)*q+h]]*(Y[[(i-1)*b+j]]-Yj_bar[j]-Yi_bar[i]+Y_bar)
        t10[i]=sum(T10)
      }
      T11=sum(t10)
      t11[j]=T11
    }
    T12=sum(t11)
    t12[h]=T12
  }
  Sxy=matrix(t12,nrow=q,ncol=1,byrow=TRUE)
  Nu=solve(Sxx)%*%Sxy
  t13=NULL;t14=NULL
  for(j in 1:b)
  {
    for(h in 1:q)
    {
      T13=(xhj_bar[(j-1)*q+h]-x_bar[h])*Nu[h,1]
      t13[h]=T13
    }
    T14=sum(t13)
    t14[j]=T14
  }
  beta=NULL
  for(j in 1:b)
  {
    T15=(Yj_bar[j]-Y_bar)-t14[j]
    beta[j]=T15
  }
  t16=NULL;t17=NULL
  for(i in 1:a)
  {
    for(h in 1:q)
    {
      T16=(xhi_bar[(i-1)*q+h])*Nu[h,1]
      t16[h]=T16
    }
    T17=sum(t16)
    t17[i]=T17
  }
  alpha=NULL
  for(i in 1:a)
  {
    T18=Yi_bar[i]-t17[i]
    alpha[i]=T18
  }
  t20=NULL;t21=NULL
  for(i in 1:a)
  {
    for(j in 1:b)
    {
      t19=NULL
      for(k in 1:N[(i-1)*b+j])
      {
        t19[k]=0
      }
      for(h in 1:q)
      {
        T19=X[[(j-1)*q+(i-1)*(b*q)+h]]-xM[(j-1)*q+(i-1)*(b*q)+h]
        T20=T19*Nu[h,1]
        t19=t19+T20
      }
      T21=(Y[[(i-1)*b+j]]-yM[(i-1)*b+j])-t19
      t21[(i-1)*b+j]=sum(T21^2)/(N[(i-1)*b+j]-q-1)
    }
  }
  S1=t21;Nu0=Nu;alpha0=alpha;beta0=beta
  repeat
  {
    U=N/S1
    t23=NULL;t25=NULL;Mi=NULL;W=NULL;wi=NULL
    for(i in 1:a)
    {
      for(j in 1:b)
      {
        t22=0
        for(h in 1:q)
        {
          T22=xM[(j-1)*q+(i-1)*(b*q)+h]
          T23=T22*Nu0[h,1]
          t22=t22+T23
        }
        T24=U[(i-1)*b+j]*(yM[(i-1)*b+j]-beta0[j]-t22)
        T25=U[(i-1)*b+j]
        wi[j]=T25
        t25[j]=T24
      }
      T26=sum(t25);T27=sum(wi)
      Mi[i]=T26/T27
      W[i]=T27
    }
    alphan=Iso::pava(Mi,W)
    t27=NULL;V=NULL
    for(j in 1:b)
    {
      for(i in 1:a)
      {
        T27=U[(i-1)*b+j]
        t27[i]=T27
      }
      V[j]=sum(t27)
    }
    d=NULL
    for(j in 1:(b-1))
    {
      T28=V[j]
      d[j]=T28
    }
    Unit=V[b]*matrix(1,nrow=b-1,ncol=b-1)
    Q=(diag(d,nrow=b-1,ncol=b-1)+Unit)
    t30=NULL;t31=NULL;yj_barn=NULL
    for(j in 1:(b-1))
    {
      for(i in 1:a)
      {
        t28=0
        t29=0
        for(h in 1:q)
        {
          T28=xM[(i-1)*(b*q)+(j-1)*q+h]
          T29=xM[(i-1)*(b*q)+(b-1)*q+h]
          T30=T28*Nu0[h,1]
          T31=T29*Nu0[h,1]
          t28=t28+T30
          t29=t29+T31
        }
        T32=U[(i-1)*b+j]*(yM[(i-1)*b+j]-alphan[i]-t28)
        t30[i]=T32
        T33=U[(i-1)*b+b]*(yM[(i-1)*b+b]-alphan[i]-t29)
        t31[i]=T33
      }
      T32=sum(t30)-sum(t31)
      yj_barn[j]=T32
    }
    B=solve(Q)%*%yj_barn
    betb=-sum(B[1:b-1])
    betan=c(B,betb)
    t7n=NULL;t8n=NULL;t9n=NULL
    for(h in 1:q)
    {
      for(l in 1:q)
      {
        for(j in 1:b)
        {
          for(i in 1:a)
          {
            T7n=X[[(i-1)*(b*q)+(j-1)*q+h]]*X[[(i-1)*(b*q)+(j-1)*q+l]]
            t7n[i]=sum(T7n)/S1[(i-1)*b+j]
          }
          T8n=sum(t7n)
          t8n[j]=T8n
        }
        T9n=sum(t8n)
        t9n[l+(h-1)*q]=T9n
      }
    }
    Sxxn=matrix(t9n,nrow=q,ncol=q,byrow=TRUE)
    t10n=NULL;t11n=NULL;t12n=NULL
    for(h in 1:q)
    {
      for(j in 1:b)
      {
        for(i in 1:a)
        {
          T10n=X[[(i-1)*(b*q)+(j-1)*q+h]]*(Y[[(i-1)*b+j]]-betan[j]-alphan[i])
          t10n[i]=sum(T10n)/S1[(i-1)*b+j]
        }
        T11n=sum(t10n)
        t11n[j]=T11n
      }
      T12n=sum(t11n)
      t12n[h]=T12n
    }
    Sxyn=matrix(t12n,nrow=q,ncol=1,byrow=TRUE)
    Nun=solve(Sxxn)%*%Sxyn
    t33n=NULL
    for(i in 1:a)
    {
      for(j in 1:b)
      {
        t32n=NULL
        for(k in 1:N[(i-1)*b+j])
        {
          t32n[k]=0
        }

        for(h in 1:q)
        {
          T31n=X[[(j-1)*q+(i-1)*(b*q)+h]]
          T32n=T31n*Nun[h,1]
          t32n=t32n+T32n
        }
        T33n=(Y[[(i-1)*b+j]]-alphan[i]-betan[j])-t32n
        t33n[(i-1)*b+j]=sum(T33n^2)/(N[(i-1)*b+j])
      }
    }
    dif3=max(abs(Nu0-Nun));dif4=max(abs(beta0-betan));dif5=max(abs(alphan-alpha0))
    if(dif3<=0.00001&dif4<=0.00001&dif5<=0.00001)
    {
      break
    }
    beta0=betan;Nu0=Nun;S1=t33n;alpha0=alphan
  }
  mle=c(alphan,beta0,Nu0,S1)
  return(mle)
}
