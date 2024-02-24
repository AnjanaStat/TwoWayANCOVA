#' It gives the test statistic value and critical value of the test TM for testing homogeneity of treatment effects against ordered alternatives in two-way ANCOVA model
#'
#' More detailed description
#'
#' @param Res a real data set
#' @param Cov a real data set
#' @param a a positive integer
#' @param b a positive integer
#' @param q a positive integer
#' @param alpha a real number which lies in between 0 and  1
#' @return numeric vector
#'
#' @examples
#' a=3;b=2;k=a*b;q=3;alpha=0.05
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
#' TM(Res,Cov,a,b,q,alpha)
#' @export
TM<-function(Res,Cov,a,b,q,alpha)
{
  fun1<-function(Res,Cov,a,b,q)
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
    Xbar=matrix(xhi_bar,nrow=a,ncol=q,byrow=FALSE)
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
    S1=t21
    t22=NULL;t23=NULL;t24=NULL
    for(h in 1:q)
    {
      for(l in 1:q)
      {
        for(j in 1:b)
        {
          for(i in 1:a)
          {
            T22=t6[[(i-1)*(b*q)+(j-1)*q+h]]*t6[[(i-1)*(b*q)+(j-1)*q+l]]
            t22[i]=S1[(i-1)*b+j]*sum(T22)
          }
          T23=sum(t22)
          t23[j]=T23
        }
        T24=sum(t23)
        t24[l+(h-1)*q]=T24
      }
    }
    S=matrix(t24,nrow=q,ncol=q,byrow=TRUE)
    t25=S1/N
    t26=NULL
    for(i in 1:a)
    {
      result=sum(t25[(b*(i-1)+1):(b*(i-1)+b)])/(b^2)
      t26[i]= result
    }
    S2=diag(t26)+Xbar%*%solve(Sxx)%*%S%*%solve(Sxx)%*%t(Xbar)
    V=NULL
    T=NULL
    for(i in 1:a-1)
    {
      V[i]=sqrt(S2[i,i]+S2[i+1,i+1]-2*S2[i,i+1])
      T[i]=(alpha[i+1]-alpha[i])/V[i]
    }
    TM=max(T)
    return(TM)
  }
  fun2<-function(data2,N,Nu0,a,b,q,S,alpha)
  {
    X=data2
    nux=NULL
    for(i in 1:a)
    {
      for(j in 1:b)
      {
        t12=NULL
        for(k1 in 1:N[(i-1)*b+j])
        {
          t12[k1]=0
        }

        for(h in 1:q)
        {
          T13=X[[(j-1)*q+(i-1)*(b*q)+h]]
          T14=T13*Nu0[h]
          t12=t12+T14
        }
        nux[[(i-1)*b+j]]=t12
      }
    }
    k=(a*b)
    B=1000
    Z<-rep(NA,B)
    for(r in 1:B)
    {
      g=NULL
      for(i in 1:k)
      {
        g[[i]]=rnorm(N[i],nux[[i]],sqrt(S[i]))
      }
      TM=fun1(g,X,a,b,q)
      Z[r]<-TM
    }
    y<-sort(Z,decreasing=FALSE)
    m=(1-alpha)*1000
    c<-y[m]
    return(c)
  }
  Y<-lapply(Res, function(col)col[!is.na(col)])
  X<-lapply(Cov, function(col)col[!is.na(col)])
  N=unlist(rbind(lapply(Y,length)))
  yM=unlist(rbind(lapply(Y,mean)))
  xM=unlist(rbind(lapply(X,mean)))
  NyM=N*yM
  tm14=NULL;tm15=NULL;tm16=NULL;tm17=NULL;Yj_bar=NULL
  for(j in 1:b)
  {
    for(i in 1:a)
    {
      T16=N[(i-1)*b+j]
      T17=NyM[(i-1)*b+j]
      tm14[i]=T16
      tm15[i]=T17
    }
    tm16[j]=sum(tm15)
    tm17[j]=sum(tm14)
    Yj_bar[j]=(tm16[j])/tm17[j]
  }
  tm18=NULL
  for(i in 1:a)
  {
    for(j in 1:b)
    {
      for(h in 1:q)
      {
        T18=(xM[(j-1)*q+(i-1)*(b*q)+h]*N[(i-1)*b+j])
        tm18[(j-1)*q+(i-1)*(b*q)+h]=T18
      }
    }
  }
  tm20=NULL
  for(jh in 1:(b*q))
  {
    tm19=0
    for(i in 1:a)
    {
      T19=tm18[(i-1)*(b*q)+jh]
      tm19=tm19+T19
    }
    tm20[jh]=tm19
  }
  xhj_bar=NULL
  for(j in 1:b)
  {
    for(h in 1:q)
    {
      T20=tm20[(j-1)*q+h]/tm17[j]
      xhj_bar[(j-1)*q+h]=T20
    }
  }
  t3=NULL
  for(i in 1:a)
  {
    for(j in 1:b)
    {
      for(h in 1:q)
      {
        T3=X[[(j-1)*q+(i-1)*(b*q)+h]]-xhj_bar[(j-1)*q+h]
        t3[[(j-1)*q+(i-1)*(b*q)+h]]=T3
      }
    }
  }
  t4=NULL;t5=NULL;t6=NULL
  for(h in 1:q)
  {
    for(l in 1:q)
    {
      for(j in 1:b)
      {
        for(i in 1:a)
        {
          T4=t3[[(i-1)*(b*q)+(j-1)*q+h]]*t3[[(i-1)*(b*q)+(j-1)*q+l]]
          t4[i]=sum(T4)
        }
        T5=sum(t4)
        t5[j]=T5
      }
      T6=sum(t5)
      t6[l+(h-1)*q]=T6
    }
  }
  Sxx=matrix(t6,nrow=q,ncol=q,byrow=TRUE)
  t7=NULL;t8=NULL;t9=NULL
  for(h in 1:q)
  {
    for(j in 1:b)
    {
      for(i in 1:a)
      {
        T7=t3[[(i-1)*(b*q)+(j-1)*q+h]]*(Y[[(i-1)*b+j]]-Yj_bar[j])
        t7[i]=sum(T7)
      }
      T8=sum(t7)
      t8[j]=T8
    }
    T9=sum(t8)
    t9[h]=T9
  }
  Sxy=matrix(t9,nrow=q,ncol=1,byrow=TRUE)
  Nu0=solve(Sxx)%*%Sxy
  t13=NULL
  for(i in 1:a)
  {
    for(j in 1:b)
    {
      t12=NULL
      for(k in 1:N[(i-1)*b+j])
      {
        t12[k]=0
      }

      for(h in 1:q)
      {
        T13=X[[(j-1)*q+(i-1)*(b*q)+h]]-xM[(j-1)*q+(i-1)*(b*q)+h]
        T14=T13*Nu0[h,1]
        t12=t12+T14
      }
      T15=(Y[[(i-1)*b+j]]-yM[(i-1)*b+j])-t12
      t13[(i-1)*b+j]=sum(T15^2)/(N[(i-1)*b+j]-q-1)
    }
  }
  S0=t13
  set.seed(4652)
  crit=fun2(X,N,Nu0,a,b,q,S0,alpha)
  TM=fun1(Y,X,a,b,q)
  result=c(TM,crit)
  print("test statistic value and critical value")
  print(result)
  r1=result[1];r2=result[2]
  if(r1>r2)
  {
    print("Null hypothesis is rejected")
  }
  else
  {
    print("Null hypothesis is not rejected")
  }
}
