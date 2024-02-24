#' This provides MLEs for the parameters of a two-way ANCOVA model when the treatment effects of row factor are identical under heterogeneous variances.
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
#' equalmle(Res,Cov,a,b,q)
#' @export
equalmle<-function(Res,Cov,a,b,q)
{
  Y=lapply(Res, function(col)col[!is.na(col)])
  X=lapply(Cov, function(col)col[!is.na(col)])
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
  t10=NULL;t11=NULL
  for(j in 1:b)
  {
    for(h in 1:q)
    {
      T10=xhj_bar[(j-1)*q+h]*Nu0[h,1]
      t10[h]=T10
    }
    T11=sum(t10)
    t11[j]=T11
  }
  beta0=NULL
  for(j in 1:b)
  {
    T12=Yj_bar[j]-t11[j]
    beta0[j]=T12
  }
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
  Nun0=Nu0
  S0=t13
  repeat
  {
    Un=N/S0
    umn=Un*yM
    t14=NULL;t15=NULL;t16=NULL;t17=NULL;Y_barn=NULL
    for(j in 1:b)
    {
      for(i in 1:a)
      {
        T16=Un[(i-1)*b+j]
        T17=umn[(i-1)*b+j]
        t14[i]=T16
        t15[i]=T17
      }
      t16[j]=sum(t15)
      t17[j]=sum(t14)
      Y_barn[j]=(t16[j])/t17[j]
    }
    t18=NULL
    for(i in 1:a)
    {
      for(j in 1:b)
      {
        for(h in 1:q)
        {
          T18=(xM[(j-1)*q+(i-1)*(b*q)+h]*Un[(i-1)*b+j])
          t18[(j-1)*q+(i-1)*(b*q)+h]=T18
        }
      }
    }
    t20=NULL
    for(jh in 1:(b*q))
    {
      t19=0
      for(i in 1:a)
      {
        T19=t18[(i-1)*(b*q)+jh]
        t19=t19+T19
      }
      t20[jh]=t19
    }
    x_barn=NULL
    for(j in 1:b)
    {
      for(h in 1:q)
      {
        T20=t20[(j-1)*q+h]/t17[j]
        x_barn[(j-1)*q+h]=T20
      }
    }
    t21=NULL;t22=NULL
    for(j in 1:b)
    {
      for(h in 1:q)
      {
        T21=x_barn[(j-1)*q+h]*Nu0[h,1]
        t21[h]=T21
      }
      T22=sum(t21)
      t22[j]=T22
    }
    betan=NULL
    for(j in 1:b)
    {
      T23=Y_barn[j]-t22[j]
      betan[j]=T23
    }
    t24=NULL;t25=NULL;t26=NULL
    for(h in 1:q)
    {
      for(l in 1:q)
      {
        for(j in 1:b)
        {
          for(i in 1:a)
          {
            T24=X[[(i-1)*(b*q)+(j-1)*q+h]]*X[[(i-1)*(b*q)+(j-1)*q+l]]
            t24[i]=sum(T24)/S0[(i-1)*b+j]
          }
          T25=sum(t24)
          t25[j]=T25
        }
        T26=sum(t25)
        t26[l+(h-1)*q]=T26
      }
    }
    Sxx=matrix(t26,nrow=q,ncol=q,byrow=TRUE)
    t27=NULL;t28=NULL;t29=NULL
    for(h in 1:q)
    {
      for(j in 1:b)
      {
        for(i in 1:a)
        {
          T27=X[[(i-1)*(b*q)+(j-1)*q+h]]*(Y[[(i-1)*b+j]]-betan[j])
          t27[i]=sum(T27)/S0[(i-1)*b+j]
        }
        T28=sum(t27)
        t28[j]=T28
      }
      T29=sum(t28)
      t29[h]=T29
    }
    Sxy=matrix(t29,nrow=q,ncol=1,byrow=TRUE)
    Nun=solve(Sxx)%*%Sxy
    t33=NULL
    for(i in 1:a)
    {
      for(j in 1:b)
      {
        t32=NULL
        for(k in 1:N[(i-1)*b+j])
        {
          t32[k]=0
        }

        for(h in 1:q)
        {
          T31=X[[(j-1)*q+(i-1)*(b*q)+h]]
          T32=T31*Nun[h,1]
          t32=t32+T32
        }
        T33=(Y[[(i-1)*b+j]]-betan[j])-t32
        t33[(i-1)*b+j]=sum(T33^2)/(N[(i-1)*b+j])
      }
    }
    dif3=max(abs(Nu0-Nun));dif4=max(abs(beta0-betan))
    if(dif3<=0.00001&dif4<=0.00001)
    {
      break
    }
    beta0=betan;Nu0=Nun;S0=t33
  }
  mle=c(beta0,Nu0,S0)
  return(mle)
}
