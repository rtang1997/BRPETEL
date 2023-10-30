###code for simulation of Wasserstein barycenter####
library(pracma)
library(MASS)
library(movMF)
library(expm)
library(base)

p=2
D=p^2
d=p*(p+1)/2
q=matrix(0,nrow=D,ncol=D-d)
t=0
theta_true=matrix(c(1.0403298456,0, 0, 1.8964114455),ncol=p)
for(j in 1:(p-1))
  for(i in (j+1):p)
  {t=t+1
   q[((j-1)*p+i),t]=1
   q[((i-1)*p+j),t]=-1}
V=nullspace(t(q))
n=500
meanerror=0
covpp90=0
covpp95=0
covpptr90=0
covpptr95=0
covppeigen95=0
covppeigen90=0
Total_time=1000
###################Define Loss function and (sub)gradient of Loss function###################

##################################################################################
er<-function(w)
{k=0
for(i in 1:n)
{x=matrix(X[i,],ncol=p)
x1=matrix(SX[i,],ncol=p)
dd=sqrt(eigen(x1%*%w%*%x1)$values)
k=k+sum(diag(x))+sum(diag(w))-2*sum(dd)
}
k=k/n 
return(k)}




gr<-function(w,i)
{G=diag(0,p)
x=matrix(X[i,],ncol=p)
x1=matrix(SX[i,],ncol=p)
a=eigen(x1%*%w%*%x1)
ad=a$values
au=a$vectors
B=diag(0,D)
s=0
for(j in 1:p)
  for(m in 1:p)
  {s=s+1
  B[,s]=matrix(x1[m,]%*%t(x1[,j]),ncol=1)}
A=rep(0,D)
for(j in 1:p)
{A1=1/sqrt(ad[j])*t(matrix(au[,j]%*%t(au[,j]),ncol=1))%*%B
A=A+t(A1)}
for (l in 1:p)
  for(j in 1:p)
    G[l,j]=(l==j) 
G=matrix(G,ncol=1)-A
return(t(V)%*%G)}

###################Projection function##################

###################################################
project<-function(w,V)
project1<-function(x)
{x1=V%*%t(V)%*%matrix(x,ncol=1)
x1=matrix(x1,ncol=p)
a=eigen(x1)
return((a$vectors)%*%diag((a$values>0)*(a$values))%*%t(a$vectors))
}



###################Compute ETEL###################

###################################################
lel<-function(w,l0)
{H=diag(0,d)
l=l0
G=rep(0,d)
iend=1
Xg=matrix(0,nrow=n,ncol=d)
for(i in 1:n)
  Xg[i,]=gr(w,i)
a0=numeric()
a1=numeric()
for(i in 1:n)
  a0[i]=exp(t(l)%*%Xg[i,])
for(k in 1:10)
{ if(iend==1)
{for(i in 1:n)
{ H=H+a0[i]*Xg[i,]%*%t(Xg[i,])
G=G+a0[i]*Xg[i,]}
  H=H/n
  G=G/n
  gamma=1
  l1=as.vector(l-gamma*solve(H)%*%G)
  for(i in 1:n)
    a1[i]=exp(t(l1)%*%Xg[i,])
  if(mean(a1)<=mean(a0))
  {l=l1
  a0=a1}else{
    iiend=1
    for(kk in 1:10)
    {if (iiend==1)
    { gamma=gamma*0.5
    l1=as.vector(l-gamma*solve(H)%*%G)
    for(i in 1:n)
      a1[i]=exp(t(l1)%*%Xg[i,])
    if(mean(a1)<=mean(a0))
      iiend=0}
    }
    l=l1
    a0=a1}
  
  if(norm(solve(H)%*%G)<=0.00001)
    iend=0
}
}

lss=numeric()
s=0
for(i in 1:n)
{lss[i]=t(l)%*%Xg[i,]
s=s+exp(lss[i])}
return(c(l,sum(lss)-n*log(s)))
}


for(time in 1:Total_time)
{###################Generate data###########################
  
  ##########################################################
  X=matrix(0,nrow=n,ncol=D)
for(i in 1:n)
{theta=rnorm(1,0,0.3)
 
epsilon=rnorm(p,0,0.3)
U=matrix(c(cos(theta),-sin(theta),sin(theta),cos(theta)),ncol=p)
d1=abs(c(1,2)-epsilon)
X[i,]=matrix(U%*%diag(d1)%*%t(U),ncol=1)}

SX=matrix(0,nrow=n,ncol=D)
for(i in 1:n)
{x=matrix(X[i,],ncol=p)
SX[i,]=matrix(sqrtm(x),ncol=1)}
##############################Generate Markov chain##########################

#######################################################################################################  

 
w=diag(1,2)
 alpha=1
 for(l in 1:50)
  {G=rep(0,d)
   for(j in 1:n)
   G=G+gr(w,j)/n
   w=w-alpha*matrix(V%*%G,ncol=p)
   w=project1(w)
    }

gg1=w
 alpha=2*log(n)
 L=3000
 accp=0
 theta0=gg1
 Theta=matrix(0,nrow=L+1,ncol=D)
 Theta[1,]=matrix(theta0,ncol=1)
 sigma=0.02 
 lt0=lel(matrix(Theta[1,],ncol=p),rep(0,d))
 lambda=lt0[1:d]
 lt=lt0[d+1]
 for(t in 1:L)
 {nu1=mvrnorm(1,rep(0,d),sigma^2*diag(1,d))
 nu=V%*%nu1
 y=Theta[t,]+nu
if ( min(eigen(matrix(y,ncol=p))$values)<0)
 {Theta[t+1,]=Theta[t,]
 }else{
  tvprime=t(V)%*%(Theta[t,]-y)
   u=runif(1,0,1)
   ly0=lel(matrix(y,ncol=p),lambda)
   ly=ly0[d+1]
   u0=exp(ly-alpha*er(matrix(y,ncol=p))-(lt-alpha*er(matrix(Theta[t,],ncol=p))))*exp((t(nu)%*%nu-t(tvprime)%*%tvprime)/(2*sigma^2))
   if(u>u0)
   {Theta[t+1,]=Theta[t,]
   }else{
 Theta[t+1,]=y
       lambda=ly0[1:d]
       lt=ly0[1+d]
       
       accp=accp+1 }
 }}
 
 
 
 
 burnin=500
 TTheta=Theta[(burnin+1):L,]
 
 ######credible region #########
 Thetap0=apply(TTheta,2,mean)
 Thetap=project1(matrix(Thetap0,ncol=p))
 ytheta=matrix(0,nrow=dim(TTheta)[1],ncol=d)
 Vthetap=V
 Thetap=matrix(Thetap,ncol=1)
 for(i in 1:dim(ytheta)[1])
 {ytheta[i,]=t(Vthetap)%*%(TTheta[i,]-Thetap)}
 Sigmap=cov(ytheta)
 iSigmap=solve(Sigmap)
 yytheta=numeric()
 for(i in 1:dim(ytheta)[1])
 {yytheta[i]=t(ytheta[i,])%*%iSigmap%*%ytheta[i,]}
 seq1=sort(yytheta)
 u1=seq1[dim(ytheta)[1]*0.95]
 true_u=t(matrix(theta_true,ncol=1)-Thetap)%*%Vthetap%*%iSigmap%*%t(Vthetap)%*%(matrix(theta_true,ncol=1)-Thetap)
 if(true_u<=u1)
   covpp95=covpp95+1
 u2=seq1[dim(ytheta)[1]*0.9]
 if(true_u<=u2)
   covpp90=covpp90+1
  meanerror=meanerror+norm((matrix(Thetap,ncol=p)-theta_true),"F")/Total_time
  ################################################
 #################################Trace##########
 Thetatr=apply(TTheta[,c(1,4)],1,sum)
 seq=sort(Thetatr)
l1=seq[length(Thetatr)*0.05]
u1=seq[length(Thetatr)*0.95]
l2=seq[length(Thetatr)*0.025]
u2=seq[length(Thetatr)*0.975]
  
  true_u=theta_true[1,1]+theta_true[2,2]
  if((true_u<=u1)&&(true_u>=l1))
    covpptr90= covpptr90+1
  
  if((true_u<=u2)&&(true_u>=l2))
    covpptr95=covpptr95+1
 
  ################################################
  ###############################################
  Thetaop=numeric()
  for (i in 1:dim(TTheta)[1])
  Thetaop[i]=max(eigen(matrix(TTheta[i,],2))$values)
  seq=sort(Thetaop)
  l1=seq[length(Thetaop)*0.05]
  u1=seq[length(Thetaop)*0.95]
  l2=seq[length(Thetaop)*0.025]
  u2=seq[length(Thetaop)*0.975]
  
  true_u=max(eigen(theta_true)$values)
  if((true_u<=u1)&&(true_u>=l1))
    covppeigen90=covppeigen90+1
  
  if((true_u<=u2)&&(true_u>=l2))
    covppeigen95=covppeigen95+1
  
  
 
 if(time%%10==0)
 {print(c(covpp95,covpp90,covpptr95,covpptr90,covppeigen95,covppeigen90,time))}
 
 

 
  }
 
 

 
 