##code for simulation of Spectral projectors of the covariance matrix####

library(pracma)
library(MASS)
library(movMF)
library(expm)
library(base)
library(ExtDist)
p=3
r=2
D=p^2
d=p*r-0.5*r*(r+1)-(0.5*r^2-0.5*r)
q=matrix(0,nrow=D,ncol=D-d)
t=0
theta_true=t(matrix(c(   0.98425068, -0.01173059, 0.12395027,
                         -0.01173059, 0.99126269, 0.09232205,
                         0.12395027, 0.09232205,0.02448663),ncol=p))

n=1000
meanerror=0
covpp90=0
covpp95=0
 
covpp95v1=0
covpp90v1=0
covpp95v2=0
covpp90v2=0
covpp95v3=0
covpp90v3=0
Total_time=1000
###################Define Loss function and (sub)gradient of Loss function###################

##################################################################################
er<-function(w)
{ A=w%*%Sigma
  k=-sum(diag(A))
return(k)}

gr<-function(w,i,V)
{G=-X[i,]%*%t(X[i,])
G=matrix(G,ncol=1)
return(t(V)%*%G)}


###################Projection function##################

###################################################
Jacobian<-function(A)
{
  q=matrix(0,nrow=D,ncol=p*(p-1)/2)
  t=0
  for(j in 1:(p-1))
    for(i in (j+1):p)
    {t=t+1
    q[((j-1)*p+i),t]=1
    q[((i-1)*p+j),t]=-1}
  q1=matrix(0,nrow=D,ncol=p*(p+1)/2)
  t=0
  for(j in 1:p)
    for(i in j:p)
    {t=t+1
    B=matrix(0,nrow=p,ncol=p)
    if(i==j)
    {for(k in 1:p)
      B[i,k]=2*A[i,k]-(i==k)
    }else{
      for (k in 1:p)
      {B[i,k]=A[j,k]-(k==j)
      B[j,k]=A[i,k]
      }
    }
    q1[,t]=matrix(B,ncol=1)
    }
  q2=matrix(diag(1,p),ncol=1)
  Q=cbind(q,q1,q2)
  return(Q)
  
}



qq<-function(A)
{ q=numeric()
t=0
for(j in 1:(p-1))
  for(i in (j+1):p)
  {t=t+1
  q[t]=A[i,j]-A[j,i]
  }
q1=numeric()
t=0
A2=A%*%A
for(j in 1:p)
  for(i in j:p)
  {t=t+1
   q1[t]=A2[i,j]-A[i,j]
  }
q2=sum(diag(A))-r
Q=c(q,q1,q2)
return(Q)
}







project<-function(w,V)
{w=matrix(w,ncol=1)
a=rep(0,D-d)
i=0
flag=0
epsilon=0.00001
nmax=10
for (i in 1:50)
{x=w+V%*%a
Jq=Jacobian(matrix(x,ncol=p))
qx=qq(matrix(x,ncol=p))
JqV=t(Jq)%*%V
delta=-1*solve(t(JqV)%*%JqV)%*%t(JqV)%*%matrix(qx,ncol=1)
a=a+delta
if(sqrt(t(qx)%*%qx)<=epsilon)
{flag=1
i=10}
}
return(c(a,flag))
}




###################Compute ETEL###################

###################################################
lel<-function(w,l0,V)
{H=diag(0,d)
l=l0
G=rep(0,d)
iend=1
Xg=matrix(0,nrow=n,ncol=d)
for(i in 1:n)
  Xg[i,]=gr(w,i,V)
a0=numeric()
a1=numeric()
for(i in 1:n)
  a0[i]=exp(t(l)%*%Xg[i,])
for(k in 1:5)
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
#lam<-lambda(w,l0)
s=0
for(i in 1:n)
{lss[i]=t(l)%*%Xg[i,]
s=s+exp(lss[i])}
return(c(l,sum(lss)-n*log(s)))
}




for(time in 1:Total_time)
{
  
###################Generate data###########################
    
    ##########################################################
  
  n=500
X=matrix(0,nrow=n,ncol=p)
for(i in 1:n)
{u=runif(1,0,1)
if(u<=0.5)
{X[i,]=mvrnorm(1,rep(0,p),matrix(c(1,0.1,0.1,0.2,1.2,0.1,0.1,0.1,0.3),ncol=p))}else{
  X[i,]=runif(3,-1,1)
}}

mu=apply(X,2,mean)

Sigma=cov(X)*(n-1)/n+mu%*%t(mu)
up=svd(Sigma)$u
w=up[,1]%*%t(up[,1])+up[,2]%*%t(up[,2]) 

V=nullspace(t(Jacobian(w)))



##############################Generate Markov chain##########################

#######################################################################################################  


L=3000
alpha=2*log(n)
theta0= matrix(w,ncol=1)
burnin=500
sigma=0.15/sqrt(n/500)
accp=0
Theta=matrix(0,nrow=L+1,ncol=D)
Theta[1,]=theta0
 
V=nullspace(t(Jacobian(matrix(theta0,ncol=p))))
Vp=matrix(nullspace(t(V)),nrow=D,ncol=D-d)
lt0=lel(matrix(theta0,ncol=p),rep(0,d),V)
lambda=lt0[1:d]
lt=lt0[d+1]
for(t in 1:L)
{v=mvrnorm(1,rep(0,d),sigma^2*diag(1,d))
tv=V%*%v
aflag=project((Theta[t,]+tv),Vp)
if (aflag[D-d+1]==0)
{Theta[t+1,]=Theta[t,]
}else{ 
  y=Theta[t,]+tv+Vp%*%aflag[1:(D-d)]
   Vy=nullspace(t(Jacobian(matrix(y,ncol=p))))
    if(dim(Vy)[2]!=d)
    {Theta[t+1,]=Theta[t,]
    }else
      {Vyp=matrix(nullspace(t(Vy)),nrow=D,ncol=D-d)
    tvprime=t(Vy)%*%(Theta[t,]-y)
    u=runif(1,0,1)
    ly0=lel(matrix(y,ncol=p),lambda,Vy)
    ly=ly0[d+1]
    
    u0=exp(ly-alpha*er(matrix(y,ncol=p))-(lt-alpha*er(matrix(Theta[t,],ncol=p))))*exp((t(v)%*%v-t(tvprime)%*%tvprime)/(2*sigma^2))
    if(u>u0)
    {Theta[t+1,]=Theta[t,]
    }else{
      raflag=project((y+Vy%*%tvprime),Vyp)
      if(raflag[D-d+1]==0)
      {Theta[t+1,]=Theta[t,]
      }else{
        Theta[t+1,]=y
        lambda=ly0[1:d]
        lt=ly0[1+d]
        V=Vy
        Vp=Vyp
        accp=accp+1
      }
    }
  }
}
}


TTheta=Theta[(burnin+1):L,]
######credible region #########
Thetap0=matrix(apply(TTheta,2,mean),ncol=p)
U=svd(Thetap0)$u
c=c(rep(1,r),rep(0,p-r))
Thetap=matrix(U%*%diag(c)%*%t(U),ncol=1)
ytheta=matrix(0,nrow=dim(TTheta)[1],ncol=d)
x=Thetap
 
Vthetap=nullspace(t(Jacobian(matrix(x,ncol=p))))
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

 

meanerror=meanerror+norm(matrix(matrix(Thetap,ncol=p)-theta_true,ncol=1),"F")/Total_time

 
Thetaop=numeric()
for (i in 1:dim(TTheta)[1])
  Thetaop[i]=TTheta[i,1]
seq=sort(Thetaop)
l1=seq[length(Thetaop)*0.05]
u1=seq[length(Thetaop)*0.95]
l2=seq[length(Thetaop)*0.025]
u2=seq[length(Thetaop)*0.975]

true_u=theta_true[1,1]
if((true_u<=u1)&&(true_u>=l1))
  covpp90v1=covpp90v1+1

if((true_u<=u2)&&(true_u>=l2))
  covpp95v1=covpp95v1+1


###############################################
Thetaop=numeric()
for (i in 1:dim(TTheta)[1])
  Thetaop[i]=TTheta[i,5]
seq=sort(Thetaop)
l1=seq[length(Thetaop)*0.05]
u1=seq[length(Thetaop)*0.95]
l2=seq[length(Thetaop)*0.025]
u2=seq[length(Thetaop)*0.975]

true_u=theta_true[2,2]
if((true_u<=u1)&&(true_u>=l1))
  covpp90v2=covpp90v2+1

if((true_u<=u2)&&(true_u>=l2))
  covpp95v2=covpp95v2+1

###############################################
Thetaop=numeric()
for (i in 1:dim(TTheta)[1])
  Thetaop[i]=TTheta[i,9]
seq=sort(Thetaop)
l1=seq[length(Thetaop)*0.05]
u1=seq[length(Thetaop)*0.95]
l2=seq[length(Thetaop)*0.025]
u2=seq[length(Thetaop)*0.975]

true_u=theta_true[3,3]
if((true_u<=u1)&&(true_u>=l1))
  covpp90v3=covpp90v3+1

if((true_u<=u2)&&(true_u>=l2))
  covpp95v3=covpp95v3+1


if(time%%10==0)
{print(c(covpp95,covpp90,time))
}
}
 