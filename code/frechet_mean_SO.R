##Code for simulation of Frechet mean estimation on Special Orthogonal Group
library(pracma)
library(MASS)

theta_true=c(0.706856, -0.7073575,  0.7073575,  0.706856)

n=500
D=4
d=1
 
covpp90=0
covpp95=0
meanerror=0
 
num_iter=1000
for(iter in 1:num_iter)
{
  
  ###################Generate data###########################
  
  ##########################################################
X=matrix(0,nrow=n,ncol=D)
theta=numeric()
for(i in 1:n)
{theta[i]=rnorm(1,pi/4,0.5)
X[i,]=c(cos(theta[i]),-sin(theta[i]),sin(theta[i]),cos(theta[i]))
}

 

 
 

###################Define Loss function and (sub)gradient of Loss function###################

##################################################################################

er<-function(w)
{k=0
for(i in 1:n)
{k=k+acos(w[1]*X[i,1]+w[3]*X[i,3])^2
}
k=k/n 
return(k)}

gr<-function(w,i,V)
{x=w[1]*X[i,1]+w[3]*X[i,3]
  if (x==1)
  {return(2*t(V)%*%c(-X[i,1],0,-X[i,3],0))}else
  {return(2*acos(x)/sqrt(1-x^2)*t(V)%*%c(-X[i,1],0,-X[i,3],0))}
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





###################Projection function##################

###################################################
project<-function(w,V)
{a=rep(0,D-d)
i=0
flag=0
epsilon=0.01
nmax=10
for (i in 1:10)
{x=w+V%*%a
Jq=matrix(c(2*x[1],2*x[2],0,0,0,0,2*x[3],2*x[4],x[3],x[4],x[1],x[2]),nrow=D,ncol=D-d)
qx=c((x[1]^2+x[2]^2-1),(x[3]^2+x[4]^2-1),(x[1]*x[3]+x[2]*x[4]))
delta=-1*solve(t(Jq)%*%V)%*%qx
a=a+delta
x=w+V%*%a 
qx=c((x[1]^2+x[2]^2-1),(x[3]^2+x[4]^2-1),(x[1]*x[3]+x[2]*x[4]))
if(sqrt(t(qx)%*%qx)<=epsilon)
{flag=1
i=10}
}
return(c(a,flag))
}



##############################Generate Markov chain##########################

#######################################################################################################  


L=3000
alpha=2*log(n)
theta0=c(cos(0),-sin(0),sin(0),cos(0))
burnin=500
sigma=0.15
accp=0
Theta=matrix(0,nrow=L+1,ncol=D)
Theta[1,]=theta0
x=Theta[1,]
Jq=matrix(c(2*x[1],2*x[2],0,0,0,0,2*x[3],2*x[4],x[3],x[4],x[1],x[2]),nrow=D,ncol=D-d)
V=matrix(nullspace(t(Jq)),nrow=D,ncol=d)
Vp=matrix(nullspace(t(V)),nrow=D,ncol=D-d)
lt0=lel(Theta[1,],rep(0,d),V)
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
  if(y[1]*y[4]-y[2]*y[3]<0)
  {Theta[t+1,]=Theta[t,]
  }else{
    Jqy=matrix(c(2*y[1],2*y[2],0,0,0,0,2*y[3],2*y[4],y[3],y[4],y[1],y[2]),nrow=D,ncol=D-d)
    Vy=matrix(nullspace(t(Jqy)),nrow=D,ncol=d)
    Vyp=matrix(nullspace(t(Vy)),nrow=D,ncol=D-d)
    tvprime=t(Vy)%*%(Theta[t,]-y)
    u=runif(1,0,1)
    ly0=lel(y,lambda,Vy)
    ly=ly0[d+1]
    
    u0=exp(ly-alpha*er(y)-(lt-alpha*er(Theta[t,])))*exp((t(v)%*%v-t(tvprime)%*%tvprime)/(2*sigma^2))
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
  }}

}



TTheta=Theta[(burnin+1):L,]
##############################credible region######################
Thetap0=apply(TTheta,2,mean)
a=Thetap0
theta=1
for(k in 1:1000)
{grad=(a[1]-cos(theta))*sin(theta)+(a[2]+sin(theta))*cos(theta)-cos(theta)*(a[3]-sin(theta))+(a[4]-cos(theta))*sin(theta)
theta=theta-0.1*grad}
Thetap=c(cos(theta),-sin(theta),sin(theta),cos(theta))
ytheta=matrix(0,nrow=dim(TTheta)[1],ncol=d)
x=Thetap
Jq=matrix(c(2*x[1],2*x[2],0,0,0,0,2*x[3],2*x[4],x[3],x[4],x[1],x[2]),nrow=D,ncol=D-d)
Vthetap=matrix(nullspace(t(Jq)),nrow=D,ncol=d)
for(i in 1:dim(ytheta)[1])
{ytheta[i,]=t(Vthetap)%*%(TTheta[i,]-Thetap)}
Sigmap=cov(ytheta)
iSigmap=solve(Sigmap)
yytheta=numeric()
for(i in 1:dim(ytheta)[1])
{yytheta[i]=t(ytheta[i,])%*%iSigmap%*%ytheta[i,]}
seq1=sort(yytheta)
u1=seq1[dim(ytheta)[1]*0.95]
true_u=t(theta_true-Thetap)%*%Vthetap%*%iSigmap%*%t(Vthetap)%*%(theta_true-Thetap)
if(true_u<=u1)
  covpp95=covpp95+1
u2=seq1[dim(ytheta)[1]*0.9]
if(true_u<=u2)
  covpp90=covpp90+1
 
meanerror=meanerror+norm(matrix(Thetap-theta_true,ncol=1),"F")/num_iter
###################################################
 

if(iter%%10==0)
{print(c(covpp95,covpp90,iter))}

}
 