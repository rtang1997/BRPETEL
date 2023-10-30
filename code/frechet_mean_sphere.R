##Code for simulation of Frechet mean estimation on Sphere
library(pracma)
library(MASS)
library(movMF)


num_iter=250
covpp90=0
covpp95=0
covpp90_1=0
covpp95_1=0
covpp90_2=0
covpp95_2=0
covpp90_3=0
covpp95_3=0
theta_true=c(0.2820149, 0.5375658,  0.7946639)
 
for(iter in 1:num_iter)
{ ###################Generate data###########################
  
  ##########################################################
  n=500

d=2
D=3
 
XX=mvrnorm(n,c(1,2,3),matrix(c(1,0.5,0.3,0.4,2,0.8,0.8,0.9,3),nrow=3))
for(i in 1:n)
{XX[i,]=XX[i,]/as.numeric(sqrt(t(XX[i,])%*%XX[i,]))}
X=XX
 

alpha=2*log(n)
theta0=rep(sqrt(D)/D,D)
L=3000
burnin=500
sigma=0.05
###################Define Loss function and (sub)gradient of Loss function###################

##################################################################################
er<-function(w)
{k=0
for(i in 1:n)
{k=k+acos(t(w)%*%X[i,])^2
}
k=k/n 
return(k)}

gr<-function(w,i,V)
{ if(t(w)%*%X[i,]==1)
{k=-2*t(V)%*%X[i,]}else
  {k=-2*as.numeric(acos(t(w)%*%X[i,])/sqrt(1-(t(w)%*%X[i,])^2))*t(V)%*%X[i,]}
return(k)}



 



###################Compute ETEL###################

###################################################
lel<-function(w,V)
{H=diag(0,d)
l=rep(0,d)
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
return(sum(lss)-n*log(s))
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
delta=-1*solve(2*t(x)%*%V)%*%(t(x)%*%x-1)
a=a+delta
x=w+V%*%a 
if((t(x)%*%x-1)<=epsilon)
{flag=1
i=10}
}
return(c(a,flag))
}


##############################Generate Markov chain##########################

#######################################################################################################  

accp=0
Theta=matrix(0,nrow=L+1,ncol=D)
Theta[1,]=theta0
for(t in 1:L)
{V=matrix(nullspace(t(Theta[t,])),nrow=D,ncol=d)
Vp=matrix(nullspace(t(V)),nrow=D,ncol=D-d)
v=mvrnorm(1,rep(0,d),sigma^2*diag(1,d))
tv=V%*%v
aflag=project((Theta[t,]+tv),Vp)
if (aflag[D-d+1]==0)
{Theta[t+1,]=Theta[t,]
}else{
  y=Theta[t,]+tv+Vp%*%aflag[1:(D-d)]
  Vy=matrix(nullspace(t(y)),nrow=D,ncol=d)
  Vyp=matrix(nullspace(t(Vy)),nrow=D,ncol=D-d)
  tvprime=t(Vy)%*%(Theta[t,]-y)
  u=runif(1,0,1)
  ly=lel(y,Vy)
  lt=lel(Theta[t,],V)
  u0=exp(ly-alpha*er(y)-(lt-alpha*er(Theta[t,])))*exp((t(v)%*%v-t(tvprime)%*%tvprime)/(2*sigma^2))
  if(u>u0)
  {Theta[t+1,]=Theta[t,]
  }else{
    raflag=project((y+Vy%*%tvprime),Vyp)
    if(raflag[D-d+1]==0)
    {Theta[t+1,]=Theta[t,]
    }else{
      Theta[t+1,]=y
      accp=accp+1
    }
  }
}

}

TTheta=Theta[(burnin+1):L,]
############################credible region##################
Thetap0=apply(TTheta,2,mean)
Thetap=as.numeric(1/sqrt(t(Thetap0)%*%Thetap0))*c(Thetap0)
ytheta=matrix(0,nrow=dim(TTheta)[1],ncol=d)
Vthetap=nullspace(t(Thetap))
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
###################################################

################################################
###########################################
Theta_1=TTheta[,1]
seq=sort(Theta_1)
l1=seq[length(Theta_1)*0.05]
u1=seq[length(Theta_1)*0.95]
l2=seq[length(Theta_1)*0.025]
u2=seq[length(Theta_1)*0.975]

true_u=theta_true[1]
if((true_u<=u1)&&(true_u>=l1))
  covpp90_1= covpp90_1+1

if((true_u<=u2)&&(true_u>=l2))
  covpp95_1=covpp95_1+1

 
 

################################################
###########################################
Theta_2=TTheta[,2]
seq=sort(Theta_2)
l1=seq[length(Theta_2)*0.05]
u1=seq[length(Theta_2)*0.95]
l2=seq[length(Theta_2)*0.025]
u2=seq[length(Theta_2)*0.975]

true_u=theta_true[2]
if((true_u<=u1)&&(true_u>=l1))
  covpp90_2= covpp90_2+1

if((true_u<=u2)&&(true_u>=l2))
  covpp95_2=covpp95_2+1

################################################
###########################################
Theta_3=TTheta[,3]
seq=sort(Theta_3)
l1=seq[length(Theta_3)*0.05]
u1=seq[length(Theta_3)*0.95]
l2=seq[length(Theta_3)*0.025]
u2=seq[length(Theta_3)*0.975]

true_u=theta_true[3]
if((true_u<=u1)&&(true_u>=l1))
  covpp90_3= covpp90_3+1

if((true_u<=u2)&&(true_u>=l2))
  covpp95_3=covpp95_3+1




if(iter%%10==0)
{print(c(covpp95,covpp90,iter))}
}
 