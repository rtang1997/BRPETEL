#Code for low rank quantile regression using Parking Birmingham Dataset
library(MASS)
library(quantreg)
library(splines)
library(stats)
#########################Load data#########################################################

#########################################################################################
A=read.csv("dataset.csv",header=TRUE)
I=read.table("covariance.txt")  
I=as.matrix(I)
T=A$LastUpdated
Year=as.numeric(substr(T,1,4))
Month=as.numeric(substr(T,6,7))
Day=as.numeric(substr(T,9,10))
Hour=as.numeric(substr(T,12,13))
minute=as.numeric(substr(T,15,16))
Time=minute+Hour*60
Time=scale(Time)
degree=2
Btime=bs(Time,degree=degree)
DD=data.frame(Btime,A$Occupancy/A$Capacity)
DD=scale(DD)
colnames(DD)[3]="Y"
Ytotal=DD[,c(3)]
Xtotal=cbind(intercept=1,DD[,c(1,2)])
theta_true=matrix(c( -0.32159704,0.3515442,0.4703422, -0.01414606,0.4121063,0.4941414, 0.29351833,0.4664097,0.5104362),ncol=3)

library(Matrix)
library(base)
library(MASS)
library(pracma)
library(quantreg)


##########Define Projection function########
Lgradient<-function(u,v,nu,theta0,theta)
{Yg=matrix(0,nrow=m,ncol=p)
S=u%*%theta+theta%*%v-u%*%theta%*%v-theta0-nu
for (i1 in 1:m)
  for(j1 in 1:p)
  {for(i in 1:m)
    for(j in 1:p)
    { Yg[i1,j1]=Yg[i1,j1]+2*S[i,j]*(u[i,i1]*(j==j1)+v[j1,j]*(i==i1)-u[i,i1]*v[j,j1])
    }
  }

return(Yg)
}



project<-function(Y)
{u=svd(Y)$u[,1:r]
v=svd(Y)$v[,1:r]
d=diag(svd(Y)$d[1:r],r)
return(u%*%d%*%t(v))}

tangentp<-function(u,v,Y)
{return(u%*%Y+Y%*%v-u%*%Y%*%v)}

 



meanerror=0
n=2000
m=3
p=3
r=2
D=m*p
d=m*r+p*r-r^2 
nn=dim(Xtotal)[1]
covpp90=rep(0,D)
covpp95=rep(0,D)
Total_time=100
for (time in 1:Total_time)
{
  ###################Subsample data###################
  
  ##################################################################################
hh=sample(seq(1,nn,1),n,replace=TRUE)
X=Xtotal[hh,]
Y=Ytotal[hh]
 

tauseq=c(0.4,0.5,0.6)

###################Define Loss function and gradient of Loss function###################

##################################################################################
er<-function(theta)
{q=0
for (j  in 1:p)
{tau=tauseq[j] 
for( i in 1:n)
{if(Y[i]<t(theta[,j])%*%X[i,])
{z=1-tau}else{
  z=tau}
  q=q+abs(Y[i]-t(theta[,j])%*%X[i,])*z}
q=q/n}
return(q)}



gr<-function(theta,i,V)
{ss=matrix(0,nrow=m*p,ncol=1)
for (w in 1:p)
{z=Y[i]-X[i,]%*%theta[,w]
if(z==0)
{k=0}else if (z>0){
  k=-1*tauseq[w]
}else{
  k=1-tauseq[w]
}
ss[(m*(w-1)+1):(m*w),1]=k*X[i,]}
return(t(V)%*%ss)

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





nullsp<-function(theta0)
{U=svd(theta0)$u[,1:r]
V=svd(theta0)$v[,1:r]
u=U%*%t(U)
v=V%*%t(V)
W=mvrnorm(p,rep(0,r),diag(1,r))
Z=mvrnorm(m,rep(0,r),diag(1,r))
a=U%*%t(W)+Z%*%t(V)
a0=matrix(a,ncol=1)
dim=d
dima0=rankMatrix(a0)
for(k in 1:1000)
{if(dima0<dim)
{ W=mvrnorm(p,rep(0,r),diag(1,r))
Z=mvrnorm(m,rep(0,r),diag(1,r))
a=U%*%t(W)+Z%*%t(V)
a=matrix(a,ncol=1)
a00=cbind(a0,a)
if(rankMatrix(a00)==dima0+1)
{a0=a00
dima0=dima0+1}

}
  
}
Vxp=nullspace(t(a0))
return(Vxp)}


solvepro<-function(u,v,nu,theta0)
{ theta=project(theta0+matrix(nu,ncol=p))
alpha=0.5
iend=1
flag=0
for(iteration in 1:30)
{if (iend==1)
{gr=Lgradient(u,v,matrix(nu,ncol=p),theta0,theta)
gr1=tangentp(u,v,gr)
theta=project(theta-alpha*gr1)
if(norm(gr1,"F")<=10^(-5))
{flag=1}
}
  
}
theta=matrix(theta,ncol=1)
return(c(theta,flag))
}
 
 
rq1 <- rq(Y1 ~.-1, data=data.frame(X,Y1=Y),tau=tauseq[1])
rq2 <- rq(Y2 ~.-1, data=data.frame(X,Y2=Y),tau=tauseq[2])
rq3 <- rq(Y3 ~.-1, data=data.frame(X,Y3=Y),tau=tauseq[3])
gg<-cbind(rq1$coefficients,rq2$coefficients,rq3$coefficients)
gg1<-project(gg)

l=gg1
for (k in 1:200)
{Vp=matrix(nullsp(matrix(l,nrow=m)),nrow=D,ncol=D-d)
V=matrix(nullspace(t(Vp)),nrow=D,ncol=d)
Xg=matrix(0,nrow=D,ncol=1)
for(i in 1:n)
  Xg=Xg+V%*%gr(l,i,V)/n
l=project(l-0.1*matrix(Xg,ncol=p))
}
gg1=l

 

##############################Generate Markov chain##########################

####################################################################################################### 

 
alpha=2*log(n)
L=4000
accp=0
theta0=matrix(gg1,ncol=1)
Theta=matrix(0,nrow=L+1,ncol=D)
Theta[1,]=matrix(theta0,ncol=1)
sigma=0.8/sqrt(n/2000)
Vp=matrix(nullsp(matrix(Theta[1,],nrow=m)),nrow=D,ncol=D-d)
V=matrix(nullspace(t(Vp)),nrow=D,ncol=d)
lt0=lel(matrix(Theta[1,],ncol=p),rep(0,d),V)
lambda=lt0[1:d]
lt=lt0[d+1]
ltseq=numeric()
for(t in 1:L)
{nu1=mvrnorm(1,rep(0,D),sigma^2*I)
nu=V%*%t(V)%*%nu1
U1=svd(matrix(Theta[t,],ncol=p))$u[,1:r]
V1=svd(matrix(Theta[t,],ncol=p))$v[,1:r]
u1=U1%*%t(U1)
v1=V1%*%t(V1)

aflag=solvepro(u1,v1,nu,matrix(Theta[t,],nrow=m))

if (aflag[m*p+1]==0)
{Theta[t+1,]=Theta[t,]
}else{
  y=aflag[1:(m*p)]
  Vyp=matrix(nullsp(matrix(y,nrow=m)),nrow=D,ncol=D-d)
  Vy=matrix(nullspace(t(Vyp)),nrow=D,ncol=d)
  tvprime=t(Vy)%*%(Theta[t,]-y)
  u=runif(1,0,1)
  ly0=lel(matrix(y,ncol=p),lambda,Vy)
  ly=ly0[d+1]
  
  a1=(t(nu)%*%V%*%solve(t(V)%*%I%*%V)%*%t(V)%*%nu)/(2*sigma^2)
  a2=(t(tvprime)%*%solve(t(Vy)%*%I%*%Vy)%*%tvprime)/(2*sigma^2)
  u0=exp(ly-alpha*er(matrix(y,ncol=p))-(lt-alpha*er(matrix(Theta[t,],ncol=p))))*exp(a1-a2)
  if(u>u0)
  {Theta[t+1,]=Theta[t,]
  }else{
    U1=svd(matrix(y,ncol=p))$u[,1:r]
    V1=svd(matrix(y,ncol=p))$v[,1:r]
    u1=U1%*%t(U1)
    v1=V1%*%t(V1)
    
    raflag=solvepro(u1,v1,Vy%*%tvprime,matrix(y,ncol=p))
    
    
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
ltseq[t]=lt
}


burnin=500
TTheta=Theta[(burnin+1):L,]

Thetap0=apply(TTheta,2,mean)
Thetap=project(matrix(Thetap0,ncol=p))
ytheta=matrix(0,nrow=dim(TTheta)[1],ncol=d)
Vthetap=nullspace(t(nullsp(Thetap)))
Thetap=matrix(Thetap,ncol=1)
 

meanerror=meanerror+norm((matrix(Thetap,ncol=p)-theta_true),"F")/Total_time

########credible interval#############
for(j in 1:D) 
{seq=sort(TTheta[,j])
l1=seq[length(seq)*0.05]
u1=seq[length(seq)*0.95]
l2=seq[length(seq)*0.025]
u2=seq[length(seq)*0.975]

true_u=theta_true[j]
if((true_u<=u1)&&(true_u>=l1))
  covpp90[j]=covpp90[j]+1

if((true_u<=u2)&&(true_u>=l2))
  covpp95[j]=covpp95[j]+1
}


if(time%%10==0)
{print(c(covpp95,covpp90,time))}

}

 