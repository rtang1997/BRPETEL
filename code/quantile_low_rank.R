##Code for simulation of Multiple quantile modeling via reduced-rank regression####
library(Matrix)
library(base)
library(MASS)
library(coda)
library(pracma)
library(quantreg)
meanerror=0
covpp90=0
covpp95=0
ovpp90=0
covpp95=0
covppfnorm95=0
covppfnorm90=0
Total_time=1000

###################Define Loss function and (sub)gradient of Loss function###################

##################################################################################
er<-function(theta)
{q=0
for (j  in 1:p)
{ind=(Y<X%*%(theta[,j]))
z=(1-tauseq[j])*ind+tauseq[j]*(1-ind)
q=q+mean(abs(Y-X%*%(theta[,j]))*z)}
return(q)
}


gr<-function(theta,i,tauseq,V)
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
  Xg[i,]=gr(w,i,tauseq,V)
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
    for(kk in 1:5)
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



###################Projection function##################

###################################################



project<-function(Y) ########retraction###########
 {u=svd(Y)$u[,1:r]
 v=svd(Y)$v[,1:r]
  d=diag(svd(Y)$d[1:r],r)
   return(u%*%d%*%t(v))}


######Find the Riemannian gradient of the objective function in finding the local parametrization######
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

 tangentp<-function(u,v,Y) 
 {return(u%*%Y+Y%*%v-u%*%Y%*%v)}

 
 #########solving unit local parametrization#############
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
 
 
 
 
 ######solve the orthonormal basis of the tangent space of the manifold####
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
 

 
 
  
 


 for(time in 1:Total_time)
 {
   ###################Generate data###########################
   
   ##########################################################
   n=500
 m=3
 p=2
 r=1
 D=m*p
 d=m*r+p*r-r^2 
 tauseq=c(0.2,0.5)
 
 
X=mvrnorm(n,rep(0,m),diag(1,m))
X=pnorm(X)
 
 epsilon=rnorm(n,0,1)
 beta=c(1,2,3)
 gamma=c(1,2,3)

 theta_true=cbind(beta+qnorm(tauseq[1])*gamma,beta+qnorm(tauseq[2])*gamma)
 Y=X%*%beta+X%*%gamma*epsilon
  
 
 
 
 
  
 


#########select the initial state in Markov chain######
rq1 <- rq(Y ~.-1, tau=tauseq, data=data.frame(X,Y))
gg<-rq1$coefficients
gg1<-project(gg)
######## select the covariance matrix I used in RRWM algorithm################## 
 

alpha=2*log(n)
L=2000
start_time <- Sys.time()
accp=0
theta0=gg1
Theta=matrix(0,nrow=L+1,ncol=D)
Theta[1,]=matrix(theta0,ncol=1)
sigma=0.32/sqrt(n/500)
Vp=matrix(nullsp(matrix(Theta[1,],nrow=m)),nrow=D,ncol=D-d)
V=matrix(nullspace(t(Vp)),nrow=D,ncol=d)
lt0=lel(matrix(Theta[1,],ncol=p),rep(0,d),V)
lambda=lt0[1:d]
lt=lt0[d+1]
for(t in 1:L)
{nu1=mvrnorm(1,rep(0,d),sigma^2*diag(1,d))
nu=V%*%nu1
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
  
  u0=exp(ly-alpha*er(matrix(y,ncol=p))-(lt-alpha*er(matrix(Theta[t,],ncol=p))))*exp((t(nu)%*%nu-t(tvprime)%*%tvprime)/(2*sigma^2))
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

}

end_time<- Sys.time()
end_time-start_time
burnin=500
TTheta=Theta[(burnin+1):L,]



Thetap0=apply(TTheta,2,mean)
Thetap=project(matrix(Thetap0,ncol=p))
ytheta=matrix(0,nrow=dim(TTheta)[1],ncol=d)
Vthetap=nullspace(t(nullsp(Thetap)))
Thetap=matrix(Thetap,ncol=1)
for(i in 1:dim(ytheta)[1])
{ytheta[i,]=t(Vthetap)%*%(TTheta[i,]-Thetap)}
Sigmap=cov(ytheta)
I=Vthetap%*%Sigmap%*%t(Vthetap)+0.01*nullspace(t(Vthetap))%*%t(nullspace(t(Vthetap)))
########################################



##############################Generate Markov chain##########################

#######################################################################################################  

 
  
alpha=2*log(n)
L=3000
accp=0
theta0=matrix(gg1,ncol=1)
Theta=matrix(0,nrow=L+1,ncol=D)
theta00=project(matrix(mvrnorm(1,theta0,0.01^2*diag(1,D)),ncol=p))
Theta[1,]=matrix(theta00,ncol=1)
sigma=1 
 
 
Vp=matrix(nullsp(matrix(Theta[1,],nrow=m)),nrow=D,ncol=D-d)
V=matrix(nullspace(t(Vp)),nrow=D,ncol=d)
lt0=lel(matrix(Theta[1,],ncol=p),rep(0,d),V)
lambda=lt0[1:d]
lt=lt0[d+1]
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

}




burnin=500
TTheta=Theta[(burnin+1):L,]

######credible region/interval #########
Thetap0=apply(TTheta,2,mean)
Thetap=project(matrix(Thetap0,ncol=p))
ytheta=matrix(0,nrow=dim(TTheta)[1],ncol=d)
Vthetap=nullspace(t(nullsp(Thetap)))
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
 ###############################################
 Thetaf=numeric()
 for (i in 1:dim(TTheta)[1])
   Thetaf[i]=norm(matrix(TTheta[i,],ncol=2),"F")
 seq=sort(Thetaf)
 l1=seq[length(Thetaf)*0.05]
 u1=seq[length(Thetaf)*0.95]
 l2=seq[length(Thetaf)*0.025]
 u2=seq[length(Thetaf)*0.975]

 true_u=norm(matrix(theta_true,ncol=2),"F")
 if((true_u<=u1)&&(true_u>=l1))
   covppfnorm90=covppfnorm90+1

 if((true_u<=u2)&&(true_u>=l2))
   covppfnorm95=covppfnorm95+1



 if(time%%10==0)
 {print(c(covpp95,covppfnorm95,time))}
 
}

 