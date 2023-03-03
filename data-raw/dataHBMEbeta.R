set.seed(123)
m=30
x1=runif(m,0,1)
v.x1=1/rgamma(m,2,1)
x1h=rnorm(m,x1,sqrt(v.x1))
x2=runif(m,0,1)
v.x2=1/rgamma(m,2,5)
x2h=rnorm(m,x2,sqrt(v.x2))
x3=runif(m,0,1)
x4=runif(m,0,1)
b0=b1=b2=b3=b4=0.5
u=rnorm(m,0,1)
pi=rgamma(1,1,0.5)
Mu <- exp(b0+b1*x1h+b2*x2h+b3*x3+b4*x4+u)/(1+exp(b0+b1*x1h+b2*x2h+b3*x3+b4*x4+u))
A=Mu*pi
B=(1-Mu) * pi
Y=rbeta(m,A,B)
MU=A/(A+B)
vardir=A*B/((A+B)^2*(A+B+1))
dataHBMEbeta=as.data.frame(cbind(Y,x1,x2,x3,x4,vardir,v.x1,v.x2))
N=nrow(dataHBMEbeta)
for(i in 1:N){
  if(dataHBMEbeta$Y[i]==1){dataHBMEbeta$Y[i]=0.9999999}
  else if(dataHBMEbeta$Y[i]==0){dataHBMEbeta$Y[i]=0.0000001}
}

usethis::use_data(dataHBMEbeta, overwrite = TRUE)
