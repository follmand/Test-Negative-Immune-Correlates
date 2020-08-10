
#
#
#  This does the example of the paper.   27 November 2017
#
#

n<-150000; IS<-1; IR<-1;
b1<- -1.12; b2<- -1.81;

set.seed(3)
par(mfrow=c(3,1))

# generate  ebola and malaria indicators for those who are UNvaccinated
#
yd0 <- rbinom(n,1,.01)
yc0 <- rbinom(n,1,.005)
IB<- ((yd0==1)&(yc0==1)); BE<-rbinom(n,1,.5); yd0[IB]<-BE[IB]; yc0[IB]<-(1-BE[IB])

#
# generate  ebola and malaria indicators for those who ARE VACCINATED. We make all Xs unique and larger than zero.  Code later requires that x is SORTED!
#
rho<-.7
x<-rnorm(n)*.5 + 3; x <-sort(x);
p<-exp(b1+b2*x)/(1+ exp(b1+b2*x))
yd1 <- rbinom(n,1,p)
yc1 <- rbinom(n,1,.005)
IB<- ((yd1==1)&(yc1==1)); BE<-rbinom(n,1,.5); yd1[IB]<-BE[IB]; yc1[IB]<-(1-BE[IB])

#
#  Now let's select those who DO come in with Ebola or malaria AND randomly assume 1/10 have vaccination 28 days ago
#
ISD0<-as.logical(as.logical(rbinom(n,1,0.25*.2))*yd0)
ISD1<-as.logical(as.logical(rbinom(n,1,0.50*.2))*yd1)
ISC0<-as.logical(as.logical(rbinom(n,1,0.25*.2))*yc0)
ISC1<-as.logical(as.logical(rbinom(n,1,0.50*.2))*yc1)
IV  <- (ISD1|ISC1)
IU  <- (ISD0|ISC0)

#
# Make w and get imputed x
#
w<-rnorm(n)*sqrt(.5^2*(1-rho^2)) + 3 + rho*(x-3);  ha<-glm(x[ISC1]~w[ISC1]);  xhat<-coef(ha)[1]+coef(ha)[2]*w;
xorg<-x

#
#   Now let's the data set for Logistic Regression  using x
#
nx<-length(x[IV])

one<-rep(1,n)
hovu<-glm( c(yd1[IV],yd0[IU])~c(x[IV],rep(0,sum(IU)))+ c(rep(0,sum(IV)), one[IU]),family=binomial)
B<-coef(hovu)
px<-1/(1+exp(-B[1]-B[2]*x[IV]))
p0<-1/(1+exp(-B[1]-B[3]))
plot(x[IV],px,xlim=c(0,5.0),ylim=c(0,1),las=1,type="n",xlab="X",ylab="Prob(Ebola Disease)" )
title("IR=X")
lines(x[IV],px)
points(0,1/(1+exp(-B[1]-B[3])),cex=2 )
lines(c(0,10),c(1/(1+exp(-B[1]-B[3])),1/(1+exp(-B[1]-B[3])) ),lty=2 )
points(x[IV],yd1[IV],pch=18)
points(jitter(rep(0,sum(IU)),20),jitter(yd0[IU],.15))
xs<-c(x[IV][nx/8],x[IV][nx*7/8])
VEX<-c(1-exp(B[2]*xs[1]-B[3]),1-exp(B[2]*xs[2]-B[3]),xs)


hox<-hovu

#
#   Now let's get the data set for regression calibration
#

x<-xhat
hovu<-glm( c(yd1[IV],yd0[IU])~c(x[IV],rep(0,sum(IU)))+ c(rep(0,sum(IV)), one[IU]),family=binomial)
B<-coef(hovu)
x<-sort(x)
px<-1/(1+exp(-B[1]-B[2]*x[IV]))
p0<-1/(1+exp(-B[1]-B[3]))
plot(x[IV],px,xlim=c(0,5.0),ylim=c(0,1),las=1,type="n",xlab=expression(widehat(X)(W)),ylab="Prob(Ebola Disease)")
title(expression(IR==widehat(X)(W)))
lines(x[IV],px)
points(0,1/(1+exp(-B[1]-B[3])),cex=2 )
lines(c(0,10),c(1/(1+exp(-B[1]-B[3])),1/(1+exp(-B[1]-B[3])) ),lty=2 )
points(x[IV],yd1[IV],pch=18)
points(jitter(rep(0,sum(IU)),20),jitter(yd0[IU],.15))

hoh<-hovu
xs<-c(x[IV][nx/8],x[IV][nx*7/8])
VEH<-c(1-exp(B[2]*xs[1]-B[3]),1-exp(B[2]*xs[2]-B[3]),xs)

#
#   Now let's do w directly
#

x<-w
hovu<-glm( c(yd1[IV],yd0[IU])~c(x[IV],rep(0,sum(IU)))+ c(rep(0,sum(IV)), one[IU]),family=binomial)
B<-coef(hovu)
x<-sort(x)
px<-1/(1+exp(-B[1]-B[2]*x[IV]))
p0<-1/(1+exp(-B[1]-B[3]))
plot(x[IV],px,xlim=c(0,5.0),ylim=c(0,1),las=1,type="n",xlab="W",ylab="Prob(Ebola Disease)" )
title("IR=W")
lines(x[IV],px)
points(0,1/(1+exp(-B[1]-B[3])),cex=2 )
lines(c(0,10),c(1/(1+exp(-B[1]-B[3])),1/(1+exp(-B[1]-B[3])) ),lty=2 )
points(x[IV],yd1[IV],pch=18)
points(jitter(rep(0,sum(IU)),20),jitter(yd0[IU],.15))



how<-hovu
xs<-c(x[IV][nx/8],x[IV][nx*7/8])
VEW<-c(1-exp(B[2]*xs[1]-B[3]),1-exp(B[2]*xs[2]-B[3]),xs)

ox<-coef(summary(hox))
oh<-coef(summary(hoh))
ow<-coef(summary(how))

#
#  Stuff to put in the text, note that the numbers change for every seed
#
T1<-table(yd0[IU])
T2<-table(yd1[IV])

#
#  Now bootstrap the regression calibration
#

YY  <- c(yd1[IV],yd0[IU])
XX1 <- c(xhat[IV],rep(0,sum(IU)))
XX2 <- c(rep(0,sum(IV)),one[IU])
XX  <- x[IV]
WW  <- w[IV]


Y1<-yd1[IV]
Y0<-yd0[IU]
X <-xorg[IV]
W <-w[IV]
NV<-sum(IV)
NU<-sum(IU)
SUM <- ISC1+IV
IG<-  c( (SUM[IV]==2))


BOOTH<-matrix(rep(0,1000*3),ncol=3);BOOTW<-matrix(rep(0,1000*3),ncol=3)
for(b in 1:1000){
  IBV  <-sample(c(1:NV),replace=TRUE)
  IBU  <-sample(c(1:NU),replace=TRUE)
#  IBV  <-seq(1:NV)
#  IBU  <-seq(1:NU)
  IGBV <- IG[IBV]
  HAB<-glm(X[IGBV]~W[IGBV]);  XHAT<-coef(HAB)[1]+coef(HAB)[2]*W;
  hhb<-glm(c(Y1[IBV],Y0[IBU]) ~ c(XHAT[IBV],rep(0,NU)) +  c(rep(0,NV),rep(1,NU)), family=binomial)
  BOOTH[b,]<- coef(summary(hhb))[,1]
  hwb<-glm(c(Y1[IBV],Y0[IBU]) ~ c(W[IBV],rep(0,NU)) +  c(rep(0,NV),rep(1,NU)), family=binomial)
  BOOTW[b,]<- coef(summary(hwb))[,1]

}
BSESH<-sqrt(apply(BOOTH,2,var))
BSESW<-sqrt(apply(BOOTW,2,var))

print(ox)
print(oh)
print(ow)
print(T1)
print(T2)
print(BSESH)
print(BSESW)

XL<-sort(x[IV])[sum(IV)*.10]
XH<-sort(x[IV])[sum(IV)*.90]
hox<-glm( c(yd1[IV],yd0[IU])~c(x[IV],rep(0,sum(IU)))+ c(rep(0,sum(IV)), one[IU]),family=binomial)
B<-coef(hox)
print(c(XL,XH,1/(1+exp(-B[1]-B[2]*XL)),1/(1+exp(-B[1]-B[2]*XH))))


XL<-sort(xhat[IV])[sum(IV)*.10]
XH<-sort(xhat[IV])[sum(IV)*.90]
hoxhat<-glm( c(yd1[IV],yd0[IU])~c(xhat[IV],rep(0,sum(IU)))+ c(rep(0,sum(IV)), one[IU]),family=binomial)
B<-coef(hoxhat)
print(c(XL,XH,1/(1+exp(-B[1]-B[2]*XL)),1/(1+exp(-B[1]-B[2]*XH))))



XL<-sort(w[IV])[sum(IV)*.10]
XH<-sort(w[IV])[sum(IV)*.90]
hoxhat<-glm( c(yd1[IV],yd0[IU])~c(w[IV],rep(0,sum(IU)))+ c(rep(0,sum(IV)), one[IU]),family=binomial)
B<-coef(hoxhat)
print(c(XL,XH,1/(1+exp(-B[1]-B[2]*XL)),1/(1+exp(-B[1]-B[2]*XH))))

