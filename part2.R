library(MASS)
#########Endogeneity######
n=1000
####cor(x1x2) = 0.2, cor(x2,Z) = 0.5
rho12 = 0.2
rho2z = 0.5
Sigma = matrix(c(1,rho12,0,0,rho12,1,0,rho2z,0,0,1,0,0,rho2z,0,1),4,4)
X = mvrnorm(n,rep(0,4),Sigma)
Z = X[,4]
X = X[,-4]

beta = c(2,-1,0.9)

Y1 = X%*%beta+rnorm(n)

lm1 = lm(Y1~X-1) 
summary(lm1)

lm2 = lm(Y1~X[,-1]-1)
summary(lm2)
lm3 = lm(Y1~X[,-3]-1)
summary(lm3)
####################2SLS

sls1 = lm(X[,2]~Z-1)
summary(sls1)
X2hat = predict(sls1)

lm2sls = lm(Y1~X2hat+X[,3]-1)
summary(lm2sls)

#######8.5
dat81=read.table("TableF8-1.csv",header=T,sep=",")
library(MASS)
library(Matrix)
attach(dat81)

N = nrow(dat81)
m0 = lm(LWAGE~EXP+I(EXP^2)+WKS+OCC+IND+SOUTH+SMSA+MS+UNION+ED+FEM+BLK)
summary(m0)

m1 = lm(WKS~LWAGE+ED+UNION+FEM)
summary(m1)
bls = coef(m1)
m01 = lm(LWAGE~IND+ED+UNION+FEM)
m11 = lm(WKS~predict(m01)+ED+UNION+FEM)
summary(m11)
Z = cbind(1,IND,  ED,UNION,FEM)
X = cbind(1,LWAGE,ED,UNION,FEM)
KX = ncol(X)
KZ = ncol(Z)
X_hat = cbind(1,Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%LWAGE,ED,UNION,FEM)


biv = solve(t(X_hat)%*%X_hat)%*%t(X_hat)%*%WKS

solve(t(Z)%*%X)%*%t(Z)%*%WKS




eiv = WKS-X%*%biv

s2iv = t(eiv)%*%eiv/(N-KZ)

asy.var.iv = s2iv[1,1]*solve(t(X)%*%Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%X)
sqrt(diag(asy.var.iv))

iv = ivreg(WKS~LWAGE+ED+UNION+FEM|IND+ED+UNION+FEM)
summary(iv)


########Hausman###########
Mls = diag(N) - X%*%solve(t(X)%*%X)%*%t(X)
Miv = diag(N) - X_hat%*%solve(t(X_hat)%*%X_hat)%*%t(X_hat)
els = Mls%*%WKS
eiv = Miv%*%WKS
s2ls = t(els)%*%els/(N-KX)
s2iv = t(eiv)%*%eiv/(N-KZ)
asy.var.ls = s2ls[1,1]*solve(t(X)%*%X)
asy.var.iv = s2iv[1,1]*solve(t(Z)%*%X)%*%(t(Z)%*%Z)%*%solve(t(X)%*%Z)
asy.var.iv = s2iv[1,1]*solve(t(X_hat)%*%X_hat)
H = t(biv-bls)%*%solve(asy.var.iv-asy.var.ls)%*%(biv-bls)
m = solve(t(X_hat)%*%X_hat)-solve(t(X)%*%X)
rankMatrix(m)
H =  (t(biv-bls)%*%ginv(m)%*%(biv-bls))/s2ls
qchisq(0.95,1)
#############Wu test#############8.6 Wu test
m12 = lm(WKS~LWAGE+ED+UNION+FEM+I(Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%LWAGE))
summary(m12)
###########################8.7
dat87 = read.table("TableF5-2.csv",header=T,sep=",")
attach(dat87)
N = nrow(dat87)
X = cbind(1,REALDPI[-1])
Y = REALCONS[-1]
Z = cbind(1,REALCONS[1:(N-1)],REALDPI[1:(N-1)])

m0 = lm(REALCONS[-1]~REALDPI[-1])
summary(m0)
m01 = lm(REALDPI[-1]~REALCONS[1:(N-1)]+REALDPI[1:(N-1)])
summary(m01)

m11 = lm(REALCONS[-1]~predict(m01))
summary(m11)

m12 = lm(REALCONS[-1]~REALDPI[-1]+predict(m01))
summary(m12)




X_hat = cbind(1,Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%REALDPI[-1])
bls = solve(t(X)%*%X)%*%t(X)%*%Y
biv = solve(t(X_hat)%*%X_hat)%*%t(X_hat)%*%Y

Mls = diag(N-1) - X%*%solve(t(X)%*%X)%*%t(X)
Miv = diag(N-1) - X_hat%*%solve(t(X_hat)%*%X_hat)%*%t(X_hat)
els = Mls%*%Y
eiv = Miv%*%Y
s2ls = t(els)%*%els/(N-1-2)
s2iv = t(eiv)%*%eiv/(N-1-2)
asy.var.ls = s2ls[1,1]*solve(t(X)%*%X)
asy.var.iv = s2iv[1,1]*solve(t(X)%*%Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%X)
#asy.var.iv = s2iv[1,1]*solve(t(X_hat)%*%X_hat)
H = t(biv-bls)%*%solve(asy.var.iv-asy.var.ls)%*%(biv-bls)

m = solve(t(X_hat)%*%X_hat)-solve(t(X)%*%X)

H =  (t(biv-bls)%*%ginv(m)%*%(biv-bls))/s2ls
qchisq(0.95,1)
##WHY£¿£¿
###########################https://bbs.pinggu.org/thread-895167-1-1.html
  
eiv = Y-X%*%biv

s2iv = t(eiv)%*%eiv/(N-1-2)
asy.var.ls = s2ls[1,1]*solve(t(X)%*%X)
asy.var.iv = s2iv[1,1]*solve(t(X)%*%Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%X)
H = t(biv-bls)%*%solve(asy.var.iv-asy.var.ls)%*%(biv-bls)

####################8.8

dat81=read.table("TableF8-1.csv",header=T,sep=",")

attach(dat81)

N = nrow(dat81)
m0 = lm(LWAGE~EXP+I(EXP^2)+WKS+OCC+IND+SOUTH+SMSA+MS+UNION+ED+FEM+BLK)
summary(m0)

m1 = lm(WKS~LWAGE+ED+UNION+FEM)
summary(m1)
bls = coef(m1)
m01 = lm(LWAGE~IND+ED+UNION+FEM+SMSA)
m11 = lm(WKS~predict(m01)+ED+UNION+FEM)
summary(m11)
Z = cbind(1,IND,ED,UNION,FEM,SMSA)
X = cbind(1,LWAGE,ED,UNION,FEM)
KX = ncol(X)
KZ = ncol(Z)
X_hat = cbind(1,Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%LWAGE,ED,UNION,FEM)

bls = solve(t(X)%*%X)%*%t(X)%*%WKS
biv = solve(t(X_hat)%*%X_hat)%*%t(X_hat)%*%WKS

Mls = diag(N) - X%*%solve(t(X)%*%X)%*%t(X)
els = Mls%*%WKS
eiv = WKS-X%*%biv
s2ls = t(els)%*%els/(N)
s2iv = t(eiv)%*%eiv/(N)
est.var.iv = s2iv[1,1]*solve(t(X)%*%Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%X)
sqrt(diag(est.var.iv))
m = t(Z)%*%eiv/N
W = t(m)%*%solve(s2iv[1,1]/N^2*t(Z)%*%Z)%*%m
W = t(m)%*%solve(s2ls[1,1]/N^2*t(Z)%*%Z)%*%m
############8.9

############9.1
############9.1

dat91 = read.table("TableF9-1.csv",header=T,sep=",")
attach(dat91)
N = nrow(dat91[AVGEXP>0,])
m1 = lm(AVGEXP~AGE+OWNRENT+INCOME+I(INCOME^2),data=dat91,subset=AVGEXP>0)
summary(m1)
b=coef(m1)
X_full = cbind(1,AGE,OWNRENT,INCOME,I(INCOME^2))
X = X_full[AVGEXP>0,]
K = ncol(X)
resi = residuals(m1)
plot(INCOME[AVGEXP>0],resi)

s0 = matrix(0,K,K)  #white hetero cov estimator
for (i in 1:N){
temp = resi[i]^2*(X[i,]%*%(t(X[i,])))
s0 = s0+temp
}
S0 = s0/(N)
est.var.white =  N*solve(t(X)%*%X)%*%S0%*%solve(t(X)%*%X)
sqrt(diag(est.var.white))



est.var.DM =  (N)*solve(t(X)%*%X)%*%(s0/(N-K))%*%solve(t(X)%*%X)
sqrt(diag(est.var.DM))

	
m2 = lm(AVGEXP~AGE+OWNRENT,data=dat91,subset=AVGEXP>0)
summary(m2)
##9-2
(0.2436-0.06393)/2/((1-0.2436)/(N-5)) #frac{(r2-r2*)/J}{(1-r2)/(n-k)}
qf(0.99,2,N-5)

f1 = function(b) b[4]
f2 = function(b) b[5]

q = matrix(c(0,0),2,1)
R = rbind(jacobian(f1,b),jacobian(f2,b))
wald = t(R%*%b-q)%*%solve(R%*%est.var.white%*%t(R))%*%(R%*%b-q)
qchisq(0.95,2)
1-pchisq(wald,3)
#############9.3
#white
income2 = INCOME^2
income2 = income2[AVGEXP>0]
age = AGE[AVGEXP>0]
income = INCOME[AVGEXP>0]
own = OWNRENT[AVGEXP>0]
m2 = lm(resi^2~age+own+income+income2+I(age^2)+I(age*own)+I(age*income)
       +I(age*income2)+I(own*income)+I(own*income2)+I(income*income2)+I(income2^2))
summary(m2)

white = N*0.199
1-pchisq(0.95,12)
qchisq(0.95,12)
#BP
ssr = t(resi)%*%resi
gi = resi^2/(ssr[1,1]/72)-1
Z = cbind(1,income,income2)
LM = 0.5*t(gi)%*%Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%gi
qchisq(0.95,2)

############9.4
dat64 = read.table("TableF6-1.csv",header=T,sep=",")
attach(dat64)
N = nrow(dat64)
nT = length(unique(T))
nI = length(unique(I))
T_matrix = rbind(diag(nT-1),0)
T_variable = T_matrix[rep(1:nT,nI),]

I_variable = matrix(0,nrow(dat64),nI-1)
Temp = 1
for (i in 1:(nI-1)){
I_variable[Temp:(Temp+nT-1),i] = i
Temp = Temp+nT
}



Y = log(dat64$C)

m_full = lm(log(C)~log(Q)+I(log(Q)^2)+log(PF),data=dat64)
summary(m_full)
### how to calculate white.asy.var?

b = coef(m_full)
X = cbind(1,log(Q),log(Q)^2,log(PF))
M = diag(nrow(dat64))-X%*%solve(t(X)%*%X)%*%t(X)
e = M%*%Y    #residual(m_full)

ssr = t(e)%*%e
gi = e^2/(ssr[1,1]/90)-1
Z = cbind(1,LF)
LM = 0.5*t(gi)%*%Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%gi
qchisq(0.95,1)

m2 = lm(log(e^2)~LF)
summary(m2)
w = exp(-coef(m2)[1]-coef(m2)[2]*LF)  #weight = 1/exp(xxx)

m_full = lm(log(C)~log(Q)+I(log(Q)^2)+log(PF),data=dat64,weight=w)  #WLS
summary(m_full)


omega = diag(exp(fitted.values(m2)))
bgls = solve(t(X)%*%solve(omega)%*%X)%*%t(X)%*%solve(omega)%*%log(C)

P = chol(solve(omega))
Xstar = P%*%X
M = diag(N)-Xstar%*%solve(t(Xstar)%*%Xstar)%*%t(Xstar)
Ystar = P%*%log(C)

e_star = Ystar-Xstar%*%bgls
s_star =  t(e_star)%*%e_star/(N-4)
sd_gls = sqrt(diag(s_star[1,1]*solve(t(X)%*%solve(omega)%*%X)))

m3 = lm(Ystar~Xstar-1)
summary(m3)

###############9.5
dat92 = read.table("TableF9-2.csv",sep=",")
dat92 = dat92[-1,]

#v1:country v2:year, v3:LGASPCAR, v4:LINCOMEP     v5:LRPMG            v6:LCARPCAP

COUNTRY = dat92[,1]
YEAR = dat92[,2]
LGASPCAR = as.numeric(dat92[,3])
LINCOMEP = as.numeric(dat92[,4])
LRPMG = as.numeric(dat92[,5])
LCARPCAP = as.numeric(dat92[,6])
N = nrow(dat92)
nT = length(unique(YEAR))
nI = length(unique(COUNTRY))
T_matrix = rbind(diag(nT-1),0)
T_variable = T_matrix[rep(1:nT,nI),]

I_variable = matrix(0,nrow(dat92),nI)
Temp = 1
for (i in 1:nI){
I_variable[Temp:(Temp+nT-1),i] = 1
Temp = Temp+nT
}
m1 = lm(LGASPCAR ~LINCOMEP +LRPMG +LCARPCAP +I_variable-1)
summary(m1)
X = cbind(LINCOMEP ,LRPMG ,LCARPCAP ,I_variable)
M = diag(N)-X%*%solve(t(X)%*%X)%*%t(X)
e = M%*%LGASPCAR 
ee = t(e)%*%e
X2 = cbind(1,LINCOMEP ,LRPMG ,LCARPCAP)
M2 =  diag(N)-X2%*%solve(t(X2)%*%X2)%*%t(X2)
e2 = M2%*%LGASPCAR 
ee2 = t(e2)%*%e2
F = (ee2-ee)/(18-1)/(ee/(N-3-18))
qf(0.95,17,N-3-18)
plot(rep(seq(1,nI,1),each = nT),e)

##########
m2 = lm(log(e^2)~I_variable-1)
summary(m2)
w = exp(-I_variable%*%coef(m2)) ##average over dummies

m_full = lm(LGASPCAR ~LINCOMEP +LRPMG +LCARPCAP +I_variable-1,weight=w)
summary(m_full)
omega = diag(exp(fitted.values(m2)))
bgls = solve(t(X)%*%solve(omega)%*%X)%*%t(X)%*%solve(omega)%*%LGASPCAR 

###########
sig_g = rep(0,nI)
temp = 1	
for (i in 1:nI){
e_g = e[temp:(temp+nT-1)]
sig_g[temp:(temp+nT-1)] = t(e_g)%*%e_g/nT
temp = temp+nT
}
w2 = 1/sig_g
m_full2 = lm(LGASPCAR ~LINCOMEP +LRPMG +LCARPCAP +I_variable-1,weight=w2)
summary(m_full2)

wtable = cbind(w/sum(w),w2/sum(w2))
####11.1
dat81 = read.table("TableF8-1.csv",header=T,sep=",")
attach(dat81)
m1 = lm(LWAGE~EXP+I(EXP^2)+WKS+OCC+IND+SOUTH+SMSA+MS+UNION+ED+FEM+BLK)
summary(m1)
N = NROW(dat81)
T = 7
Ni = NROW(dat81)/T
X = cbind(1,EXP,I(EXP^2),WKS,OCC,IND,SOUTH,SMSA,MS,UNION,ED,FEM,BLK)
K = NCOL(X)
resi = residuals(m1)
s0_white = matrix(0,K,K)
for (i in 1:N){
temp = resi[i]^2*(X[i,]%*%(t(X[i,])))
s0_white = s0_white+temp
}
S0_white = s0_white/(N)
est.var.white =  N*solve(t(X)%*%X)%*%S0_white%*%solve(t(X)%*%X)
sqrt(diag(est.var.white))


s0_robust = matrix(0,K,K)
XX = matrix(0,K,K)
ind = 1
for (i in 1:Ni){
xe = colSums(resi[ind:(ind+T-1)]*X[ind:(ind+T-1),])
temp1 = xe%*%t(xe)
s0_robust = s0_robust+temp1
temp2 = t(X[ind:(ind+T-1),])%*%X[ind:(ind+T-1),]
XX = XX+temp2
ind = ind+T
}
est.var.robust=  solve(XX)%*%s0_robust%*%solve(XX)
#est.var.robust=  solve(t(X)%*%X)%*%s0_robust%*% solve(t(X)%*%X)
sqrt(diag(est.var.robust))
##########11.2
dat41 = read.table("TableF4-1.csv",header=T,sep=",")
attach(dat41)
N = NROW(dat41)
surface = HEIGHT*WIDTH
aspect = HEIGHT/WIDTH
m1 = lm(log(PRICE)~log(surface)+SIGNED+log(aspect))
summary(m1)
resi = residuals(m1)
G = length(unique(PICTURE))
X = cbind(1,log(surface),SIGNED,log(aspect))
K = NCOL(X)
s0_robust = matrix(0,K,K)
ind = 1
for (i in 1:G){
ng = length(PICTURE[PICTURE==i])
if (ng==1){
xe = resi[ind:(ind+ng-1)]*X[ind:(ind+ng-1),]
} else xe = colSums(resi[ind:(ind+ng-1)]*X[ind:(ind+ng-1),])
temp1 = xe%*%t(xe)
s0_robust = s0_robust+temp1
ind = ind+ng
}
est.var.robust=  solve(t(X)%*%X)%*%((G/(G-1))*s0_robust)%*% solve(t(X)%*%X)
sqrt(diag(est.var.robust))
?????why?????
s0_robust = matrix(0,K,K)
XX = matrix(0,K,K)
ind = 1
for (i in 1:G){
ng = length(PICTURE[unique(PICTURE)[i]])
if (ng==1){
xe = resi[ind:(ind+ng-1)]*X[ind:(ind+ng-1),]
} else xe = colSums(resi[ind:(ind+ng-1)]*X[ind:(ind+ng-1),])
temp1 = xe%*%t(xe)
s0_robust = s0_robust+temp1
X_temp = as.matrix(X[ind:(ind+ng-1),])
temp2 = X_temp%*%t(X_temp)
XX = XX+temp2
ind = ind+ng
}
est.var.robust=  solve(XX)%*%((G/(G-1))*s0_robust)%*% solve(XX)
sqrt(diag(est.var.robust))
##############11.4
dat114 = read.table("TableF6-3.csv",header=T,sep=",")
attach(dat114)
G = length(unique(COUNTRY))
......
##################11.5
dat81 = read.table("TableF8-1.csv",header=T,sep=",")
attach(dat81)
m1 = lm(LWAGE~EXP+I(EXP^2)+WKS+OCC+IND+SOUTH+SMSA+MS+UNION)
summary(m1)
resi0 = residuals(m1)
N = NROW(dat81)
T = 7
Ni = NROW(dat81)/T
X = cbind(EXP,I(EXP^2),WKS,OCC,IND,SOUTH,SMSA,MS,UNION)
X1 = cbind(1,X)
Z = cbind(1,ED,FEM,BLK)
K = NCOL(X)

s0_white = matrix(0,K,K)
for (i in 1:N){
temp = resi0[i]^2*(X[i,]%*%(t(X[i,])))
s0_white = s0_white+temp
}
I_variable = matrix(0,nrow(dat81),Ni)
Temp = 1
for (i in 1:(Ni)){
I_variable[Temp:(Temp+T-1),i] = 1
Temp = Temp+T
}
#I_variable = I_variable[,-1]
m2 = lm(LWAGE~EXP+I(EXP^2)+WKS+OCC+IND+SOUTH+SMSA+MS+UNION+I_variable)
summary(m2)
D = I_variable
MD = diag(N)-D%*%solve(t(D)%*%D)%*%t(D)
b_within = solve(t(X)%*%MD%*%X)%*%t(X)%*%MD%*%LWAGE
ee =  t(MD%*%LWAGE-MD%*%X%*%b_within)%*%(MD%*%LWAGE-MD%*%X%*%b_within)
s2 = ee/(N-K-ncol(D))
est.var = s2[1,1]*solve(t(X)%*%MD%*%X)
sd = sqrt(diag(est.var))
#####11.11_fixed
s0_robust = matrix(0,K,K)
Xd = MD%*%XolSums(resi[ind:(ind+T-1)]*Xd[ind:(ind+T-1),]
resi = MD%*%LWAGE-MD%*%X%*%b_within
XX = matrix(0,K,K)
ind = 1
for (i in 1:Ni){
xe = c)
temp1 = xe%*%t(xe)
s0_robust = s0_robust+temp1
temp2 = t(Xd[ind:(ind+T-1),])%*%Xd[ind:(ind+T-1),]
XX = XX+temp2
ind = ind+T
}
#est.var.robust=  solve(XX)%*%s0_robust%*%solve(XX)  ##need corrections! see footnote
est.var.robust=  solve(XX)%*%((N-1)*Ni/((Ni-1)*(N-K-Ni))*s0_robust)%*%solve(XX)

sqrt(diag(est.var.robust))

#####
a = solve(t(D)%*%D)%*%t(D)%*%(LWAGE-X%*%b_within)
alpha = rep(a,each=T)
m3 = lm(alpha~Z-1)
summary(m3)


m4 = lm(a~Z[seq(1,4159,7),]-1)
summary(m4)
hn = residuals(m4)
h = D%*%hn

m5 = lm(LWAGE~X+Z[,-1]+h)
summary(m5)

###########11.6
m6 = lm(LWAGE~X)
summary(m6)
e2 = residuals(m6)
s22 = t(e2)%*%e2/(N-K-1)
sigu2 = s22-s2
i = matrix(1,T,1)

theta = 1-sqrt(s2[1,1])/sqrt(s2[1,1]+T*sigu2[1,1])


SIGMA = s2[1,1]*diag(T)+sigu2[1,1]*i%*%t(i)
SIGMA_inv = solve(SIGMA)
OMEGA_INV = diag(Ni)%x%SIGMA_inv
bgls = solve(t(X1)%*%OMEGA_INV%*%X1)%*%t(X1)%*%OMEGA_INV%*%LWAGE

#ev = eigen(SIGMA)
#tt = t(ev$vec%*%diag(ev$val^(-0.5))%*%t(ev$vec))
#t(tt)%*%tt

P = chol(OMEGA_INV)

Xstar = P%*%X1
M = diag(N)-Xstar%*%solve(t(Xstar)%*%Xstar)%*%t(Xstar)
Ystar = P%*%LWAGE

e_star = Ystar-Xstar%*%bgls
s_star =  t(e_star)%*%e_star/(N-K-1)
sd_gls = sqrt(diag(s_star[1,1]*solve(t(X1)%*%OMEGA_INV%*%X1)))

m7 = lm(Ystar~Xstar-1)
summary(m7)
resip = residuals(m7)



############9-33/11.11

s0_robust2 =  matrix(0,K+1,K+1)
XX = matrix(0,K+1,K+1)
ind = 1
for (i in 1:Ni){
xe = colSums(resip[ind:(ind+T-1)]*Xstar[ind:(ind+T-1),])
temp1 = xe%*%t(xe)
s0_robust2 = s0_robust2+temp1
temp2 = t(Xstar[ind:(ind+T-1),])%*%Xstar[ind:(ind+T-1),]
XX = XX+temp2
ind = ind+T
}


est.var.robust2 = solve(XX)%*%s0_robust2%*%solve(XX)
sqrt(diag(est.var.robust2))

####################################
ev = eigen(SIGMA)
P2 = t(ev$vec%*%diag(ev$val^(-0.5))%*%t(ev$vec))
P2 = diag(Ni)%x%P2
(P2%*%X1)[1,]


Xtheta = (X1+theta*(MD%*%X1-X1))#/sqrt(s2[1,1])
Ytheta = (LWAGE+theta*(MD%*%LWAGE-LWAGE))#/sqrt(s2[1,1])
mm1 = lm(Ytheta~Xtheta-1)
r1 = residuals(mm1)

s0_robust2 =  matrix(0,K+1,K+1)
XX = matrix(0,K+1,K+1)
ind = 1
for (i in 1:Ni){
xe = colSums(r1[ind:(ind+T-1)]*Xtheta[ind:(ind+T-1),])
temp1 = xe%*%t(xe)
s0_robust2 = s0_robust2+temp1
temp2 = t(Xtheta[ind:(ind+T-1),])%*%Xtheta[ind:(ind+T-1),]
XX = XX+temp2
ind = ind+T
}


est.var.robust2 = solve(XX)%*%s0_robust2%*%solve(XX)
sqrt(diag(est.var.robust2))

##
e2_bar = rep(0,Ni)
ind = 1
for (i in 1:Ni){
temp = e2[ind:(ind+T-1)]
e2_bar[i] = mean(temp)
ind = ind+T
}
Te2 = sum((T*e2_bar)^2)
LM = N/(2*T-1)*(Te2/t(e2)%*%e2-1)^2
################
install.packages("plm")
library(plm)
ID =  rep(seq(1,Ni,1,),each=T)
YEAR = rep(seq(1,7,1),Ni)
dat81$ID = ID
dat81$YEAR = YEAR

m8 = plm(LWAGE~EXP+I(EXP^2)+WKS+OCC+IND+SOUTH+SMSA+MS+UNION,data=dat81,
        index=c("ID","YEAR"),effect="individual",model = "pooling")
summary(m8)


############11.8
H = t(b_within-bgls[-1])%*%solve(est.var-vcov(m7)[-1,-1])%*%(b_within-bgls[-1])

qchisq(0.95,9)
###13.2
n=10000
x = rnorm(n,5,2)

m1 = mean(x)
m2 = mean(x^2)

var = m2-m1^2
var(x)*(n-1)/n
########13.3
library(statmod)
n=1000
mu=2
lambda=5
x = rinvgauss(n,mean=mu,shape = lambda)
invgauss = function(a){

mean = a[1]
shape = a[2]
density = sqrt(shape/(2*pi*x^3))*exp(-shape*(x-mean)^2/(2*mean^2*x))
return(sum(log(density)))
}
a = c(0.5,1)
mle = maxLik(logLik=invgauss , print.level=2,start=a,method="BFGS") #BHHH
summary(mle)

a = mean(x)
b = mean(1/x)
c = mean(x^2)
variance = c-a^2
lambda = a^3/variance

#########13.4
n=1000
x1 = rnorm(n,0,1.5)
x2 = rnorm(n,3,0.8)
x = c(x1,x2)
plot(density(x))

Lambda = function(lambda,mu1,mu2,sig1,sig2,t){
cal = lambda*exp(t*mu1+t^2*sig1^2/2)+(1-lambda)*exp(t*mu2+t^2*sig2^2/2)
return(cal)
}

M = function(t){
cal = mean(exp(t*x))
return(cal)
}

k=50
m = 10
j1 = seq(-25,-1,1)
j2 = seq(1,25,1)
j = c(j1,j2)
theta = 2*pi*j/(25*m)

moment_f = function(b){
f = NULL
lambda = b[1]
mu1 = b[2]
mu2 = b[3]
sig1 = b[4]
sig2 = b[5]
for (t in 1:k){
f[t] = Lambda(b[1],b[2],b[3],b[4],b[5],theta[t])-M(theta[t])
}
return(sum(f^2))
}


coeff = c(0.4,0,2,1,1)
result = optim(coeff,moment_f,control = list(maxit = 20000))
result
#########13.5

library(BB)
library(nleqslv)
library(numDeriv)

rm(list=ls(all=TRUE))
datfc = read.table("TableFC-1.csv",header=T,sep=",")
attach(datfc)
n = nrow(datfc)
m = colSums(cbind(Y,Y^2,log(Y),1/Y))/n

moment_f = function(x,list){
f = NULL
p = x[1]
lambda = x[2]
f[1] = p/lambda-m[1]
f[2] = p*(p+1)/lambda^2-m[2]
f[3] = digamma(p)-log(lambda)-m[3]
f[4] = lambda/(p-1)-m[4]
return(f[list])
}

result12 =  nleqslv(x=c(2,1), fn=moment_f,list=c(1,2))
result12
result23 =  nleqslv(x=c(2,1), fn=moment_f,list=c(2,3))
result23
result13 =  nleqslv(x=c(2,0.1), fn=moment_f,list=c(1,3))
result13


f13 = function(x){
f = NULL
p = x[1]
lambda = x[2]
f[1] = p/lambda
f[2] = digamma(p)-log(lambda)
return(f)
}

J = jacobian(f13,result13$x)
data = cbind(Y,log(Y))
mx = diag(1,n)-matrix(1,n,n)/n
covA = t(data)%*%mx%*%data/(n-1)
covA



delta = solve(t(J)%*%solve(covA)%*%(J))/n


#####################
library(reshape)
library(car)
dta = read.table("exp7_4.txt",sep="")
dta = rename(dta,c(V1="ID",V2 = "female",V3 = "year",V4 = "age", V5 = "hsat",V6="handicap",
      V7="handper",V8="hhninc",V9="hhkids",V10="educ",V11="married",V12 = "haupts",V13="reals",V14="fachhs",
      V15="abitur",V16="univ",V17="working",V18="bluec",V19="whitec",V20="self",V21="beamt",V22="docvis",
      V23="hospvis",V24="public",V25="addon"))

dta = dta[-which(dta[,8]==0),]
n = NROW(dta)
attach(dta)
X = cbind(1,age,age^2,educ,female,I(female*educ),I(age*educ))

nlse = function(b){
beta = b
resi = hhninc-exp(X%*%beta)
S = 0.5*sum(resi^2)
}
m0 = lm(log(hhninc)~age+I(age^2)+educ+female+I(female*educ)+I(age*educ))
summary(m0)
coeff = coef(m0)
result0 = optim(coeff,nlse)
result0
#car
m1 = nls(hhninc~exp(a+b*age+c*I(age^2)+d*educ+e*female+f*I(female*educ)+g*I(age*educ)),data=dta,
      start= list(a=coeff[1],b=coeff[2],c=coeff[3],d=coeff[4],e=coeff[5],f=coeff[6],g=coeff[7]))
summary(m1)


dta2 = dta[year==1988,]
n2 = NROW(dta2)
detach(dta)
attach(dta2)
X = cbind(1,age,educ,female)

m0 = lm(log(hhninc)~age+educ+female)
summary(m0)
coeff = coef(m0)
result0 = optim(coeff,nlse)
result0
mui = exp(X%*%result0$par)
epsi = hhninc-mui
X0 = as.numeric(mui)*X
XX0 = t(X0)%*%X0
est.var.nls = sum(epsi^2)/(n2-4)*solve(XX0)
sqrt(diag(est.var.nls))

m1 = nls(hhninc~exp(a+b*age+d*educ+e*female),data=dta2,
      start= list(a=coeff[1],b=coeff[2],d=coeff[3],e= coeff[4]))
summary(m1)

k=4
mm = function(b){
m = NULL
beta1 = b[1]
beta2 = b[2]
beta3 = b[3]
beta4 = b[4]
epsi = hhninc-exp(X%*%b)
for (i in 1:k){
m[i] = mean(as.numeric(epsi)*X[,k])
}
return(sum(m^2))
}
result1 = optim(coeff,mm)
result1

h = exp(X%*%result1$par)
epsi = hhninc-h
Phin = matrix(0,k,k)
Gn = matrix(0,k,k)
for (i in 1:n2){
temp1 = (as.numeric(epsi[i])*X[i,])%*%t(as.numeric(epsi[i])*X[i,])
temp2 = X[i,]%*%t(-as.numeric(h)[i]*X[i,])
Phin = Phin+temp1
Gn = Gn+temp2
}
Phi = Phin/n2
G = Gn/n2
est.var.mm = solve(G%*%solve(Phi)%*%t(G))/n2
sqrt(diag(est.var.mm))



k=6
X2 = cbind(1,age,educ,female,hsat,married)
mm2 = function(b){
m = NULL
beta1 = b[1]
beta2 = b[2]
beta3 = b[3]
beta4 = b[4]
epsi = hhninc-exp(X%*%b)
for (i in 1:k){
m[i] = mean(as.numeric(epsi)*X2[,k])
}
return(sum(m^2))
}
result2 = optim(coeff,mm2)
result2

