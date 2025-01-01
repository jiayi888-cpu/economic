X = matrix(c(15,120,19.31,111.79,99.77,120,1240,164.3,1035.9,875.6,19.31,
          164.3,25.218,148.98,131.22,111.79,1035.9,148.98,953.86,799.02,
           99.770,875.60,131.22,799.02,716.67),5,5)
Y = matrix(c(3.0500,26.004,3.9926,23.521,20.732),5,1)
solve(t(X)%*%X)%*%t(X)%*%Y


##############
test_dta = read.table("tst.csv",sep=",",header=T) #data.frame

cor(test_dta[,-2])
attach(test_dta)

#### 
X1 = cbind(1,T,rGNP,R) #X2£ºinflation; rbind£ºrow  bind
X2 = cbind(P)
y = matrix(test_dta[,1],NROW(test_dta),1)

MX1 = diag(NROW(test_dta))-X1%*%solve(t(X1)%*%X1)%*%t(X1) #I - X%*%(t(X)%*%X)^{-1}%*%t(X); residual maker matrix
y_star = MX1%*%y
x2_star = MX1%*%X2
cor(y_star,x2_star)

X1 = cbind(1,rGNP,R,P)
X2 = cbind(T)
y = matrix(test_dta[,1],NROW(test_dta))

MX1 = diag(NROW(test_dta))-X1%*%solve(t(X1)%*%X1)%*%t(X1)
y_star = MX1%*%y
x2_star = MX1%*%X2
cor(y_star,x2_star)

##########centring matrix

mean_reg = lm(y~1)
summary(mean_reg)
derivation = residuals(mean_reg) #y-y^{hat}; y_hat = xb; e = y-xb
N = length(y)
i = matrix(1,N,1)
M0 = diag(1,N)-1/N*i%*%t(i) #centring matrix
derivation2 = M0%*%y
derivation3 = y-mean(y)
cbind(derivation,derivation2,derivation3)
###########4.1
n=10000
omega = rnorm(n,0,1)
x = rnorm(n,0,1)
sigma = 0.5*omega
y = 0.5+0.5*x+sigma   #DGP: data generation process
data = data.frame(y,x)
lm1 = lm(y~x)
summary(lm1)

b_store = NULL
for (t in 1:1000){                                   #loop 
index = sample(seq(1,n,1),1000)
#sub_reg = lm(y~x-1,data=data,subset =index)
#b_store[t] = coef(sub_reg)
x_sub = x[index]
y_sub = y[index]
b_store[t] = solve(t(x_sub)%*%x_sub)%*%t(x_sub)%*%y_sub #(X'X)-1X'y
}
plot(density(b_store))
hist(b_store)


######
n=10000
var = 5
X = rnorm(n,0,sqrt(var))
t(X)%*%X/n
####
n=10000
X1 = rnorm(n)
X2 = rnorm(n)
X = cbind(X1,X2)
beta = c(0.3,0.6)
Y = X%*%beta+rnorm(n)
result = lm(Y~X)
summary(result)
#var(b|x) = sigma^2*(X'X)^(-1) 
M = (diag(n))-X%*%solve(t(X)%*%X)%*%t(X)
e = M%*%Y

s2 = t(e)%*%e/(n-2)  
s2_unscaled = t(e)%*%e/n
lm_vcov = matrix(rep(s2,4),2,2)* solve(t(X)%*%X) #est.var
sqrt(lm_vcov)

lm_vcov_unscaled = matrix(rep(s2_unscaled,4),2,2)* solve(t(X)%*%X)
sqrt(lm_vcov_unscaled)

################

#install.packages("numDeriv")

require(numDeriv) #library(numDeriv)
require(quantreg)  
dat = read.table("TableF2-2.csv",header=T,sep=",")
attach(dat)
t = NROW(dat)
G = 1000000*GASEXP/(GASP*POP)
model = lm(log(G[2:t])~log(GASP[2:t])+log(INCOME[2:t])+log(PNC[2:t])+log(PUC[2:t])+log(G[1:(t-1)]),data=dat)
summary(model)
s2 = sum(residuals(model)^2)/(t-7)
sqrt(0.00021705)

#model = lm(log(G)~log(GASP)+log(INCOME)+log(PNC)+log(PUC),data=dat)
#summary(model)



res = residuals(model)
s2 = t(res)%*%res/(t-7)
s2 = s2[1,1]
X = cbind(1,log(GASP[2:t]),log(INCOME[2:t]),log(PNC[2:t]),log(PUC[2:t]),log(G[1:(t-1)]))
est.var = solve(t(X)%*%X)*s2 #est.var
sd = sqrt(diag(est.var))

f1 = function(b){
b[2]/(1-b[6])
}
f2 = function(b){
b[3]/(1-b[6])
}

g2 = jacobian(f1,coef(model))
g3 = jacobian(f2,coef(model))

delta1 = g2%*%est.var%*%t(g2)
delta2 = g3%*%est.var%*%t(g3)


#############QR##########
rq_result <- rq(log(G[2:t])~log(GASP[2:t])+log(INCOME[2:t])+log(PNC[2:t])+log(PUC[2:t])+log(G[1:(t-1)]),data=dat, tau=0.5)
summary(rq_result)
#######################

dat2 = read.table("TableF4-1.csv",header=T,sep=",")
attach(dat2)
surface = HEIGHT*WIDTH
aspect = HEIGHT/WIDTH

X = cbind(1,log(surface),aspect)
Y = log(PRICE)
m1 = lm(log(PRICE)~log(surface)+aspect)
summary(m1)

quanreg = function(b){
mu = X%*%b
LAD = sum(abs(Y-mu))
return(LAD)
}
m2 = optim(par=c(0,0,0),quanreg )


rq_result <- rq(log(PRICE)~log(surface)+log(aspect), tau=0.5)
summary(rq_result)
########################################
K = 100
R = 100
b_ols = matrix(0,3,K)
b_lad = matrix(0,3,K)
for (k in 1:K){
obs_list = sample(seq(1,420,1),R,replace=F)
b_ols[,k] = coef(lm(log(PRICE[obs_list])~log(surface[obs_list])+aspect[obs_list]))
b_lad[,k] = coef(rq(log(PRICE[obs_list])~log(surface[obs_list])+aspect[obs_list]), tau=0.5)
}
rowMeans(b_ols)
rowMeans(b_lad)

apply(b_ols,1,sd) #ÏòÁ¿»¯ÔËËã
apply(b_lad,1,sd)


###############P92
n = 10000
x = rnorm(n,0,1)
plot(density(x)) #PDF
critical_value1 = qnorm(0.025,0,1)
critical_value2 = qnorm(0.95+0.025,0,1)
abline(v=critical_value1,col=2)
abline(v=critical_value2,col=2)
#############4.8
qt(0.025,47)
qt(0.975,47)
#############4.10

plot(hist(PRICE,30))
X11()
plot(hist(log(PRICE),30))
m1 = lm(log(PRICE)~log(surface)+aspect)
summary(m1)
vcov(m1)

x0 = c(1,log(25.6*31.9),25.6/31.9)
prediction = coef(m1)%*%x0
mean(log(PRICE))



res = residuals(m1)
s2 = t(res)%*%res/(430-3)
s2 = s2[1,1]
X = cbind(1,log(surface),aspect)
est.var = solve(t(X)%*%X)*s2

var_x0 = s2+x0%*%est.var%*%x0  #4-46
sd_x0 = sqrt(var_x0)
x0_upper = prediction+1.96*sd_x0
x0_lower = prediction-1.96*sd_x0
mean(exp(res)) #smearing estimator
exp(prediction)*mean(exp(res))
exp(x0_lower)
exp(x0_upper)


######grid search
mu0 = prediction
sigma0 = sd_x0
 
LO = exp(x0_lower)
decrement = 0.005
#K = 100
denL=1
denU=0
while(denL>denU){
LO = LO-decrement
denL = dlnorm(LO,mu0,sigma0)
pL = pnorm((log(LO)-mu0)/sigma0) #plnorm(LO,mu0,sigma0),¶ÔÊýÕýÌ¬·Ö²¼¾­¶ÔÊý±ä»»ºó¼´ÎªÕýÌ¬·Ö²¼
UO = exp(sigma0*qnorm(pL+0.95)+mu0) #qlnorm(pL+0.95,mu0,sigma0)

denU = dlnorm(UO,mu0,sigma0)

}
#######Lonley data


dat3 = read.table("TableF4-2.csv",header=T,sep=",")
attach(dat3)
m1 = lm(EMPLOY~I(YEAR-1947)+PRICE+GNP+ARMED)
summary(m1)

R1 = lm(PRICE~I(YEAR-1947)+GNP+ARMED)
1/(1-0.9868) #VIF
#########PCA
datpca = read.table("TableF4-3.csv",header=T,sep=",")
attach(datpca)
m1 = lm(log(BOX)~ACTION+COMEDY+ANIMATED+HORROR+I(MPRATING==1)+
        I(MPRATING==2)+I(MPRATING==3)+log(BUDGET)+SEQUEL+STARPOWR)
summary(m1)
buzz1 = log(ADDICT)
buzz2 = log(CMNGSOON)
buzz3 = log(FANDANGO)
buzz4 = CNTWAIT3
N = nrow(datpca)
M0 = diag(N)-1/N
buzz = cbind(buzz1,buzz2,buzz3,buzz4)
# colMeans(M0%*%buzz)
Z = apply(M0%*%buzz, 2, function(x) x/(sd(x)))
V = 1/(N-1)*t(Z)%*%Z   #why??
#cor(Z)
C = eigen(V)
sum(eigen(V)$values/4)
cor(eigen(V)$vectors)
c1 = C$vectors[,1]
c2 = C$vectors[,2]
Zc1 = Z%*%c1
Zc2 = Z%*%c2
cor(Zc1,Zc2)
cor(c1,c2)
m2 = lm(log(BOX)~ACTION+COMEDY+ANIMATED+HORROR+I(MPRATING==1)+
        I(MPRATING==2)+I(MPRATING==3)+log(BUDGET)+SEQUEL+STARPOWR+Zc1)
summary(m2)
########5.1
qnorm(0.95)
qt(0.95,430-3)
pnorm(1,1.33372,0.09072)
pnorm(1,1.33372,0.09072)
#########5.2
dat4 = read.table("TableF5-1.csv",header=T,sep=",")
attach(dat4)
table(LFP)
table(I(KL6+K618>0))
m1 = lm(log(WHRS*WW)~WA+I(WA^2)+WE+I(KL6+K618>0),data=dat4,subset=LFP==1)
X = cbind(1,WA,I(WA^2),WE,I(KL6+K618>0))[LFP==1,]
Y = log(WHRS*WW)[LFP==1]
M = diag(428)-X%*%solve(t(X)%*%X)%*%t(X)
e = M%*%Y
s2 = t(e)%*%e/(428-5)
s2 = s2[1,1]
est.var = s2*solve(t(X)%*%X)
t_tests = coef(m1)/sqrt(diag(est.var))

f1 = function(b){
return(b[2] - 3*b[4])
}

f2 = function(b){
return(b[2] + b[5])
} 



b = coef(m1)
R = rbind(jacobian(f1,b),jacobian(f2,b))
q = matrix(c(0,-0.15),2,1)
wald = t(R%*%b-q)%*%solve(R%*%solve(t(X)%*%X)%*%t(R))%*%(R%*%b-q)/(s2)
qchisq(0.95,2)
1-pchisq(wald,2)

F = t(R%*%b-q)%*%solve(R%*%est.var%*%t(R))%*%(R%*%b-q)/2
1-pf(F,2,428-5)
############5.3
dat53 = read.table("TableF5-2.csv",header=T,sep=",")
attach(dat53)
#nominal interest rates, it , the rate of inflation,
#pt , (the log of) real output, lnYt , and other factors that trend upward through time,
#embodied in the time trend, t.
t = seq(2,NROW(dat53),1)
m1 = lm(log(REALINVS[-1])~TBILRATE[-1]+INFL[-1]+log(REALGDP[-1])+t)
summary(m1)
X = cbind(1,TBILRATE[-1],INFL[-1],log(REALGDP[-1]),t)
Y = log(REALINVS[-1])
M = diag(NROW(dat53)-1)- X%*%solve(t(X)%*%X)%*%t(X)
e = M%*%Y
s2 = t(e)%*%e/(NROW(dat53)-1-5)
s2 = s2[1,1]
est.var = s2*solve(t(X)%*%X)
b = coef(m1)
plus = b[2]+b[3]
f1 = function(b){
b[2]+b[3]
}
g1 = jacobian(f1,b)
g1_sd = sqrt(g1%*%est.var%*%t(g1))
g1_ttest = plus/g1_sd

f2 = function(b){
b[4]
}

f3 = function(b){
b[5]
}
q = matrix(c(0,1,0),3,1)


R = rbind(jacobian(f1,b),jacobian(f2,b),jacobian(f3,b))

wald = t(R%*%b-q)%*%solve(R%*%solve(t(X)%*%X)%*%t(R))%*%(R%*%b-q)/(s2)
qchisq(0.95,3)
1-pchisq(wald,3)	

F = t(R%*%b-q)%*%solve(R%*%est.var%*%t(R))%*%(R%*%b-q)/3
1-pf(F,3,203-5)
###########5.4
dat54 = read.table("TableF5-3.csv",header=T,sep=",")
attach(dat54)
m54 = lm(log(VALUEADD)~log(LABOR)+log(CAPITAL)+I(0.5*log(LABOR)^2)+I(0.5*log(CAPITAL)^2)
      +I(log(LABOR)*log(CAPITAL)))
N = NROW(dat54)
M0 = diag(N)-1/N
X = cbind(1,log(LABOR),log(CAPITAL),I(0.5*log(LABOR)^2),I(0.5*log(CAPITAL)^2)
      ,I(log(LABOR)*log(CAPITAL)))
M = diag(N) - X%*%solve(t(X)%*%X)%*%t(X)
Y = log(VALUEADD)
e = M%*%Y
ee = t(e)%*%e
R2 = 1-t(e)%*%e/t(Y)%*%M0%*%Y
s2 = ee/(N-6)
est.var = s2[1,1]*solve(t(X)%*%X)


m54_cd = lm(log(VALUEADD)~log(LABOR)+log(CAPITAL))
summary(m54_cd)
X_cd = cbind(1,log(LABOR),log(CAPITAL))
M_cd = diag(N) - X_cd%*%solve(t(X_cd)%*%X_cd)%*%t(X_cd)
e_cd = M_cd%*%Y
ee_cd= t(e_cd)%*%e_cd
R2_cd = 1-t(e_cd)%*%e_cd/t(Y)%*%M0%*%Y
s2_cd = ee_cd/(N-3)
est.var_cd = s2_cd[1,1]*solve(t(X_cd)%*%X_cd)

F = (R2-R2_cd)/3/((1-R2)/(N-6))
F = (ee_cd-ee)/3/(ee/(N-6))
1-pf(F,3,N-6)
qf(0.95,3,N-6)

f1 = function(a){
a[2]+a[3]
}

b = coef(m54_cd)
q = 1
R = jacobian(f1,b)

F = t(R%*%b-q)%*%solve(R%*%est.var_cd%*%t(R))%*%(R%*%b-q)/1
qf(0.95,1,N-3)
1-pf(F,1,N-3)

t2 = (b[2]+b[3]-1)^2/(est.var_cd[2,2]+est.var_cd[3,3]+2*est.var_cd[2,3])


f2 = function(a){
a[4]+a[5]+2*a[6]
}
b = coef(m54)
q = matrix(c(1,0),2,1)
R = rbind(jacobian(f1,b),jacobian(f2,b))

F = t(R%*%b-q)%*%solve(R%*%est.var%*%t(R))%*%(R%*%b-q)/2
qf(0.95,2,N-6)

###########5.5
dat4 = read.table("TableF5-1.csv",header=T,sep=",")
attach(dat4)
table(LFP)
table(I(KL6+K618>0))
m1 = lm(log(WHRS*WW)~WA+I(WA^2)+WE+I(KL6+K618>0),data=dat4,subset=LFP==1)
X = cbind(1,WA,I(WA^2),WE,I(KL6+K618>0))[LFP==1,]
Y = log(WHRS*WW)[LFP==1]
M = diag(428)-X%*%solve(t(X)%*%X)%*%t(X)
e = M%*%Y
s2 = t(e)%*%e/(428-5)
s2 = s2[1,1]
est.var = s2*solve(t(X)%*%X)
t_tests = coef(m1)/sqrt(diag(est.var))

f1 = function(b){
return(b[2])
}

f2 = function(b){
return(b[3])
} 

f3 = function(b){
return(b[4])
} 

f4 = function(b){
return(b[5])
} 


b = coef(m1)
R = rbind(jacobian(f1,b),jacobian(f2,b),jacobian(f3,b),jacobian(f4,b))
q = matrix(rep(0,4),4,1)	


F = t(R%*%b-q)%*%solve(R%*%est.var%*%t(R))%*%(R%*%b-q)/4
qf(0.95,4,428-5)
1-pf(F,4,428-5)

###########5.6

dat56 = read.table("TableF5-2.csv",header=T,sep=",")
attach(dat56)
t = nrow(dat56)
m56 = lm(log(REALCONS[2:t])~log(REALDPI[2:t])+log(REALCONS[1:(t-1)]),data=dat56)
summary(m56)
X = cbind(1,log(REALDPI[2:t]),log(REALCONS[1:(t-1)]))
Y = log(REALCONS[2:t])
M = diag(t-1)-X%*%solve(t(X)%*%X)%*%t(X)
e = M%*%Y
resi = residuals(m56)
s2 = t(e)%*%e/(t-4)
s2 = s2[1,1]
s22 = t(resi)%*%resi/(t-3)
est.var = s2*solve(t(X)%*%X)
sqrt(diag(est.var))
b=coef(m56)
d = function(b){
b[2]/(1-b[3])}
g = jacobian(d,b)
sd = sqrt(g%*%est.var%*%t(g))  #0.0002585,why?
((coef(m56)[2]/(1-coef(m56)[3]))-1)/sd #<1.96 not rejected.
###########6.2

dat2 = read.table("TableF4-1.csv",header=T,sep=",")
attach(dat2)
surface = HEIGHT*WIDTH
aspect = HEIGHT/WIDTH

X = cbind(1,log(surface),aspect,SIGNED)
Y = log(PRICE)
m1 = lm(Y~X-1)
summary(m1)
smearing = mean(exp(residuals(m1)))
Exp_price = exp(fitted.values(m1))*smearing

###########6.3

datpca = read.table("TableF4-3.csv",header=T,sep=",")
attach(datpca)
other = 1-ACTION-ANIMATED-COMEDY-HORROR
X = cbind(1,ACTION,ANIMATED,COMEDY,HORROR,other,I(MPRATING==1),
        I(MPRATING==2),I(MPRATING==3),log(BUDGET),SEQUEL,STARPOWR)

solve(t(X)%*%X)

buzz1 = log(ADDICT)
buzz2 = log(CMNGSOON)
buzz3 = log(FANDANGO)
buzz4 = CNTWAIT3
N = nrow(datpca)
M0 = diag(N)-1/N
buzz = cbind(buzz1,buzz2,buzz3,buzz4)
Z = apply(M0%*%buzz, 2, function(x) x/(sd(x)))
V = 1/(N-1)*t(Z)%*%Z   #why??
C = eigen(V)
sum(eigen(V)$values/4)
cor(eigen(V)$vectors)
c1 = C$vectors[,1]
c2 = C$vectors[,2]
Zc1 = Z%*%c1
Zc2 = Z%*%c2

m2 = lm(log(BOX)~ACTION+COMEDY+ANIMATED+HORROR+I(MPRATING==1)+
        I(MPRATING==2)+I(MPRATING==3)+log(BUDGET)+SEQUEL+STARPOWR+Zc1)
summary(m2)
exp(coef(m2)[2])-1
###########6.4
library(dplyr)
dat64 = read.table("TableF6-1.csv",header=T,sep=",")
attach(dat64)
m1 = lm(log(C)~log(Q)+log(Q)^2+log(PF)+LF+T+I,data=dat64)
summary(m1)

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

m_full = lm(log(C)~log(Q)+I(log(Q)^2)+log(PF)+LF+T_variable+I_variable,data=dat64)
summary(m_full)
b = coef(m_full)
X = cbind(1,log(Q),log(Q)^2,log(PF),LF,T_variable,I_variable)
M = diag(nrow(dat64))-X%*%solve(t(X)%*%X)%*%t(X)
e = M%*%Y
t(e)%*%e

K = ncol(X)
s2 = t(e)%*%e/(nrow(dat64)-K)
s2 = s2[1,1]

m_timeonly = lm(log(C)~log(Q)+I(log(Q)^2)+log(PF)+LF+T_variable,data=dat64)
summary(m_timeonly)
X_timeonly = cbind(1,log(Q),log(Q)^2,log(PF),LF,T_variable)
M = diag(nrow(dat64))-X_timeonly%*%solve(t(X_timeonly)%*%X_timeonly)%*%t(X_timeonly)
e_timeonly = M%*%Y
t(e_timeonly)%*%e_timeonly

F =  ((t(e_timeonly)%*%e_timeonly-t(e)%*%e)/5)/(t(e)%*%e/(nrow(dat64)-K))
qf(0.95,5,nrow(dat64)-K)
1-pf(F,5,nrow(dat64)-K)
R = cbind(matrix(0,5,K-5),diag(5))
q  = matrix(rep(0,5),5,1)
wald = t(R%*%b-q)%*%solve(R%*%solve(t(X)%*%X)%*%t(R))%*%(R%*%b-q)/(s2)
wald/F


m_firmonly = lm(log(C)~log(Q)+I(log(Q)^2)+log(PF)+LF+I_variable,data=dat64)
summary(m_firmonly)
X_firmonly = cbind(1,log(Q),log(Q)^2,log(PF),LF,I_variable)
M = diag(nrow(dat64))-X_firmonly%*%solve(t(X_firmonly)%*%X_firmonly)%*%t(X_firmonly)
e_firmonly = M%*%Y
t(e_firmonly)%*%e_firmonly
F =  ((t(e_firmonly)%*%e_firmonly-t(e)%*%e)/14)/(t(e)%*%e/(nrow(dat64)-K))
qf(0.95,14,nrow(dat64)-K)
1-pf(F,14,nrow(dat64)-K)

m_no = lm(log(C)~log(Q)+I(log(Q)^2)+log(PF)+LF,data=dat64)
summary(m_no)
X_no = cbind(1,log(Q),log(Q)^2,log(PF),LF)
ncol(X)
M = diag(nrow(dat64))-X_no%*%solve(t(X_no)%*%X_no)%*%t(X_no)
e_no = M%*%Y
t(e_no)%*%e_no
F =  ((t(e_no)%*%e_no-t(e)%*%e)/19)/(t(e)%*%e/(nrow(dat64)-K))
qf(0.95,19,nrow(dat64)-K)
1-pf(F,19,nrow(dat64)-K)


##6.6
dat66 = read.table("TableF6-2.csv",header=T,sep=",")
attach(dat66)
N = NROW(dat66)
m660 = lm(log(COST)~log(OUTPUT)+log(KPRICE)+log(LPRICE)+log(FPRICE))
summary(m660)
X = cbind(1,log(OUTPUT),log(KPRICE),log(LPRICE),log(FPRICE))
Y = log(COST)
K = ncol(X)
M0 = diag(N)-1/N
e = M0%*%Y
s2 = t(e)%*%e/(N-K)
est.var = s2[1,1]*solve(t(X)%*%X)
R = matrix(c(0,0,1,1,1),1,K)
q = 1
b = coef(m660)
wald = t(R%*%b-q)%*%solve(R%*%solve(t(X)%*%X)%*%t(R))%*%(R%*%b-q)/(s2)
qchisq(0.95,1)
1-pchisq(wald,1)

m66 = lm(log(COST/FPRICE,10)~log(OUTPUT,10)+log(KPRICE/FPRICE,10)+log(LPRICE/FPRICE,10))
summary(m66)
coef(m66)[2]-1.96*0.016859
coef(m66)[2]+1.96*0.016859
###########6.7
rm(list=ls(all=TRUE))
datfc = read.table("TableFC-1.csv",header=T,sep=",")
attach(datfc)

gamma_den = function(a){
beta = a[1]
rho = a[2]
den = (beta+x)^(-rho)/gamma(rho)*y^(rho-1)*exp(-y/(beta+x))
logL = log(den)
return(sum(logL))
}


x = E
y = Y
n = NROW(datfc)
MLE = maxLik(gamma_den,start=c(1,1),method="BFGS")
summary(MLE)


m = lm(Y~E)
summary(m)
f = function(a){
beta = a[1]/a[2]
return(beta)
}
res = residuals(m)
s2 = t(res)%*%res/(n-2)
s2 = s2[1,1]
X = cbind(1,E)
est.var = solve(t(X)%*%X)*s2 

g = jacobian(f,coef(m))

delta = g%*%est.var%*%t(g)

###########6.9



dat = read.table("TableF2-2.csv",header=T,sep=",")
attach(dat)
t = NROW(dat)
G = 1000000*GASEXP/(POP)
model = lm(log(G)~log(INCOME)+log(GASP)+log(PNC)+log(PUC)+I(YEAR-1952),data=dat)
summary(model)
res = residuals(model)
ee = t(res)%*%res
pre = I(YEAR<1974)
post = I(YEAR>1973)

model1 = lm(log(G)~log(INCOME)+log(GASP)+log(PNC)+log(PUC)+I(YEAR-1952),data=dat,subset = YEAR<1974)
summary(model1)
res1 = residuals(model1)
ee1 = t(res1)%*%res1

model2 = lm(log(G)~log(INCOME)+log(GASP)+log(PNC)+log(PUC)+I(YEAR-1952),data=dat,subset = YEAR>1973)
summary(model2)
res2 = residuals(model2)
ee2 = t(res2)%*%res2	

F1 = (ee-ee1-ee2)/6/((ee1+ee2)/(t-2*6))
qf(0.95,6,t-2*6)


model3 = lm(log(G)~log(INCOME)+log(GASP)+log(PNC)+log(PUC)+I(YEAR-1952)+I(YEAR==1974)+I(YEAR==1975)+I(YEAR==1980)+I(YEAR==1981),data=dat)
res3 = residuals(model3)
ee3 = t(res3)%*%res3
F2 = (ee-ee3)/4/(ee3/(t-6-4)) #why



model4 = lm(log(G)~log(INCOME)+log(GASP)+log(PNC)+log(PUC)+I(YEAR-1952)+I(YEAR>1973)+I(YEAR<1974)-1,data=dat)
res4 = residuals(model4)
ee4 = t(res4)%*%res4
F3 = (ee4-ee1-ee2)/5/((ee1+ee2)/(t-2*6))


model5 = lm(log(G)~log(INCOME)*pre+log(GASP)*pre+log(PNC)+log(PUC)+I(YEAR-1952),data=dat)
res5 = residuals(model5)
ee5 = t(res5)%*%res5
F3 = (ee5-ee1-ee2)/3/((ee1+ee2)/(t-2*6))


s2 = ee/(t-6)
s2 = s2[1,1]
X = cbind(1,log(INCOME),log(GASP),log(PNC),log(PUC),I(YEAR-1952))
est.var = solve(t(X)%*%X)*s2
sd = sqrt(diag(est.var))
f1 = function(b){
b[2]/(1-b[6])
}
f2 = function(b){
b[3]/(1-b[6])
}

g2 = jacobian(f1,coef(model))
g3 = jacobian(f2,coef(model))

delta1 = g2%*%est.var%*%t(g2)
delta2 = g3%*%est.var%*%t(g3)


