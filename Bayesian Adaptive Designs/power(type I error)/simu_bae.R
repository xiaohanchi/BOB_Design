################################################################################
# BOB: Bayesian Optimal Design for Biosimilar Trials with Co-Primary Endpoints #
#                                                                              #
# -Numerical Studies                                                           #
#  |-Bayesian Adaptive Designs                                                 #
#    |-Empirical power or type I error rate                                    #
#      |-BAE                                                                   #
################################################################################

###===========================Simulation Settings============================###
maxnsample=160#max sample size
Tmax=4 #stages
nsample=maxnsample/Tmax#sample size per stage
pR=0.5
pT=0.3
overallmuR=0
overallmuT=0
tau2=0.8^2
rho=0
lambda=0.9482
gamma=1.06
sn=20000

solvemu<-function(p,x,tau2,rho){
  A=matrix(c(p,(1-p),1,-1),nrow = 2,byrow = T)
  b=c(x,sqrt(tau2)*rho/sqrt(p*(1-p)))
  mu<-solve(A,b)
  return(mu)
}
muR<-solvemu(pR,overallmuR,tau2,rho)
sigmaR<-sqrt(tau2-pR*(1-pR)*(muR[1]-muR[2])**2)
muT<-solvemu(pT,overallmuT,tau2,rho)
sigmaT<-sqrt(tau2-pT*(1-pT)*(muT[1]-muT[2])**2)

###===========================Posterior Function============================###
##asymptotic cdf of muT-muR
Pr_dmu<-function(x){
  mean=xbarT - xbarR
  sd=sqrt((nR-1)*(varT+varR)/nR^2)
  r<-pnorm(q = x,mean = mean,sd = sd)
  return(r)
}

##pdf of (delta*sigmaR)
dist_delta_sigmaR<-function(delta,nR,varR,y){
  alpha=nR/2
  beta=(nR-1)*varR/2
  loga<-alpha*log(beta)
  logb<-2*alpha*log(abs(delta))+log(2)
  if(delta>0){
    logc<-(2*alpha+1)*log(y)
  }else if(delta<0){
    logc<-(2*alpha+1)*log(-y)
  }
  log_val<-loga+logb-logc-log(gamma(alpha))-beta*delta^2/y^2
  val<-exp(log_val)
  return(val)
}

##integrand of upper or lower limit
#delta here = (muT-muR)/sigma0
f_upper<-function(y2){
  (1-Pr_dmu(y2))*dist_delta_sigmaR(delta = 0.4,nR = nR,varR = varR,y = y2)
}

f_lower<-function(y1){
  Pr_dmu(y1)*dist_delta_sigmaR(delta = -0.4,nR = nR,varR = varR,y = y1)
}

###========================Simulation Implementation=========================###
nextsn<-sn
pstop<-c(0,0,0)
for(t in 1:Tmax){
  addsample=nsample*t#added up sample size at the current stage
  Cf<-lambda*(addsample/maxnsample)**gamma
  set.seed(233+10*t)
  if(t==1){
    sampleTt<-sapply(1:sn, function(r) rbinom(nsample,1,pT))
    sampleTe<-sapply(1:sn, function(r) rnorm(nsample,(sampleTt[,r]*muT[1]+(1-sampleTt[,r])*muT[2]),sigmaT))
    sampleRt<-sapply(1:sn, function(r) rbinom(nsample,1,pR))
    sampleRe<-sapply(1:sn, function(r) rnorm(nsample,(sampleRt[,r]*muR[1]+(1-sampleRt[,r])*muR[2]),sigmaR))
  }else{
    sampleTt2<-sapply(1:sn, function(r) rbinom(nsample,1,pT))
    sampleTe2<-sapply(1:sn, function(r) rnorm(nsample,(sampleTt2[,r]*muT[1]+(1-sampleTt2[,r])*muT[2]),sigmaT))
    sampleRt2<-sapply(1:sn, function(r) rbinom(nsample,1,pR))
    sampleRe2<-sapply(1:sn, function(r) rnorm(nsample,(sampleRt2[,r]*muR[1]+(1-sampleRt2[,r])*muR[2]),sigmaR))
    sampleTe<-rbind(sampleTe,sampleTe2)
    sampleRe<-rbind(sampleRe,sampleRe2)
  }
  mu_bios<-vector(mode="numeric",length=sn)
  for(i in 1:sn){
    dataxT <- sampleTe[,i]
    nT<-length(dataxT)
    xbarT<-mean(dataxT)
    varT<-var(dataxT)
        
    dataxR <- sampleRe[,i]
    nR<-length(dataxR)
    xbarR<-mean(dataxR)
    varR<-var(dataxR)
        
    mu_bios[i]<- 1 - integrate(f_upper,0,1)$value - integrate(f_lower,-1,0)$value
  }
  if(t!=1){
    mu_bios[addf]<- -1
  }
  stagef<-which(mu_bios < (Cf-0.000001) & mu_bios != -1)
  mu_bios[stagef]<- -1
  addf<-which(mu_bios == -1)
  if(t!=Tmax){
    pstop[t]<-length(stagef)/sn
  }
  nextstagetrials<-setdiff(c(1:sn),addf)# #of trials
  nextsn<-length(nextstagetrials)
  if(t==Tmax|nextsn==0){
    power<-nextsn/sn
    EN0<-sum(pstop*c(nsample,2*nsample,3*nsample))+t*nsample*(1-sum(pstop))
  }
}

cat("The power (or type I error rate) of the design BAE:",power,"\n",sep = '')
cat("Expected Sample Size(EN):",EN0,"\n",sep = '')
