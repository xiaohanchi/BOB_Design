################################################################################
# BOB: Bayesian Optimal Design for Biosimilar Trials with Co-Primary Endpoints #
#                                                                              #
# -Numerical Studies                                                           #
#  |-Fixed Designs                                                             #
#    |-FE                                                                      #
################################################################################

###===========================Simulation Settings============================###
maxnsample=160#max sample size
Tmax=1 #stages
nsample=maxnsample/Tmax#sample size per stage
pR=0.5
pT=0.5
overallmuR=0
overallmuT=0
tau2=0.8^2
rho=0
sn=20000

solvemu<-function(p,x,tau2,rho){
  A=matrix(c(p,(1-p),1,-1),nrow = 2,byrow = T)
  b=c(x,sqrt(tau2)*rho/sqrt(p*(1-p)))
  mu<-solve(A,b)
  return(mu)
}

s_ttest<-function(alpha,nT,nR,meanT,meanR,varT,varR,delta){
  #scaled two sample t test for symmetric delta
  lambda<-nT*nR*delta**2/(nT+nR)
  T_stat<-sqrt(nT*nR*(nT+nR-2)/(nT+nR))*(meanT-meanR)/sqrt(varT*(nT-1)+varR*(nR-1))
  C_value<-sqrt(qf(alpha,1,(nT+nR-2),lambda))
  if(abs(T_stat) < C_value){
    return(1)
  }else{
    return(0)
  }
}

muR<-solvemu(pR,overallmuR,tau2,rho)
sigmaR<-sqrt(tau2-pR*(1-pR)*(muR[1]-muR[2])**2)
muT<-solvemu(pT,overallmuT,tau2,rho)
sigmaT<-sqrt(tau2-pT*(1-pT)*(muT[1]-muT[2])**2)

###========================Simulation Implementation=========================###
for(t in 1:Tmax){
  addsample=nsample*t#added up sample size at the current stage
  set.seed(233+10*t)
  sampleTt<-sapply(1:sn, function(r) rbinom(nsample,1,pT))
  sampleTe<-sapply(1:sn, function(r) rnorm(nsample,(sampleTt[,r]*muT[1]+(1-sampleTt[,r])*muT[2]),sigmaT))
  sampleRt<-sapply(1:sn, function(r) rbinom(nsample,1,pR))
  sampleRe<-sapply(1:sn, function(r) rnorm(nsample,(sampleRt[,r]*muR[1]+(1-sampleRt[,r])*muR[2]),sigmaR))    
  #mu
  datamuT <- lapply(1:sn,function(r) sampleTe[,r])
  datamuR <- lapply(1:sn,function(r) sampleRe[,r])
  xbarT <- sapply(1:sn,function(r) mean(datamuT[[r]]))
  xbarR <- sapply(1:sn,function(r) mean(datamuR[[r]]))
  varT <- sapply(1:sn,function(r) var(datamuT[[r]]))
  varR <- sapply(1:sn,function(r) var(datamuR[[r]]))
  bios_mu<-sapply(1:sn,function(r) s_ttest(0.05,nsample,nsample,xbarT[r],xbarR[r],varT[r],varR[r],0.2/0.5))
  bios<-bios_mu
  power<-sum(bios)/sn
}

cat("The power (type I error rate) of the design FE:",power,"\n",sep = '')
