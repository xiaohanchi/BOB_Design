################################################################################
# BOB: Bayesian Optimal Design for Biosimilar Trials with Co-Primary Endpoints #
#                                                                              #
# -Numerical Studies                                                           #
#  --Design Calibration                                                        #
#    ---BOB.avg                                                                #
################################################################################

###===========================Simulation Settings============================###
library(mvtnorm)
maxnsample=300#max sample size
Tmax=4 #stages
nsample=maxnsample/Tmax#sample size per stage
pT=0.5
pR=0.5
overallmuT<- 0
overallmuR<- 0
tau2=0.8^2
rho=0
pn=20000
sn=20000

solvemu<-function(p,x,tau2,rho){
  A=matrix(c(p,(1-p),1,-1),nrow = 2,byrow = T)
  b=c(x,sqrt(tau2)*rho/sqrt(p*(1-p)))
  mu<-solve(A,b)
  return(mu)
}

muR<-solvemu(pR,overallmuR,tau2,rho)
muT<-solvemu(pT,overallmuT,tau2,rho)
sigmaR<-sqrt(tau2-pR*(1-pR)*(muR[1]-muR[2])**2)
sigmaT<-sqrt(tau2-pT*(1-pT)*(muT[1]-muT[2])**2)

###========================Simulation Implementation=========================###
mup_bios<-vector(mode="numeric",length=sn)
for(t in 1:Tmax){
  addsample=nsample*t#added up sample
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
    sampleTt<-rbind(sampleTt,sampleTt2)
    sampleTe<-rbind(sampleTe,sampleTe2)
    sampleRt<-rbind(sampleRt,sampleRt2)
    sampleRe<-rbind(sampleRe,sampleRe2)
  }
  sumyT<-colSums(sampleTt)
  sumyR<-colSums(sampleRt)
  pThat<-sumyT/addsample
  pRhat<-sumyR/addsample
  mean_deltap<-(sumyT-sumyR)/(addsample+2)
  var_deltap<-((sumyT+1)*(addsample-sumyT+1)+(sumyR+1)*(addsample-sumyR+1))/((addsample+2)^2*(addsample+3))
  #mu
  xbarT<-sapply(1:sn,function(r) mean(sampleTe[,r]))
  varT<-sapply(1:sn,function(r) var(sampleTe[,r]))
  xbarR<-sapply(1:sn,function(r) mean(sampleRe[,r]))
  varR<-sapply(1:sn,function(r) var(sampleRe[,r]))
  mean_deltamu<-xbarT-xbarR
  var_deltamu<-(addsample-1)/(addsample^2)*(varT+varR)
  
  #mu1T
  datax1T <- lapply(1:sn,function(r) sampleTe[which(sampleTt[,r]==1),r])
  xbar1T<-sapply(1:sn,function(r) mean(datax1T[[r]]))
  #mu2T
  datax2T <- lapply(1:sn,function(r) sampleTe[which(sampleTt[,r]==0),r])
  xbar2T<-sapply(1:sn,function(r) mean(datax2T[[r]]))
  #mu1R
  datax1R <- lapply(1:sn,function(r) sampleRe[which(sampleRt[,r]==1),r])
  xbar1R<-sapply(1:sn,function(r) mean(datax1R[[r]]))
  #mu2R
  datax2R <- lapply(1:sn,function(r) sampleRe[which(sampleRt[,r]==0),r])
  xbar2R<-sapply(1:sn,function(r) mean(datax2R[[r]]))
  rho_deltahat<-(pThat*(1-pThat)*(xbar1T-xbar2T)+pRhat*(1-pRhat)*(xbar1R-xbar2R))/sqrt((pThat*(1-pThat)+pRhat*(1-pRhat))*(varT+varR))
  for(i in 1:sn){
    Sig_joint<-matrix(c(var_deltap[i],rho_deltahat[i]*sqrt(var_deltap[i]*var_deltamu[i]),
                        rho_deltahat[i]*sqrt(var_deltap[i]*var_deltamu[i]),var_deltamu[i]),
                      nrow=2,ncol=2,byrow=T)
    mup_bios[i]<-pmvnorm(lower=c(-0.15,-0.223), upper=c(0.15,0.223), 
                         mean=c(mean_deltap[i],mean_deltamu[i]), 
                         sigma=Sig_joint, alg=Miwa())[1]
  }
  write.table(mup_bios,file = paste("powers",t,".txt",sep = ''))#Save results
}




  