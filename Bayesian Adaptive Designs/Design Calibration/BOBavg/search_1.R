################################################################################
# BOB: Bayesian Optimal Design for Biosimilar Trials with Co-Primary Endpoints #
#                                                                              #
# -Numerical Studies                                                           #
#  |-Bayesian Adaptive Designs                                                 #
#    |-Design Calibration                                                      #
#      |-BOB.avg                                                               #
################################################################################

###===========================Simulation Settings============================###
library(parallel)
library(invgamma)
maxnsample=160#max sample size
Tmax=4 #stages
nsample=maxnsample/Tmax#sample size per stage
pT=0.5
pR=0.5
overallmuT<- 0
overallmuR<- 0
deltamu <- -0.32 #or 0.32
tau2=0.8^2
rho=0
pn=20000
sn=20000
x<-seq(-3,3,by=0.0001)

solvemu<-function(p,x,tau2,rho){
  A=matrix(c(p,(1-p),1,-1),nrow = 2,byrow = T)
  b=c(x,sqrt(tau2)*rho/sqrt(p*(1-p)))
  mu<-solve(A,b)
  return(mu)
}

muR<-solvemu(pR,overallmuR,tau2,rho)
sigmaR<-sqrt(tau2-pR*(1-pR)*(muR[1]-muR[2])**2)
###===========================Posterior Function============================###
posterior_mu1T<-function(mu,n){
  xstd<-(mu-xbar1T[n])/sdT[n]
  r<-dt(xstd,addsample)/sdT[n]
  return(r)
}
posterior_mu2T<-function(mu,n){
  xstd<-(mu-xbar2T[n])/sdT[n]
  r<-dt(xstd,addsample)/sdT[n]
  return(r)
}
posterior_mu1R<-function(mu,n){
  xstd<-(mu-xbar1R[n])/sdR[n]
  r<-dt(xstd,addsample)/sdR[n]
  return(r)
}
posterior_mu2R<-function(mu,n){
  xstd<-(mu-xbar2R[n])/sdR[n]
  r<-dt(xstd,addsample)/sdR[n]
  return(r)
}
###===========================Posterior Function============================###

####Pre-allocated Memory
distpT<-distpR<-pd<-matrix(NA,pn,sn)
dist1T<-dist2T<-dist1R<-dist2R<-distsigmaR<-matrix(NA,pn,sn)
pdbios<-cmud<-list(NA)
length(pdbios)<-length(cmud)<-sn
p_bios<-vector(mode="numeric",length=sn)

####parallel
cores <- detectCores()

###========================Simulation Implementation=========================###
deltap<-runif(sn,-0.2,0.2)
muT<-sapply(1:sn, function(r) solvemu((pT+deltap[r]),(overallmuT+deltamu),tau2,rho))
sigmaT<-sqrt(tau2-(pT+deltap)*(1-(pT+deltap))*(muT[1,]-muT[2,])**2)
for(t in 1:Tmax){
  addsample=nsample*t #added up sample
  set.seed(233+10*t)
  if(t==1){
    sampleTt<-sapply(1:sn, function(r) {
      s<-rbinom(nsample,1,(pT+deltap[r]))
      if(sum(s)==0|sum(s)==1){s[1]<-s[2]<-1}
      else if(sum(s)==nsample|sum(s)==(nsample-1)){s[1]<-s[2]<-0}
      return(s)
    })
    sampleTe<-sapply(1:sn, function(r) rnorm(nsample,(sampleTt[,r]*muT[1,r]+(1-sampleTt[,r])*muT[2,r]),sigmaT[r]))
    sampleRt<-sapply(1:sn, function(r) {
      s<-rbinom(nsample,1,pR)
      if(sum(s)==0|sum(s)==1){s[1]<-s[2]<-1}
      else if(sum(s)==nsample|sum(s)==(nsample-1)){s[1]<-s[2]<-0}
      return(s)
    })
    sampleRe<-sapply(1:sn, function(r) rnorm(nsample,(sampleRt[,r]*muR[1]+(1-sampleRt[,r])*muR[2]),sigmaR))
  }else{
    sampleTt2<-sapply(1:sn, function(r) rbinom(nsample,1,(pT+deltap[r])))
    sampleTe2<-sapply(1:sn, function(r) rnorm(nsample,(sampleTt2[,r]*muT[1,r]+(1-sampleTt2[,r])*muT[2,r]),sigmaT[r]))
    sampleRt2<-sapply(1:sn, function(r) rbinom(nsample,1,pR))
    sampleRe2<-sapply(1:sn, function(r) rnorm(nsample,(sampleRt2[,r]*muR[1]+(1-sampleRt2[,r])*muR[2]),sigmaR))
    sampleTt<-rbind(sampleTt,sampleTt2)
    sampleTe<-rbind(sampleTe,sampleTe2)
    sampleRt<-rbind(sampleRt,sampleRt2)
    sampleRe<-rbind(sampleRe,sampleRe2)
  }
  sumyT<-colSums(sampleTt)
  sumyR<-colSums(sampleRt)
  
  cl <- makeCluster(cores)
  clusterExport(cl, c('pn', 'sumyT', 'sumyR', 'addsample','sn'))
  distpT<-parSapply(cl,1:sn,function(r) rbeta(pn,(sumyT[r]+1),(addsample-sumyT[r]+1)))
  distpR<-parSapply(cl,1:sn,function(r) rbeta(pn,(sumyR[r]+1),(addsample-sumyR[r]+1)))
  stopCluster(cl)
  
  pd<-distpT-distpR
  pdbios<-lapply(1:sn, function(r) which(abs(pd[,r])<0.20))
  p_bios<-sapply(1:sn, function(r) length(pdbios[[r]])/pn)
  
  #mu1T
  datax1T <- lapply(1:sn,function(r) sampleTe[which(sampleTt[,r]==1),r])
  xbar1T<-sapply(1:sn,function(r) mean(datax1T[[r]]))
  var1T<-sapply(1:sn,function(r) var(datax1T[[r]]))
  #mu2T
  datax2T <- lapply(1:sn,function(r) sampleTe[which(sampleTt[,r]==0),r])
  xbar2T<-sapply(1:sn,function(r) mean(datax2T[[r]]))
  var2T<-sapply(1:sn,function(r) var(datax2T[[r]]))
  
  sdT <- sqrt(((sumyT-1)*var1T+(addsample-sumyT-1)*var2T)/(addsample*(addsample-sumyT)))
  #mu1R
  datax1R <- lapply(1:sn,function(r) sampleRe[which(sampleRt[,r]==1),r])
  xbar1R<-sapply(1:sn,function(r) mean(datax1R[[r]]))
  var1R<-sapply(1:sn,function(r) var(datax1R[[r]]))
  #mu2R
  datax2R <- lapply(1:sn,function(r) sampleRe[which(sampleRt[,r]==0),r])
  xbar2R<-sapply(1:sn,function(r) mean(datax2R[[r]]))
  var2R<-sapply(1:sn,function(r) var(datax2R[[r]]))
  
  sdR <- sqrt(((sumyR-1)*var1R+(addsample-sumyR-1)*var2R)/(addsample*(addsample-sumyR)))
  scaleR<- 1/2*((sumyR-1)*var1R+(addsample-sumyR-1)*var2R)
  #sample
  cl <- makeCluster(cores)
  clusterEvalQ(cl,library(invgamma))
  clusterExport(cl, c('pn','sn','x','addsample',
                      'posterior_mu1T','posterior_mu2T','posterior_mu1R','posterior_mu2R',
                      'xbar1T','sdT','xbar2T','xbar1R','sdR','xbar2R','scaleR'))
  dist1T<-parSapply(cl,1:sn,function(r) sample(x,pn,replace = T,prob=posterior_mu1T(x,r)))
  dist2T<-parSapply(cl,1:sn,function(r) sample(x,pn,replace = T,prob=posterior_mu2T(x,r)))
  dist1R<-parSapply(cl,1:sn,function(r) sample(x,pn,replace = T,prob=posterior_mu1R(x,r)))
  dist2R<-parSapply(cl,1:sn,function(r) sample(x,pn,replace = T,prob=posterior_mu2R(x,r)))
  distsigmaR2<-parSapply(cl,1:sn,function(r) rinvgamma(n=pn,shape=(addsample/2),rate=scaleR[r]))
  stopCluster(cl)
  
  cmud<-lapply(1:sn, function(r){
    d<-(distpT[pdbios[[r]],r]*dist1T[pdbios[[r]],r] + dist2T[pdbios[[r]],r] - dist2T[pdbios[[r]],r]*distpT[pdbios[[r]],r]) - 
      (distpR[pdbios[[r]],r]*dist1R[pdbios[[r]],r] + dist2R[pdbios[[r]],r] - dist2R[pdbios[[r]],r]*distpR[pdbios[[r]],r])
    tauR2<-distsigmaR2[pdbios[[r]],r]+distpR[pdbios[[r]],r]*(1-distpR[pdbios[[r]],r])*(dist1R[pdbios[[r]],r]-dist2R[pdbios[[r]],r])^2
    d_scale<-d/sqrt(tauR2)
    return(d_scale)
  })
  mu_bios<-sapply(1:sn, function(r) length(which(abs(cmud[[r]]) < 0.4))/length(pdbios[[r]]))
  mup_bios<-p_bios*mu_bios
  write.table(mup_bios,file = paste("deltamu",deltamu,"s",t,".txt",sep = ''))
}
