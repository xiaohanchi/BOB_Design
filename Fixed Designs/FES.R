################################################################################
# BOB: Bayesian Optimal Design for Biosimilar Trials with Co-Primary Endpoints #
#                                                                              #
# -Numerical Studies                                                           #
#  |-Fixed Designs                                                             #
#    |-FES                                                                     #
################################################################################

###===========================Simulation Settings============================###
maxnsample=200#max sample size
Tmax=1 #stages
nsample=maxnsample/Tmax#sample size per stage
pR=0.5
pT=0.35
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
muR<-solvemu(pR,overallmuR,tau2,rho)
sigmaR<-sqrt(tau2-pR*(1-pR)*(muR[1]-muR[2])**2)
muT<-solvemu(pT,overallmuT,tau2,rho)
sigmaT<-sqrt(tau2-pT*(1-pT)*(muT[1]-muT[2])**2)

propt<-function(p1,p2,deltap,n){
  var<-p1*(1-p1)/n+p2*(1-p2)/n
  t1=(p1-p2+deltap)/sqrt(var)
  t2=(deltap-(p1-p2))/sqrt(var)
  if(t1 >= qt(p = 0.95,df = n-2) && t2 >= qt(p = 0.95,df = n-2)){
    return(1)
  }else{return(0)}
}

###========================Simulation Implementation=========================###
for(t in 1:Tmax){
  addsample=nsample*t#added up sample size at the current stage
  set.seed(233+10*t)
  sampleTt<-sapply(1:sn, function(r) rbinom(nsample,1,pT))
  sampleTe<-sapply(1:sn, function(r) rnorm(nsample,(sampleTt[,r]*muT[1]+(1-sampleTt[,r])*muT[2]),sigmaT))
  sampleRt<-sapply(1:sn, function(r) rbinom(nsample,1,pR))
  sampleRe<-sapply(1:sn, function(r) rnorm(nsample,(sampleRt[,r]*muR[1]+(1-sampleRt[,r])*muR[2]),sigmaR))
  #eva of p
  peT<-apply(sampleTt/addsample,2,sum)
  peR<-apply(sampleRt/addsample,2,sum)
  resultp<-sapply(1:sn, function(r) propt(peT[r],peR[r],deltap=0.15,n=addsample))
  bios_p<-which(resultp==1)
  #mu
  datamuT <- lapply(1:sn,function(r) sampleTe[,r])
  datamuR <- lapply(1:sn,function(r) sampleRe[,r])
  pmu1<-sapply(1:sn,function(r) t.test(datamuT[[r]],datamuR[[r]],alternative = "greater",mu = -0.223)$`p.value`)
  pmu2<-sapply(1:sn,function(r) t.test(datamuT[[r]],datamuR[[r]],alternative = "less",mu = 0.223)$`p.value`)
  bios_mu<-intersect(which(pmu1<0.05),which(pmu2<0.05))

  bios<-intersect(bios_p,bios_mu)
  power<-length(bios)/sn
}
cat("The power (type I error rate) of the design FES:",power,"\n",sep = '')