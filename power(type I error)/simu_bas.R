################################################################################
# BOB: Bayesian Optimal Design for Biosimilar Trials with Co-Primary Endpoints #
#                                                                              #
# -Numerical Studies                                                           #
#  |-Bayesian Adaptive Designs                                                 #
#    |-Empirical power or type I error rate                                    #
#      |-BAS                                                                   #
################################################################################

###===========================Simulation Settings============================###
maxnsample=400#max sample size
Tmax=4 #stages
nsample=maxnsample/Tmax#sample size per stage
pR=0.5
pT=0.35
overallmuT=0
overallmuR=0
lambda=0.9527
gamma=1.04
sn=20000

###===========================Posterior Function============================###
Prpd<-function(x){
  b1<-x-0.15
  b2<-x+0.15
  r1<-dbeta(x,(sumyR[i]+1),(addsample-sumyR[i]+1))
  r2<-pbeta(b2,(sumyT[i]+1),(addsample-sumyT[i]+1)) - pbeta(b1,(sumyT[i]+1),(addsample-sumyT[i]+1))
  r<-r1*r2
  return(r)
}

###========================Simulation Implementation=========================###
nextsn<-sn
pstop<-c(0,0,0)
for(t in 1:Tmax){
  addsample=nsample*t#added up sample size at the current stage
  Cf<-lambda*(addsample/maxnsample)**gamma
  set.seed(2333+10*t)
  if(t==1){
    sampleTt<-sapply(1:sn, function(r) rbinom(nsample,1,pT))
    sampleRt<-sapply(1:sn, function(r) rbinom(nsample,1,pR))
  }else{
    sampleTt2<-sapply(1:sn, function(r) rbinom(nsample,1,pT))
    sampleRt2<-sapply(1:sn, function(r) rbinom(nsample,1,pR))
    sampleTt<-rbind(sampleTt,sampleTt2)
    sampleRt<-rbind(sampleRt,sampleRt2)
  }
  #posterior of diff p
  sumyT<-colSums(sampleTt)
  sumyR<-colSums(sampleRt)
  p_bios<-vector(mode="numeric",length=sn)
  for(i in 1:sn){
    p_bios[i]<-integrate(Prpd,0,1)$value
  }
  if(t!=1){
    p_bios[addf]<- -1
  }
  stagef<-which(p_bios < (Cf-0.000001) & p_bios != -1)
  p_bios[stagef]<- -1
  addf<-which(p_bios == -1)
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

cat("The power (or type I error rate) of the design BAS:",power,"\n",sep = '')
cat("Expected Sample Size(EN):",EN0,"\n",sep = '')
