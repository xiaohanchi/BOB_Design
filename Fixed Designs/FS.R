################################################################################
# BOB: Bayesian Optimal Design for Biosimilar Trials with Co-Primary Endpoints #
#                                                                              #
# -Numerical Studies                                                           #
#  |-Fixed Designs                                                             #
#    |-FS                                                                      #
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
  sampleRt<-sapply(1:sn, function(r) rbinom(nsample,1,pR))
  #eva of p
  peT<-apply(sampleTt/addsample,2,sum)
  peR<-apply(sampleRt/addsample,2,sum)
  resultp<-sapply(1:sn, function(r) propt(peT[r],peR[r],deltap=0.20,n=addsample))
  bios_p<-which(resultp==1)
      
  bios<-bios_p
  power<-length(bios)/sn
}
cat("The power (type I error rate) of the design FS:",power,"\n",sep = '')
