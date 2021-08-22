################################################################################
# BOB: Bayesian Optimal Design for Biosimilar Trials with Co-Primary Endpoints #
#                                                                              #
# -Numerical Studies                                                           #
#  |-Bayesian Adaptive Designs                                                 #
#    |-Design Calibration                                                      #
#      |-BOB.s                                                                 #
################################################################################

###===========================Simulation Settings============================###
library(mvtnorm)
maxnsample=300#max sample size
Tmax=4 #stages
nsample=maxnsample/Tmax#sample size per stage
tau2=0.8^2
rho=0
sn=20000
pR=0.5
pT=0.5
overallmuT<- 0.223
overallmuR<- 0

solvemu<-function(p,x,tau2,rho){
  A=matrix(c(p,(1-p),1,-1),nrow = 2,byrow = T)
  b=c(x,sqrt(tau2)*rho/sqrt(p*(1-p)))
  mu<-solve(A,b)
  return(mu)
}

muR<-solvemu(pR,overallmuR,tau2,rho)
sigmaR<-sqrt(tau2-pR*(1-pR)*(muR[1]-muR[2])**2)

###========================Simulation Implementation=========================###
mup_bios<-vector(mode="numeric",length=sn)
### Fix muT at +/-0.223
for(i in 1:31){
  pT=0.35+0.01*(i-1)
  muT<-solvemu(pT,overallmuT,tau2,rho)
  sigmaT<-sqrt(tau2-pT*(1-pT)*(muT[1]-muT[2])**2)
  for(t in 1:Tmax){
    addsample=nsample*t#added up sample size at the current stage
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
      sampleRt<-rbind(sampleRt,sampleRt2)
      sampleTe<-rbind(sampleTe,sampleTe2)
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
    write.table(mup_bios,file = paste("muT",overallmuT,"pT",pT,"s",t,".txt",sep = ''))
  }
}
### Fix pT at 0.35 or 0.65
for(i in 1:45){
  overallmuT<- -0.22+0.01*(i-1)
  muT<-solvemu(pT,overallmuT,tau2,rho)
  sigmaT<-sqrt(tau2-pT*(1-pT)*(muT[1]-muT[2])**2)
  for(t in 1:Tmax){
    addsample=nsample*t#added up sample size at the current stage
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
      sampleRt<-rbind(sampleRt,sampleRt2)
      sampleTe<-rbind(sampleTe,sampleTe2)
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
    write.table(mup_bios,file = paste("muT",overallmuT,"pT",pT,"s",t,".txt",sep = ''))
  }
}

###Simulate power
mup_bios<-vector(mode="numeric",length=sn)
pT=0.5
overallmuT=0
muT<-solvemu(pT,overallmuT,tau2,rho)
sigmaT<-sqrt(tau2-pT*(1-pT)*(muT[1]-muT[2])**2)
for(t in 1:Tmax){
  addsample=nsample*t#added up sample size at the current stage
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
    sampleRt<-rbind(sampleRt,sampleRt2)
    sampleTe<-rbind(sampleTe,sampleTe2)
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
  write.table(mup_bios,file = paste("muT",overallmuT,"pT",pT,"s",t,".txt",sep = ''))
}

###===========================Parameter Searching============================###
###Find maximum type I error
#Fix muT at +/-0.223
power<-matrix(NA,31,4)
for(m in 1:31){
  for(n in 1:2){
    overallmuT=c(-0.223,0.223)[n]
    pT=0.35+0.01*(m-1)
    muT<-solvemu(pT,overallmuT,tau2,rho)
    sigmaT<-sqrt(tau2-pT*(1-pT)*(muT[1]-muT[2])**2)
    nextsn<-sn
    EN0<-0
    p4<-1
    for(t in 1:Tmax){
      addsample=(maxnsample/Tmax)*t#added up sample size at the current stage
      Cf<-lambda*(addsample/maxnsample)**gamma
      if(Tmax==1){
        mup_bios0<-read.table(paste("muT",overallmuT,"pT",pT,"s4.txt",sep = ''))
      }else{
        mup_bios0<-read.table(paste("muT",overallmuT,"pT",pT,"s",t,".txt",sep = ''))
      }
      if(t!=1){
        mup_bios0[addf,1]<- -1
      }
      mup_bios<-mup_bios0[,1]
      stagef<-which(mup_bios < (Cf-0.000001) & mup_bios != -1)
      mup_bios[stagef]<- -1
      addf<-which(mup_bios == -1)
      if(t!=Tmax){
        EN0<-EN0+(length(stagef)/sn)*addsample
        p4<- p4-length(stagef)/sn
      }else if(t==Tmax){
        EN0<-EN0+p4*addsample
        power[m,2*n]<-EN0
      }
      nextstagetrials<-setdiff(c(1:sn),addf)
      nextsn<-length(nextstagetrials)
      if(t==Tmax|nextsn==0){
        power[m,2*n-1]<-100*nextsn/sn
      }
    }
  }
}
colnames(power)<-c("TIE(pT=0.35)","EN","TIE(pT=0.65)","EN")
power_muT<-power
max_TIE1<-max(power_muT[,c(1,3)])

#Fix pT at 0.35 or 0.65
power<-matrix(NA,45,4)
for(m in 1:45){
  for(n in 1:2){
    pT=c(0.35,0.65)[n]
    overallmuT=-0.22+0.01*(m-1)
    muT<-solvemu(pT,overallmuT,tau2,rho)
    sigmaT<-sqrt(tau2-pT*(1-pT)*(muT[1]-muT[2])**2)
    nextsn<-sn
    EN0<-0
    p4<-1
    for(t in 1:Tmax){
      addsample=(maxnsample/Tmax)*t#added up sample size at the current stage
      Cf<-lambda*(addsample/maxnsample)**gamma
      if(Tmax==1){
        mup_bios0<-read.table(paste("muT",overallmuT,"pT",pT,"s4.txt",sep = ''))
      }else{
        mup_bios0<-read.table(paste("muT",overallmuT,"pT",pT,"s",t,".txt",sep = ''))
      }
      if(t!=1){
        mup_bios0[addf,1]<- -1
      }
      mup_bios<-mup_bios0[,1]
      stagef<-which(mup_bios < (Cf-0.000001) & mup_bios != -1)
      mup_bios[stagef]<- -1
      addf<-which(mup_bios == -1)
      if(t!=Tmax){
        EN0<-EN0+(length(stagef)/sn)*addsample
        p4<- p4-length(stagef)/sn
      }else if(t==Tmax){
        EN0<-EN0+p4*addsample
        power[m,2*n]<-EN0
      }
      nextstagetrials<-setdiff(c(1:sn),addf)
      nextsn<-length(nextstagetrials)
      if(t==Tmax|nextsn==0){
        power[m,2*n-1]<-100*nextsn/sn
      }
    }
  }
}
colnames(power)<-c("TIE(pT=0.35)","EN","TIE(pT=0.65)","EN")
power_pT<-power
max_TIE2<-max(power_pT[,c(1,3)])

if(max_TIE1>=max_TIE2){
  muT_c<-c(-0.223,0.223)[ceiling(which.max(power_muT[,c(1,3)])/31)]
  pT_c<-matrix(rep(seq(0.35,0.65,by = 0.01),4),31,4)[which.max(power_muT[,c(1,3)])]
}else{
  muT_c<-matrix(rep(seq(-0.22,0.22,by = 0.01),4),45,4)[which.max(power_pT[,c(1,3)])]
  pT_c<-c(0.35,0.65)[ceiling(which.max(power_pT[,c(1,3)])/45)]
}

###Grid Searching
lambda<-seq(0.939,0.944,by=0.0002)
len1<-length(lambda)
gamma<-seq(1,1.2,by=0.02)
len2<-length(gamma)
for(c in 1:2){
  totals<-matrix(NA,len1,len2)
  EN<-matrix(NA,len1,len2)
  for(m in 1:len1){
    for(n in 1:len2){
      nextsn<-sn
      EN0<-0
      p4<-1
      for(t in 1:Tmax){
        addsample=(maxnsample/Tmax)*t#added up sample size at the current stage
        Cf<-lambda[m]*(addsample/maxnsample)**gamma[n]
        if(c==1){mup_bios0<-read.table(paste("muT",muT_c,"pT",pT_c,"s",t,".txt",sep = ''))}
        else if(c==2){mup_bios0<-read.table(paste("muT0pT0.5s",t,".txt",sep = ''))}
        if(t!=1){
          mup_bios0[addf,1]<- -1
        }
        mup_bios<-mup_bios0[,1]
        stagef<-which(mup_bios < (Cf-0.000001) & mup_bios != -1)
        mup_bios[stagef]<- -1
        addf<-which(mup_bios == -1)
        if(t!=Tmax){
          EN0<-EN0+(length(stagef)/sn)*addsample
          p4<- p4-length(stagef)/sn
        }else if(t==Tmax){
          EN0<-EN0+p4*addsample
          EN[m,n]<-EN0
        }
        nextstagetrials<-setdiff(c(1:sn),addf)
        nextsn<-length(nextstagetrials)
        if(t==Tmax|nextsn==0){
          totals[m,n]<-nextsn/sn
        }
      }
    }
  }
  if(c==1){
    TIE_s<-totals
    EN_s<-EN
  }else if(c==2){
    EP<-totals
  }
}

###Calibration
alpha_s<-which(TIE_s<=0.05) #alpha_mu
bEP<-max(EP[alpha_s])-0.01
power<-which(EP>=bEP)
r<-intersect(alpha_s,power)
power_n=floor(length(r)/2)
filter_r<-order(EP[r],decreasing = T)[1:power_n]
finalset<-r[filter_r]
value<-finalset[which.min(EN_s[finalset])]
### Output
lambda_n<-value%%length(lambda)
gamma_n<-ceiling(value/length(lambda))
cat("Optimal parameters for BOBs:(",lambda[lambda_n],",",gamma[gamma_n],")",sep = '')
