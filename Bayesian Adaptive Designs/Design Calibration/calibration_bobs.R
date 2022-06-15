################################################################################
# BOB: Bayesian Optimal Design for Biosimilar Trials with Co-Primary Endpoints #
#                                                                              #
# -Numerical Studies                                                           #
#  |-Bayesian Adaptive Designs                                                 #
#    |-Design Calibration                                                      #
#      |-BOB.s                                                                 #
################################################################################

###===========================Simulation Settings============================###
maxnsample=160#max sample size
Tmax=4 #stages
nsample=maxnsample/Tmax#sample size per stage
tau2=0.8^2
rho=0
sn=20000
pR=0.5
pT=0.5
overallmuT<- 0.32
overallmuR<- 0

solvemu<-function(p,x,tau2,rho){
  A=matrix(c(p,(1-p),1,-1),nrow = 2,byrow = T)
  b=c(x,sqrt(tau2)*rho/sqrt(p*(1-p)))
  mu<-solve(A,b)
  return(mu)
}
muR<-solvemu(pR,overallmuR,tau2,rho)
sigmaR<-sqrt(tau2-pR*(1-pR)*(muR[1]-muR[2])**2)

###===============================Function================================###
##asymptotic cdf of muT-muR
Pr_dmu<-function(x){
  mean=xbarT - xbarR
  sd=sqrt((nR-1)*(varT+varR)/nR^2)
  r<-pnorm(q = x,mean = mean,sd = sd)
  return(r)
}
##posterior of p
Prpd<-function(x){
  b1<-x-0.20
  b2<-x+0.20
  r1<-dbeta(x,(sumyR[i]+1),(addsample-sumyR[i]+1))
  r2<-pbeta(b2,(sumyT[i]+1),(addsample-sumyT[i]+1)) - pbeta(b1,(sumyT[i]+1),(addsample-sumyT[i]+1))
  r<-r1*r2
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
mup_bios<-vector(mode="numeric",length=sn)
### Fix muT at +/-0.32
for(i in 1:41){
  pT=0.3+0.01*(i-1)
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
  #p
  sumyT<-colSums(sampleTt)
  sumyR<-colSums(sampleRt)
  p_bios<-vector(mode="numeric",length=sn)
  ##mu
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
    p_bios[i]<-integrate(Prpd,0,1)$value
  }
  mup_bios<-mu_bios*p_bios
  write.table(mup_bios,file = paste("muT",overallmuT,"pT",pT,"s",t,".txt",sep = ''))
  }
}
### Fix pT at 0.30 or 0.70
for(i in 1:41){
  overallmuT<- -0.32+0.02*(i-1)
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
  #p
  sumyT<-colSums(sampleTt)
  sumyR<-colSums(sampleRt)
  p_bios<-vector(mode="numeric",length=sn)
  ##mu
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
    p_bios[i]<-integrate(Prpd,0,1)$value
  }
  mup_bios<-mu_bios*p_bios
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
  addsample=nsample*t
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
  #p
  sumyT<-colSums(sampleTt)
  sumyR<-colSums(sampleRt)
  p_bios<-vector(mode="numeric",length=sn)
  ##mu
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
    p_bios[i]<-integrate(Prpd,0,1)$value
  }
  mup_bios<-mu_bios*p_bios
  write.table(mup_bios,file = paste("muT",overallmuT,"pT",pT,"s",t,".txt",sep = ''))
}


###===========================Parameter Searching============================###
###Find maximum type I error
#Fix muT at +/-0.32
power<-matrix(NA,41,4)
for(m in 1:41){
  for(n in 1:2){
    overallmuT=c(-0.32,0.32)[n]
    pT=0.3+0.01*(m-1)
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
colnames(power)<-c("TIE(muT=-0.32)","EN","TIE(muT=0.32)","EN")
power_muT<-power
max_TIE1<-max(power_muT[,c(1,3)])

#Fix pT at 0.30 or 0.70
power<-matrix(NA,33,4)
for(m in 1:33){
  for(n in 1:2){
    pT=c(0.3,0.7)[n]
    overallmuT=-0.32+0.02*(m-1)
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
colnames(power)<-c("TIE(pT=0.3)","EN","TIE(pT=0.7)","EN")
power_pT<-power
max_TIE2<-max(power_pT[,c(1,3)])

if(max_TIE1>=max_TIE2){
  muT_c<-c(-0.32,0.32)[ceiling(which.max(power_muT[,c(1,3)])/41)]
  pT_c<-matrix(rep(seq(0.3,0.7,by = 0.01),4),41,4)[which.max(power_muT[,c(1,3)])]
}else{
  muT_c<-matrix(rep(seq(-0.32,0.32,by = 0.02),4),33,4)[which.max(power_pT[,c(1,3)])]
  pT_c<-c(0.3,0.7)[ceiling(which.max(power_pT[,c(1,3)])/33)]
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
