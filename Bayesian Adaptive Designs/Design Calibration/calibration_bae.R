################################################################################
# BOB: Bayesian Optimal Design for Biosimilar Trials with Co-Primary Endpoints #
#                                                                              #
# -Numerical Studies                                                           #
#  |-Bayesian Adaptive Designs                                                 #
#    |-Design Calibration                                                      #
#      |-BAE                                                                   #
################################################################################

###===========================Simulation Settings============================###
maxnsample=160#max sample size
Tmax=4 #stages
nsample=maxnsample/Tmax#sample size per stage
pT=0.5
pR=0.5
overallmuR<-0
tau2=0.8^2
rho=0
pn=20000
sn=20000#simulation times

solvemu<-function(p,x,tau2,rho){
  A=matrix(c(p,(1-p),1,-1),nrow = 2,byrow = T)
  b=c(x,sqrt(tau2)*rho/sqrt(p*(1-p)))
  mu<-solve(A,b)
  return(mu)
}
muR<-solvemu(pR,overallmuR,tau2,rho)
sigmaR<-sqrt(tau2-pR*(1-pR)*(muR[1]-muR[2])**2)

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
for(j in 1:3){
  overallmuT<- c(-0.32,0.32,0)[j]
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
    write.table(mu_bios,file = paste("n",maxnsample,"rho",rho,"muT",overallmuT,"s",t,"mu_bios.txt",sep = ''))
  }
}

###===========================Parameter Searching============================###
###Grid Searching
lambda=seq(0.948,0.951,by = 0.0001)
len1<-length(lambda)#m
gamma<-seq(1,1.2,by = 0.02)
len2<-length(gamma)#n
totals<-EN<-matrix(NA,len1,len2)
for(c in 1:3){
  for(n in 1:len2){
    for(m in 1:len1){
      nextsn<-sn
      EN0<-0
      p4<-1
      for(t in 1:Tmax){
        addsample=nsample*t
        Cf<-lambda[m]*(addsample/maxnsample)**gamma[n]
        if(c==1){
          overallmuT=-0.32
          mu_bios0<-read.table(paste("n",maxnsample,"rho",rho,"muT",overallmuT,"s",t,"mu_bios.txt",sep = ''))
        }else if(c==2){
          overallmuT=0.32
          mu_bios0<-read.table(paste("n",maxnsample,"rho",rho,"muT",overallmuT,"s",t,"mu_bios.txt",sep = ''))
        }else if(c==3){
          overallmuT=0
          mu_bios0<-read.table(paste("n",maxnsample,"rho",rho,"muT",overallmuT,"s",t,"mu_bios.txt",sep = ''))
        }
        if(t!=1){
          mu_bios0[addf,1]<- -1
        }
        mu_bios<-mu_bios0[,1]
        stagef<-which(mu_bios < (Cf-0.000001) & mu_bios != -1)
        mu_bios[stagef]<- -1
        addf<-which(mu_bios == -1)
        if(t!=Tmax){
          EN0<-EN0+(length(stagef)/sn)*addsample
          p4<- p4-length(stagef)/sn
        }else if(t==Tmax){
          EN0<- EN0+p4*addsample
          EN[m,n]<- EN0
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
    errorl<-totals
    muEN_l<-EN
  }else if(c==2){
    errorr<-totals
    muEN_r<-EN
  }else if(c==3){
    EP<-totals
  }
}

###Calibration
TIE<-pmax(errorl,errorr)
EN<-pmax(muEN_l,muEN_r)
alpha<-which(TIE<=0.05) #alpha_mu
bEP<-max(EP[alpha])-0.01
power<-which(EP>=bEP)
r<-intersect(alpha,power)
power_n=floor(length(r)/2)
filter_r<-order(EP[r],decreasing = T)[1:power_n]
finalset<-r[filter_r]
value<-finalset[which.min(EN[finalset])]

### Output
lambda_n<-value%%length(lambda)
gamma_n<-ceiling(value/length(lambda))
cat("Optimal parameters for BAE:(",lambda[lambda_n],",",gamma[gamma_n],")",sep = '')
