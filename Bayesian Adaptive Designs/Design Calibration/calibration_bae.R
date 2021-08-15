################################################################################
# BOB: Bayesian Optimal Design for Biosimilar Trials with Co-Primary Endpoints #
#                                                                              #
# -Numerical Studies                                                           #
#  |-Bayesian Adaptive Designs                                                 #
#    |-Design Calibration                                                      #
#      |-BAE                                                                   #
################################################################################

###===========================Simulation Settings============================###
maxnsample=400#max sample size
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
posterior_muR<-function(mu){
  xstd<-(mu-xbarR)/sdR
  r<-dt(xstd,nR)/sdR
  return(r)
}

Prmud<-function(x){
  b1 <- (x-0.223-xbarT)/sdT
  b2 <- (x+0.223-xbarT)/sdT
  r1 <- posterior_muR(x)
  r2 <- pt(b2,nT)-pt(b1,nT)
  r <- r1*r2
  return(r)
}

###========================Simulation Implementation=========================###
for(i in 1:3){
  overallmuT<- c(-0.223,0.223,0)[i]
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
      sdT<-sqrt((nT-1)*var(dataxT)/nT^2)
      dataxR <- sampleRe[,i]
      nR<-length(dataxR)
      xbarR<-mean(dataxR)
      sdR<-sqrt((nR-1)*var(dataxR)/nR^2)
      mu_bios[i]<-integrate(Prmud,-5,5)$value
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
          overallmuT=-0.223
          mu_bios0<-read.table(paste("n",maxnsample,"rho",rho,"muT",overallmuT,"s",t,"mu_bios.txt",sep = ''))
        }else if(c==2){
          overallmuT=0.223
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
alpha<-which(TIE<=0.05) 
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
cat("(",lambda[lambda_n],",",gamma[gamma_n],")",sep = '')
