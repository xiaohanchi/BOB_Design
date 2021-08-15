maxnsample=400#max sample size
Tmax=4 #stages
nsample=maxnsample/Tmax#sample size per stage
pR=0.5
sn=20000
setwd("D:/Master/summer/0619simu/results/BFS")

Prpd<-function(x){
  b1<-x-0.15
  b2<-x+0.15
  r1<-dbeta(x,(sumyR[i]+1),(addsample-sumyR[i]+1))
  r2<-pbeta(b2,(sumyT[i]+1),(addsample-sumyT[i]+1)) - pbeta(b1,(sumyT[i]+1),(addsample-sumyT[i]+1))
  r<-r1*r2
  return(r)
}

###save
for(i in 1:3){
  pT=c(0.35,0.65,0.5)[i]
  for(t in 1:Tmax){
    nsample=maxnsample/Tmax
    addsample=(maxnsample/Tmax)*t
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
    sumyT<-colSums(sampleTt)
    sumyR<-colSums(sampleRt)
    p_bios<-vector(mode="numeric",length=sn)
    for(i in 1:sn){
      p_bios[i]<-integrate(Prpd,0,1)$value
    }
    write.table(p_bios,file = paste("n",maxnsample,"pT",pT,"s",t,"p_bios.txt",sep = ''))
  }
}

###search
lambda=seq(0.952,0.953,by = 0.0001)
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
          pT=0.35
          p_bios0<-read.table(paste("n",maxnsample,"pT",pT,"s",t,"p_bios.txt",sep = ''))
        }else if(c==2){
          pT=0.65
          p_bios0<-read.table(paste("n",maxnsample,"pT",pT,"s",t,"p_bios.txt",sep = ''))
        }else if(c==3){
          pT=0.5
          p_bios0<-read.table(paste("n",maxnsample,"pT",pT,"s",t,"p_bios.txt",sep = ''))
        }
        if(t!=1){
          p_bios0[addf,1]<- -1
        }
        p_bios<-p_bios0[,1]
        stagef<-which(p_bios < (Cf-0.000001) & p_bios != -1)
        p_bios[stagef]<- -1
        addf<-which(p_bios == -1)
        if(t!=Tmax){
          EN0<-EN0+(length(stagef)/nextsn)*addsample
          p4<- p4-length(stagef)/nextsn
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

lambda_n<-value%%length(lambda)
gamma_n<-ceiling(value/length(lambda))
cat("(",lambda[lambda_n],",",gamma[gamma_n],")",sep = '')


