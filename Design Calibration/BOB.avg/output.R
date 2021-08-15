################################################################################
# BOB: Bayesian Optimal Design for Biosimilar Trials with Co-Primary Endpoints #
#                                                                              #
# -Numerical Studies                                                           #
#  --Design Calibration                                                        #
#    ---BOB.avg                                                                #
################################################################################

###===========================Parameter Searching============================###
###Grid Searching
lambda<-seq(0.894,0.901,by=0.0002)
len1<-length(lambda)
gamma<-seq(1,1.2,by=0.02)
len2<-length(gamma)
for(c in 1:5){
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
        if(c==1){
          deltamu=-0.223
          mup_bios0<-read.table(paste("deltamu",deltamu,"s",t,".txt",sep = ''))
        }else if(c==2){
          deltamu=0.223
          mup_bios0<-read.table(paste("deltamu",deltamu,"s",t,".txt",sep = ''))
        }else if(c==3){
          #power
          mup_bios0<-read.table(paste("powers",t,".txt",sep = ''))
        }else if(c==4){
          deltap=-0.15
          mup_bios0<-read.table(paste("deltap",deltap,"s",t,".txt",sep = ''))
        }else if(c==5){
          deltap=0.15
          mup_bios0<-read.table(paste("deltap",deltap,"s",t,".txt",sep = ''))
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
    TIE_el<-totals
    muEN_l<-EN
  }else if(c==2){
    TIE_er<-totals
    muEN_r<-EN
  }else if(c==3){
    EP<-totals
  }else if(c==4){
    TIE_sl<-totals
    pEN_l<-EN
  }else if(c==5){
    TIE_sr<-totals
    pEN_r<-EN
  }
}


###Calibration
TIE<-pmax(TIE_el,TIE_er,TIE_sl,TIE_sr)
EN<-pmax(muEN_l,muEN_r,pEN_l,pEN_r)
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
cat("(",lambda[lambda_n],",",gamma[gamma_n],")",sep = '')









