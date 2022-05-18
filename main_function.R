# rm(list=ls())
#library(spef)
library(tidyverse)
library(survival)
#data(skinTumor)
library(nleqslv)
#library(openxlsx)
library(extraDistr)
#library(fda)
library(splines)
library(splines2)
#library(ggpubr)
#library(matrixcalc)
library(doParallel)

if (exists("cl")) rm(cl)
cl <- makeCluster(40)

#cl2 <- doAzureParallel::makeCluster("simulation/cluster.json")
registerDoParallel(cl)
getDoParWorkers()

nrep=1000

p1 <- 2
p2 <- 1

gamma1<- 0.5
gamma2 <- 0.5

tau<-1
nknots=3
ndegree=3

###event time
f.mu1<-function(t){
  2*t
  #(sin(4*pi*t)+4*pi*t)/(2*pi)
}

f.mu2<-function(t){
  2*t^2
  #(cos(4*pi*t)+4*pi*t)/(2*pi)
}

##observation time
f.mutilde1<-function(t){
  (3*t+4)/2
}

f.mutilde2<-function(t){
  (2*t+3)/2
}

######################setting 1
f.beta1<-function(t,a){
  #t^2/2 #setting2
  t #setting1
  #(sin(4*pi*t)+4*pi*t)/24
}

###var1
data.gen <-
  function(id,                                #---just for 1 subject
           tau                              
  ){
    covar <- generate.covariate(tau)
    #B11,B12,V1,B21,B22,V2,B31,B32,V3,B41,B42,V4
    B11 <- covar[1]
    B12 <- covar[2]
    V1 <- covar[3]
    B21 <- covar[4]
    B22 <- covar[5]
    V2 <- covar[6]
    B31 <- covar[7]
    B32 <- covar[8]
    V3 <- covar[9]
    # B41 <- covar[10]
    # B42 <- covar[11]
    # V4 <- covar[12]
    # ##Censor
    Censor <-
      runif(1, tau * 0.5, tau * 0.99) ###one for each subject, to ensure censoring rate is reasonable
    
    P <- 1
    #rgamma(1,1,1)
    #---To generate NHPP
    ##for the 1st type event
    obs.bar1 <-
      optimize(function(x) P*f.mutilde1(x),
               interval = c(0, Censor),
               maximum = TRUE
      )$objective
    obs.Time1 <- numeric()
    
    
    while(length(obs.Time1)==0){
      obs.Time1 <-
        NHPP.gen(
          lambdat = function(x) P*f.mutilde1(x),
          tau = Censor,
          Censor,
          lambda.bar = obs.bar1
        )
      #browser()
      obs.Time1 <- unique(floor(obs.Time1*100)/100)
      obs.Time1 <- obs.Time1[obs.Time1>0]
    }
    nn1 <- length(obs.Time1)
    
    ##for the 2nd type event
    obs.bar2 <-
      optimize(function(x) P*f.mutilde2(x),
               interval = c(0, Censor),
               maximum = TRUE
      )$objective
    obs.Time2 <- numeric()
    
    while(length(obs.Time2)==0){
      obs.Time2 <-
        NHPP.gen(
          lambdat=function(x) P*f.mutilde2(x),
          tau = Censor,
          Censor,
          lambda.bar = obs.bar2
        )
      #browser()
      obs.Time2 <- unique(floor(obs.Time2*100)/100)
      obs.Time2 <- obs.Time2[obs.Time2>0]
    }
    nn2 <- length(obs.Time2)
    
    #browser()
    ####generent count
    Q <- 1
    #Q <- rgamma(1,1,1)
    #rgamma(1,shape=0.5,scale=2)
    #rgamma(1,1,1)
    #### Ni1
    mean_val1 <- Q*sapply(obs.Time1, function(tt) event.mean1(t = tt, tau, N, B11,B12,V1,B21,B22,V2
                                                              ,B31,B32,V3
                                                              #,B41,B42,V4
    ))
    mean_val11 <- Q*sapply(tau/2, function(tt) event.mean1(t = tt, tau, N, B11,B12,V1,B21,B22,V2
                                                           ,B31,B32,V3
                                                           #,B41,B42,V4
    ))
    mean_val12 <- Q*sapply(tau, function(tt) event.mean1(t = tt, tau, N, B11,B12,V1,B21,B22,V2
                                                         ,B31,B32,V3
                                                         #,B41,B42,V4
    ))
    
    mean_diff1 <- diff(c(0,mean_val1))
    if(sum(mean_diff1<0)) browser()
    dcount1 <- rpois(nn1,mean_diff1)
    if(is.na(sum(dcount1))) browser()
    count1 <- cumsum(dcount1)
    countprocess1 <- stepfun(obs.Time1,c(0,count1))
    atime1 <- c(obs.Time1, Censor, V1,V2
                ,V3
                #,V4
    )
    eventind1 <- c(rep(1, nn1), rep(0,(1+p1+p2)))#
    #### Ni2
    mean_val2 <- Q*sapply(obs.Time2, function(tt) event.mean2(t = tt, tau, N, B11,B12,V1,B21,B22,V2
                                                              ,B31,B32,V3
                                                              #,B41,B42,V4
    ))
    mean_val21 <- Q*sapply(tau/2, function(tt) event.mean2(t = tt, tau, N, B11,B12,V1,B21,B22,V2
                                                           ,B31,B32,V3
                                                           #,B41,B42,V4
    ))
    mean_val22 <- Q*sapply(tau, function(tt) event.mean2(t = tt, tau, N, B11,B12,V1,B21,B22,V2
                                                         ,B31,B32,V3
                                                         #,B41,B42,V4
    ))
    
    mean_diff2 <- diff(c(0,mean_val2))
    if(sum(mean_diff2<0)) browser()
    dcount2 <- rpois(nn2,mean_diff2)
    if(is.na(sum(dcount2))) browser()
    count2 <- cumsum(dcount2)
    countprocess2 <- stepfun(obs.Time2,c(0,count2))
    atime2 <- c(obs.Time2, Censor, V1,V2
                ,V3
                #,V4
    )
    eventind2 <- c(rep(1, nn2), rep(0,(1+p1+p2)))#
    #########
    ###to find recurrent event process N(t)
    #browser()
    ####
    temp1 <- cbind(id,
                   count1=countprocess1(atime1), atime1, eventind1)
    temp1 <- temp1[order(temp1[, "atime1"]), ]
    res1 <-
      cbind(temp1[, c("id", "count1")],
            c(0, temp1[, "atime1"][-(nn1+1+p1+p2)]),  ##tstart
            temp1[, c("atime1", "eventind1")], #time, event
            ifelse(temp1[, "atime1"] <= V1, B11, B12),ifelse(temp1[, "atime1"] <= V2, B21, B22)
            ,ifelse(temp1[, "atime1"] <= V3, B31, B32),
            mean_val11,
            mean_val12
            #,ifelse(temp1[, "atime1"] <= V4, B41, B42)
      )
    #####
    temp2 <- cbind(id,
                   count2=countprocess2(atime2), atime2, eventind2)
    temp2 <- temp2[order(temp2[, "atime2"]), ]
    res2 <-
      cbind(temp2[, c("id", "count2")],
            c(0, temp2[, "atime2"][-(nn2+1+p1+p2)]),  ##tstart
            temp2[, c("atime2", "eventind2")], #time, event
            ifelse(temp2[, "atime2"] <= V1, B11, B12),ifelse(temp2[, "atime2"] <= V2, B21, B22)
            ,ifelse(temp2[, "atime2"] <= V3, B31, B32),
            mean_val21,
            mean_val22
            #,ifelse(temp2[, "atime2"] <= V4, B41, B42)
      )
    if(res2[dim(res2)[1],5]==1) browser()
    #colnames(res) <- c("id","count","tstart","time","event",paste("W",1:p1,sep=""),paste("Z",1:p2,sep=""))
    ###combine
    return(list(res1,res2))
  }

nsub=300
simres1 <- foreach(i = 1:nrep,
                   .packages = c("MASS", "extraDistr", "tidyverse","nleqslv","foreach","survival","lubridate","splines","splines2")) %dopar% {
                     #ptbegin <- Sys.time()
                     data.all=data.all.gen(nsub,                                #---Sample size
                                           tau
                                           #,N                                #---Length of study = 6
                     )                                  #how many subintervals of tau
                     #test<-data.all.gen(10,obs.rate,6)
                     #browser()
                     data_tran <- data.tran(data.all,tau)
                     
                     # Add spline
                     ###1 st type
                     data1 <- data_tran[[1]]
                     STN1 <- cbind(data1$stn11,data1$stn12)
                     data1 <- data1 %>%
                       select(-c(stn11,stn12))
                     yy1 <- data1 %>%  dplyr::select(id, count1, tstart1, time1,event1,timegroup1)
                     ww1<-data1 %>%  dplyr::select(W1,W2)
                     zz1 <- data1 %>%  dplyr::select(-id, -count1, -tstart1, -time1,-event1, -timegroup1,-W1,-W2
                                                     #,-W2
                     )
                     #browser()
                     knots1<-unname(quantile(yy1[,"time1"], seq(0, 1, length=nknots+2)))[-c(1,(nknots+2))]
                     basis_val1 <- ns(yy1[,"time1"],knots=knots1,df = ndegree,intercept = T,Boundary.knots = c(0,1))
                     ztilde1 <- do.call(cbind,lapply(array_branch(as.matrix(zz1),2),function(zzz) zzz*basis_val1)) #p2*(nknots+ndegree-1)
                     xx1<-cbind(ww1,ztilde1)
                     colnames(xx1) <- paste("X",1:(dim(xx1)[2]),sep="")
                     data1 <- cbind(yy1,xx1)
                     
                     ###2 nd type
                     data2 <- data_tran[[2]]
                     STN2 <- cbind(data2$stn21,data2$stn22)
                     data2 <- data2 %>%
                       select(-c(stn21,stn22))
                     yy2 <- data2 %>%  dplyr::select(id, count2, tstart2, time2,event2,timegroup2)
                     ww2<-data2 %>%  dplyr::select(W1,W2)
                     zz2 <- data2 %>%  dplyr::select(-id, -count2, -tstart2, -time2,-event2, -timegroup2,-W1,-W2)
                     #browser()
                     knots2<-unname(quantile(yy2[,"time2"], seq(0, 1, length=nknots+2)))[-c(1,(nknots+2))]
                     basis_val2 <- ns(yy2[,"time2"],knots=knots2,df = ndegree,intercept = T,Boundary.knots = c(0,1))
                     ztilde2 <- do.call(cbind,lapply(array_branch(as.matrix(zz2),2),function(zzz) zzz*basis_val2)) #p2*(nknots+ndegree-1)
                     xx2<-cbind(ww2,ztilde2)
                     colnames(xx2) <- paste("X",1:(dim(xx2)[2]),sep="")
                     data2 <- cbind(yy2,xx2)
                     
                     ##estimating equation
                     #browser()
                     alpha_hat<-nleqslv(rep(0.5,p1+p2*(nknots+ndegree-1)),function(ee) panelEE(ee,data1,data2))$x
                     #####
                     #browser()
                     sand.temp<-get.sand(alpha_hat,data1,data2)
                     sand.temp1<-sand.temp[[1]]
                     sand.temp2<-sand.temp[[2]]
                     ###1st type
                     mu1hat <- lapply(1:(range(data1$timegroup1)[2]-1), function(i) sand.temp1[[i]][[1]])
                     S1.0<-lapply(1:(range(data1$timegroup1)[2]-1), function(i) sand.temp1[[i]][[2]])
                     S1.1<-lapply(1:(range(data1$timegroup1)[2]-1), function(i) sand.temp1[[i]][[3]])
                     S1.2<-lapply(1:(range(data1$timegroup1)[2]-1), function(i) sand.temp1[[i]][[4]])
                     S1.10.1<-lapply(1:length(S1.0), function(i) matrix(S1.2[[i]],nrow=(p1+p2*(nknots+ndegree-1)))/as.numeric(S1.0[[i]]))
                     X.bar1<- (lapply(1:length(S1.0), function(i) S1.1[[i]]/as.numeric(S1.0[[i]])))
                     S1.10.2<-lapply(1:length(S1.0), function(s) outer(unlist(X.bar1[[s]]),unlist(X.bar1[[s]]))) #ntimegroup,16*16
                     #browser()
                     X.deriv1<- lapply(1:length(S1.0),function(i) S1.10.1[[i]]-S1.10.2[[i]]) #ntimegroup,16*16,
                     
                     ###2nd type
                     mu2hat <- lapply(1:(range(data2$timegroup2)[2]-1), function(i) sand.temp2[[i]][[1]])
                     S2.0<-lapply(1:(range(data2$timegroup2)[2]-1), function(i) sand.temp2[[i]][[2]])
                     S2.1<-lapply(1:(range(data2$timegroup2)[2]-1), function(i) sand.temp2[[i]][[3]])
                     S2.2<-lapply(1:(range(data2$timegroup2)[2]-1), function(i) sand.temp2[[i]][[4]])
                     S2.10.1<-lapply(1:length(S2.0), function(i) matrix(S2.2[[i]],nrow=(p1+p2*(nknots+ndegree-1)))/as.numeric(S2.0[[i]]))
                     X.bar2<- (lapply(1:length(S2.0), function(i) S2.1[[i]]/as.numeric(S2.0[[i]])))
                     S2.10.2<-lapply(1:length(S2.0), function(s) outer(unlist(X.bar2[[s]]),unlist(X.bar2[[s]]))) #ntimegroup,16*16
                     #browser()
                     X.deriv2<- lapply(1:length(S2.0),function(i) S2.10.1[[i]]-S2.10.2[[i]]) #ntimegroup,16*16,
                     
                     #browser()
                     sanda<-get.sanda(alpha_hat,data1,data2,X.deriv1,X.deriv2)
                     sanda1 <- sanda[[1]]# ignore (-)
                     sanda2 <- sanda[[2]]
                     #browser()
                     #mu0hat <- get.mu0hat(alpha_hat,data.tran)
                     #browser()
                     sandb<-get.sandb(alpha_hat,data1,data2,X.bar1,X.bar2,mu1hat,mu2hat)
                     ###
                     #browser()
                     #if(!is.singular.matrix(sanda1+sanda2)) browser()
                     var.sandwich<- solve((sanda1+sanda2)) %*% (sandb) %*% solve((sanda1+sanda2))
                     BIC <- dim(sandb)[1]*log(nsub)-2*sum(diag(sandb))
                     AIC <- dim(sandb)[1]*2-2*sum(diag(sandb))
                     rescol<-list(alpha_est=alpha_hat,
                                  mu1hat =mu1hat,
                                  mu2hat =mu2hat,
                                  var.sandwich=var.sandwich,
                                  B = sandb,
                                  AIC=AIC,
                                  BIC=BIC,
                                  #bsres=bsres,
                                  mean_val11 = unique(STN1[,1]), #match to sample
                                  mean_val12 = unique(STN1[,2]),
                                  mean_val21 = unique(STN2[,1]), #match to sample
                                  mean_val22 = unique(STN2[,2])
                     )
                   }
