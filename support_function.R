#######covariate
generate.covariate<-function(tau){
  ### for W1
  B11 <- runif(1,0,0.5)
  B12 <- runif(1,0.5,1)
  V1 <- runif(1,0,tau)
  ### for W2
  B21 <- runif(1,0,0.5)
  B22 <- runif(1,0.5,1)
  V2 <- runif(1,0,tau)
  ### for Z1
  B31 <- runif(1,0,0.5)
  B32 <- runif(1,0.5,1)
  V3 <- runif(1,0,tau)
  # ### for Z2
  # B41 <- runif(1,0,0.5)
  # B42 <- runif(1,0.5,1)
  # V4 <- runif(1,0,tau)
  return(c(B11,B12,V1,
           B21,B22,V2,
           B31,B32,V3
           #B41,B42,V4
  ))
}

event.mean1<-function(t,tau,N,B11,B12,V1,B21,B22,V2,B31,B32,V3,B41,B42,V4){
  f.mu1(t)*exp(gamma1*ifelse(t>V1,B12,B11) 
               + gamma2*ifelse(t>V2,B22,B21) 
               + f.beta1(t)*ifelse(t>V3,B32,B31)
               #+ f.beta2(t)*ifelse(t>V4,B42,B41)
  )
}

event.mean2<-function(t,tau,N,B11,B12,V1,B21,B22,V2,B31,B32,V3,B41,B42,V4){
  f.mu2(t)*exp(gamma1*ifelse(t>V1,B12,B11) 
               + gamma2*ifelse(t>V2,B22,B21) 
               + f.beta1(t)*ifelse(t>V3,B32,B31) 
               #+ f.beta2(t)*ifelse(t>V4,B42,B41)
  )
}

NHPP.gen <- function(lambdat, tau, Censor, lambda.bar) {
  m <- 0
  while (m == 0) {
    nn <- rtpois(1, Censor * lambda.bar, 0)
    tt <- sort(runif(nn, min = 0, max = 1)*Censor)
    ind <- (sapply(tt, lambdat) /  lambda.bar) > runif(nn)
    tt <- tt[ind]
    m <- length(tt)
  }
  return(tt)
}

event.process<-function(obs.Time,event.time){
  count <- sapply(obs.Time,function(ot) sum(event.time <= ot))
  return(count)
}

##generate data for 1 subject!!!!
#data.gen(id=1,tau=1)

##################
# variable columns: id, cumcount, start time, observed time, event indcator, all other covariates 
#
#
data.all.gen<- function(nsub,                                #---Sample size
                        tau                               #---Length of study = 6
)                                  #how many subintervals of tau                    
{
  res <- lapply(1:nsub,function(id) data.gen(id=id,tau)) 
  #browser()
  res1 <- list()
  res1 <- lapply(1:nsub, function(i) res[[i]][[1]])
  res1 <- do.call(rbind,res1)
  #p1 <- 1
  #p2 <- 2
  colnames(res1) <- c("id","count1","tstart1","time1","event1",paste("W",1:p1,sep=""),paste("Z",1:p2,sep=""),"stn11","stn12")
  #
  res2 <- list()
  res2 <- lapply(1:nsub, function(i) res[[i]][[2]])
  res2 <- do.call(rbind,res2)
  colnames(res2) <- c("id","count2","tstart2","time2","event2",paste("W",1:p1,sep=""),paste("Z",1:p2,sep=""),"stn21","stn22")
  return(list(as.data.frame(res1),as.data.frame(res2)))
}  

#set.seed(2)
#data1 <- data.all.gen(nsub=10,tau=1)

data.tran<-function(data, tau){
  u <- tau
  # c<-0
  
  ####1st type
  data1 <- data[[1]]
  time_grid1 <- sort(unique(data1[data1[,"event1"]==1,][,"time1"]))
  #transformed_data <- as.data.frame(transformed_data)
  transformed_data1 <-
    survSplit(
      Surv(tstart1, time1, event1) ~ .,
      data1,
      cut = time_grid1,
      episode = "timegroup1"
    )
  ###tocheck the last time group
  check1 <- transformed_data1[transformed_data1$timegroup1==max(transformed_data1$timegroup1),]
  if(any(check1$event1!=0)) browser()
  
  ####2nd type
  data2 <- data[[2]]
  time_grid2 <- sort(unique(data2[data2[,"event2"]==1,][,"time2"]))
  #transformed_data <- as.data.frame(transformed_data)
  transformed_data2 <-
    survSplit(
      Surv(tstart2, time2, event2) ~ .,
      data2,
      cut = time_grid2,
      episode = "timegroup2"
    )
  ###tocheck the last time group
  check2 <- transformed_data2[transformed_data2$timegroup2==max(transformed_data2$timegroup2),]
  if(any(check2$event2!=0)) browser()
  
  #browser()
  return(list(as.data.frame(transformed_data1),
              as.data.frame(transformed_data2)))
}

#aa<-data.tran(data1,1)

panelEE <- function(alpha,data1,data2) {
  ###1st type
  zbar_mat1 <- data1 %>% filter(event1 == 1) %>% 
    group_by(timegroup1) %>%
    dplyr::select(-id,-count1,-tstart1,-time1,-event1
    ) %>%
    nest(-timegroup1, .key = Z1) %>%
    mutate(
      Z1 = map(Z1,  ~ as.matrix(as.data.frame(.))),
      expZ1 = map(Z1,  ~ exp(as.vector(. %*% alpha))),
      Zbar1 = map2(
        Z1,
        expZ1,
        .f = function(Zf, expZf)
          colSums(Zf * expZf) / sum(expZf)
      )
    ) %>%
    dplyr::select(timegroup1, Zbar1)
  
  eetemp1 <- data1 %>%
    filter(event1 == 1) %>%
    dplyr::select(-tstart1,-time1,-event1
    ) %>%
    nest(-id, -count1, -timegroup1, .key = Z1) %>%
    inner_join(zbar_mat1, by = "timegroup1") %>%
    mutate(
      Z1 = map(Z1,  ~ as.matrix(as.data.frame(.))),
      Zdiff1 = map2(Z1, Zbar1,  ~ as.vector(.x - .y)),
      temp1 = map2(count1, Zdiff1,  ~ as.vector(.x * .y))
    ) %>%
    pull(temp1)
  res1 <- Reduce("+", eetemp1)
  ###2nd type
  zbar_mat2 <- data2 %>% filter(event2 == 1) %>% 
    group_by(timegroup2) %>%
    dplyr::select(-id,-count2,-tstart2,-time2,-event2
    ) %>%
    nest(-timegroup2, .key = Z2) %>%
    mutate(
      Z2 = map(Z2,  ~ as.matrix(as.data.frame(.))),
      expZ2 = map(Z2,  ~ exp(as.vector(. %*% alpha))),
      Zbar2 = map2(
        Z2,
        expZ2,
        .f = function(Zf, expZf)
          colSums(Zf * expZf) / sum(expZf)
      )
    ) %>%
    dplyr::select(timegroup2, Zbar2)
  
  eetemp2 <- data2 %>%
    filter(event2 == 1) %>%
    dplyr::select(-tstart2,-time2,-event2
    ) %>%
    nest(-id, -count2, -timegroup2, .key = Z2) %>%
    inner_join(zbar_mat2, by = "timegroup2") %>%
    mutate(
      Z2 = map(Z2,  ~ as.matrix(as.data.frame(.))),
      Zdiff2 = map2(Z2, Zbar2,  ~ as.vector(.x - .y)),
      temp2 = map2(count2, Zdiff2,  ~ as.vector(.x * .y))
    ) %>%
    pull(temp2)
  res2 <- Reduce("+", eetemp2)
  res <- res1 + res2
}

#panelEE(c(0.5,0.5,0.5,0.5),data.tran(data.all.gen(nsub=100,obs.rate,tau=1),1))

############################sandwich######################################
get.sand<-function(alpha,data1,data2){
  ###1st type
  res1<-list()
  for(i in 1:(range(data1$timegroup1)[2]-1)){
    #print(i)
    #if(i==46) browser()
    countij1 <- data1[(data1$timegroup1==i & data1$event1 ==1),]$count1
    newcovarij1<-dplyr::select(data1[(data1$timegroup1==i & data1$event1 ==1),],-id,-count1,-tstart1,-time1,-event1,-timegroup1)
    res1.S0 <- exp(alpha%*%t(as.vector(newcovarij1)))
    res1.S0 <- sum(res1.S0)
    mu1hat <- sum(countij1)/res1.S0
    res1.S1<-colSums(exp(alpha%*%t(as.vector(newcovarij1)))*as.vector(newcovarij1))
    res1.S2<- lapply(1:dim(newcovarij1)[1], function(i) exp(alpha%*%t(as.vector(newcovarij1)))[i]*outer(unlist(newcovarij1[i,]),unlist(newcovarij1[i,])))
    res1.S2 <- Reduce('+',res1.S2)
    res1[[i]] <-list(mu1hat,res1.S0,res1.S1,res1.S2)
    if(is.na(mean(unlist(res1)))) browser()
  }
  ###2nd stype
  res2<-list()
  for(i in 1:(range(data2$timegroup2)[2]-1)){
    #print(i)
    #if(i==46) browser()
    countij2 <- data2[(data2$timegroup2==i & data2$event2 ==1),]$count2
    newcovarij2<-dplyr::select(data2[(data2$timegroup2==i & data2$event2 ==1),],-id,-count2,-tstart2,-time2,-event2,-timegroup2)
    res2.S0 <- exp(alpha%*%t(as.vector(newcovarij2)))
    res2.S0 <- sum(res2.S0)
    mu2hat <- sum(countij2)/res2.S0
    res2.S1<-colSums(exp(alpha%*%t(as.vector(newcovarij2)))*as.vector(newcovarij2))
    res2.S2<- lapply(1:dim(newcovarij2)[1], function(i) exp(alpha%*%t(as.vector(newcovarij2)))[i]*outer(unlist(newcovarij2[i,]),unlist(newcovarij2[i,])))
    res2.S2 <- Reduce('+',res2.S2)
    res2[[i]] <-list(mu2hat,res2.S0,res2.S1,res2.S2)
    if(is.na(mean(unlist(res1)))) browser()
  }
  
  return(list(res1,res2))
}

get.muhat <- function(alpha,data1,data2){
  ###1st type
  res1<-list()
  for(j in 1:(range(data1$timegroup1)[2]-1)){
    countij1<- data1[(data1$timegroup1==j & data1$event1 ==1),]$count1
    res1[[j]] <-sum(countij1)
  }
  time_grid1 <- sort(unique(data1[data1[,"event1"]==1,][,"time1"]))
  res1 <- cbind(time_grid1,res1,1:(range(data1$timegroup1)[2]-1))
  ###deleta the last one, where event==1, mu0t==0
  colnames(res1) <- c("time1","mu1hat","timegroup1")
  ###2nd type
  res2<-list()
  for(j in 1:(range(data2$timegroup2)[2]-1)){
    countij2<- data2[(data2$timegroup2==j & data2$event2 ==1),]$count2
    res2[[j]] <-sum(countij2)
  }
  time_grid2 <- sort(unique(data2[data2[,"event2"]==1,][,"time2"]))
  res2 <- cbind(time_grid2,res2,1:(range(data2$timegroup2)[2]-1))
  ###deleta the last one, where event==1, mu0t==0
  colnames(res2) <- c("time2","mu2hat","timegroup2")
  return(list(res1,res2))
}

######
get.sanda<-function(alpha,data1,data2,X.deriv1,X.deriv2){
  #1st type
  res.temp1<-list()
  for(i in 1:(range(data1$timegroup1)[2]-1)){
    countij1<- data1[(data1$timegroup1==i & data1$event1 ==1),]$count1
    res1<-sum(countij1)*X.deriv1[[i]] # MODIFIED add sum to countij
    res.temp1[[i]]<-res1
  }
  #res<-Reduce('+',res)
  res1 <- (-1)*reduce(res.temp1[-length(res.temp1)],`+`)
  ##2nd type
  res.temp2<-list()
  for(i in 1:(range(data2$timegroup2)[2]-1)){
    countij2<- data2[(data2$timegroup2==i & data2$event2 ==1),]$count2
    res2<-sum(countij2)*X.deriv2[[i]] # MODIFIED add sum to countij
    res.temp2[[i]]<-res2
  }
  #res<-Reduce('+',res)
  res2 <- (-1)*reduce(res.temp2[-length(res.temp2)],`+`)
  return(list(res1,res2))
}

#
get.sandb <- function(alpha, data1,data2, X.bar1,X.bar2, mu1hat,mu2hat) {
  #1st type
  #browser()
  res1 <- list()
  #data1 <- data1[data1$timegroup1<=length(X.bar1),]
  ###2nd type
  res2 <- list()
  #data2 <- data2[data2$timegroup2<=length(X.bar2),]
  res <- list()
  for (i in 1:nsub) {
    #print(i)
    #if(i==60) browser()
    data.temp1 <- data1[(data1$id == i & data1$event1 == 1), ]
    countij1 <- data.temp1$count1
    newcovarij1 <-
      dplyr::select(data.temp1, -id, -count1, -tstart1, -time1, -event1, -timegroup1)
    timegroup_ind1 <- data.temp1$timegroup1
    comp1 <- unlist(mu1hat[timegroup_ind1]) * exp(alpha %*% t(newcovarij1))
    if(max(data.temp1$timegroup1)>length(X.bar1)) browser()
    res1.1 <-
      lapply(1:dim(data.temp1)[1], function(i)
        (countij1[i] - comp1[i]) * (newcovarij1[i, ] - X.bar1[[timegroup_ind1[i]]]))
    
    res1.2 <-
      lapply(1:dim(data.temp1)[1], function(i)
        (countij1[i] ) * (newcovarij1[i, ] - X.bar1[[timegroup_ind1[i]]]))
    
    res1temp <-
      Reduce('+', res1.1) ###for every subject, sum for all time group
    
    #print(i)
    #if(i==nsub) browser()
    data.temp2 <- data2[(data2$id == i & data2$event2 == 1), ]
    countij2 <- data.temp2$count2
    newcovarij2 <-
      dplyr::select(data.temp2, -id, -count2, -tstart2, -time2, -event2, -timegroup2)
    timegroup_ind2 <- data.temp2$timegroup2
    comp2 <- unlist(mu2hat[timegroup_ind2]) * exp(alpha %*% t(newcovarij2))
    if(max(data.temp2$timegroup2)>length(X.bar2)) browser()
    
    res2.1 <-
      lapply(1:dim(data.temp2)[1], function(i)
        (countij2[i] - comp2[i]) * (newcovarij2[i, ] - X.bar2[[timegroup_ind2[i]]]))
    
    res2.2 <-
      lapply(1:dim(data.temp2)[1], function(i)
        (countij2[i] ) * (newcovarij2[i, ] - X.bar2[[timegroup_ind2[i]]]))
    
    res2temp <-
      Reduce('+', res2.1) ###for every subject, sum for all time group
    
    ####add together
    restemp <- res1temp + res2temp
    
    res[[i]] <- outer(unlist(restemp), unlist(restemp)) #for every subject
  }
  #browser()
  res <- Reduce('+', res) #sum for all subjects
  return(res)
}
