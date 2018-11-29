library(MASS) #mvrnorm
library(nnet) #which.is.max
#library(doMC) #parallel
#library(foreach) #parallel
#set.seed(1000)
#############################################
#    Initialize the space for chromosomes   # 
#############################################

Init <- function(n.chor,n.Island){
  MDL <<- NULL # discription length
  Tau <<- NULL
  temp <- NULL
  for (j in 1:n.chor) {temp[[j]] <- rep(NA, K)}
  for (i in 1:n.Island) {
    Tau[[i]] <<- temp # [[n.Island]] [[n.chor]] [k] (k = 1, 2, ..., K)
    MDL[[i]] <<- rep(NA, n.chor) # [[n.Island]] [k] (k = 1, 2, ..., n.chor)
  }
  Sorted.MDL <<- MDL # To store sorted discription length
  Elite.Tau <<- NULL 
  Elite.MDL <<- NULL
  for (i in 1:n.Island) {
    Elite.Tau[[i]] <<- rep(NA, K)
    Elite.MDL[i] <<- NA
  }
}

#################################
#     Sorting MDL in Island i   # 
#################################

SIG<-function(i) {
  Sorted.MDL[[i]] <<- sort(MDL[[i]])
}

##### Calculate MDL for single series#####

MDL_Single <- function(tau,y) { # tau: real number, y: univariate series of length n
  if (tau==0) {
    return(- log2(n) + 0.5 * log2(n) + 0.5 * (n * log(2*pi) + sum ((y - mean(y))^2)) * log2(exp(1))) # To be discussed: should we minus log2(n) for tau=0?
  } else {
    ybar_1 <- mean(y[1:tau]) 
    ybar_2 <- mean(y[(tau+1):n])
    ybar <- c(rep(ybar_1,tau),rep(ybar_2,n-tau))
    return(0.5 * (log2(tau) + log2(n-tau)) + 0.5 * (n * log(2*pi) + sum ((y - ybar)^2)) * log2(exp(1)))
  }
}

##### Choose the tau that minimize MDL for each single sequence ####

Choose_Tau <- function(Tau_Set, y){ # Tau_Set: all possible tau's for choice, y: univariate series
  mdl_temp <- sapply(Tau_Set,MDL_Single, y=y)
  return(list(tau = which.is.max(-mdl_temp), partial_mdl = min(mdl_temp)))
  #return(which.is.max(-mdl_temp))
}

##### Calculate MDL for given sample and change points estimators #####
MDL_Total <- function(Tau_Set){ # Tau_Set: all possible tau's for choice Y: K sequences with length n
  if ((length(Tau_Set) == 1) & (Tau_Set[1] == 0)) {CL_Part1 <- (log2(K) + K)} else {
    m <- length(Tau_Set) # number of different change points
    CL_Part1 <- (log2(K) + K * log2(m) + m * log2(n)) # Codelength for change point parameters
  }
    CL_Part2 <- sum(sapply(1:K,function(i){
      Choose_Tau(Tau_Set,Y[,i])$partial_mdl
    }))#Codelength for mu's and residuals
  return(CL_Part1 + CL_Part2)
}

##### Calculate MDL for l-th chromosome #####
Cal.MDL <- function(l) {
  i<-l%%n.Island+1
  j<-l%/%n.Island+1
  MDL[[i]][j] <<- MDL_Total(Tau[[i]][[j]])
} 
##### Check if Tau_Set has minimum mesh greater than delta #####
#Valid_Tau <- function(Tau_Set,delta){
#  if (length(Tau_Set)<=1) {return(TRUE)} else {
#  s <- sort(Tau_Set)
#  return(min(diff(s)) >= delta)
#  }
#}

##### Generate a new random chromosome #####
CP.chromosome <- function(delta){
  m <- min(c(rpois(n=1,lambda=2)+1,n,K)) # No. of different change points
#  p <- FALSE
#  while (!p) {
    Tau_Set <- sample(0:n,m)
    Tau_Set <- Tau_Set[unique(sapply(1:K, function(i) {Choose_Tau(Tau_Set,Y[,i])$tau}))] # Remove Redundant tau's in Tau_Set
    Tau_Set[Tau_Set>=n] <- 0
#    p <- Valid_Tau(Tau_Set,delta)
#  }
#  Tau_Set[(1:length(Tau_Set)) * ((Tau_Set >= n - delta) + (Tau_Set <= delta))] <- 0 # to be discussed
  return(Tau_Set)
}

##### Store Chromosomes #####
chromosome <- function(l,delta){
  i<-l%%n.Island+1
  j<-l%/%n.Island+1
  Tau[[i]][[j]] <<- CP.chromosome(delta)
}

###########################################
#     CM Parents Selection/Generation     #
###########################################
CM <- function(l,pc,delta){
  i<-l %% n.Island + 1
  j<-l %/% n.Island + 1  
  
  rnumber <- 1/(1:n.chor) #prob for parent choosing
  pnumber <- sample(n.chor,1,prob=rnumber)
  index <- (MDL[[i]]==Sorted.MDL[[i]][pnumber])*(1:n.chor)
  C1 <- NULL # The first parent
  C1$Tau <-Tau[[i]][index][[1]]
 
    check1 <- runif(1, 0, 1)
    if (check1 < pc) { 
      pnumber <- sample(n.chor, 1, prob = rnumber)			# rnumber is 1/n.chor for decreasing prob
      index <- (MDL[[i]] == Sorted.MDL[[i]][pnumber]) * (1:n.chor)
      C2 <- NULL # The second parent
      C2$Tau <- Tau[[i]][index][[1]]
    } else {
      C2 <- NULL # The second parent
      C2$Tau <-CP.chromosome(delta)
    }
  temp.TAR <- CO.TAR(C1,C2,delta)
  Tau[[i]][[j]] <- temp.TAR
}

#####################
#     Crossover     #
#####################

CO.TAR <- function(C1,C2,delta){ # AR order part can be finetuned
  
  Tau1 <- C1$Tau
  Tau2 <- C2$Tau
  
  Tau_combine <- c(Tau1, Tau2)
  p <- FALSE
#  while (!p) {
    index <- (1:length(Tau_combine)) * sample(c(0,1),length(Tau_combine),replace=T)
    Tau_Set <- unique(Tau_combine[index])
    Tau_Set[Tau_Set>=n] <- 0
#    p <- Valid_Tau(Tau_Set,delta)
#  }
  return(Tau_Set)
}

###############################################
#    Save elite chromosome for all islands    #
###############################################

elitesaving <- function(){
  sapply(1:n.Island, function(i) {
    Elite.MDL[i] <<- Sorted.MDL[[i]][1]
    Elite.Tau[[i]] <<- Tau[[i]][MDL[[i]] == Elite.MDL[i]][[1]]
  })
}

#########################
#     Elitest Step      #
#########################

ELI<-function(l) {
  Tau[[l]][MDL[[l]]==Sorted.MDL[[l]][n.chor]][[1]] <<- Elite.Tau[[l]]
  MDL[[l]][MDL[[l]]==Sorted.MDL[[l]][n.chor]][1] <<- Elite.MDL[[l]]
  for (j in 1:round(n.chor/3)){
    # find the index of j-th worst chromosome, if there are more than one, select one randomly
    temp=(MDL[[l]]==Sorted.MDL[[l]][n.chor-j])*(1:length(MDL[[l]]))
    temp0=temp[temp>0]
    temp=temp0[sample(length(temp0),1)]
    
    tau <- Elite.Tau[[l]] + sample(-5:5,length(Elite.Tau[[l]]),replace=T)
    tau[tau<0] <- 0
    tau[tau>=n] <- 0
    Tau[[l]][[temp]] <<- tau
  }
}

#############################
#     Migration Function    #
#############################

#1st island good-> 2nd island bad ,..., n.Islandth island good -> 1st island bad

MIG<-function(i) {
  if (i>1) {
    Tau[[i]][MDL[[i]]==Sorted.MDL[[i]][n.chor]][[1]]<<-Tau[[i-1]][MDL[[i-1]]==Sorted.MDL[[i-1]][1]][[1]]
    Tau[[i]][MDL[[i]]==Sorted.MDL[[i]][n.chor-1]][[1]]<<-Tau[[i-1]][MDL[[i-1]]==Sorted.MDL[[i-1]][2]][[1]]
  } else {
    Tau[[i]][MDL[[i]]==Sorted.MDL[[i]][n.chor]][[1]]<<-Tau[[n.Island]][MDL[[n.Island]]==Sorted.MDL[[n.Island]][1]][[1]]
    Tau[[i]][MDL[[i]]==Sorted.MDL[[i]][n.chor-1]][[1]]<<-Tau[[n.Island]][MDL[[n.Island]]==Sorted.MDL[[n.Island]][2]][[1]]
  }
}

####################
#    Main Flow     #
####################

GA <- function(n.Island,n.chor,pc,delta) {
  
  # Initialzie storage space
  Init(n.chor,n.Island)
  
  # Generate the population
#  print("Generate Chromosome")
  sapply(0:(n.Island*n.chor-1),chromosome,delta)
#  print(c("Find MDL",date()))
  sapply(0:(n.Island*n.chor-1),Cal.MDL)
#  print(c("Finish MDL",date()))
  sapply(1:n.Island,SIG)
  elitesaving() 
  count <- 0
  last.MDL <- 0
  # Main process
  for (b in 1:300) {
#    print(c("Generation",b,"start",date()))
    sapply(0:(n.Island*n.chor-1),Cal.MDL)
#     print(c("Generation",b,"MDL Calculated",date()))
    sapply(1:n.Island,SIG)
    # Doing Elite Step for Monotonic increasing of Chromosome
    sapply(1:n.Island,ELI)
    # Do SIG function again as the sigma value changes due to the elite step
    sapply(0:(n.Island*n.chor-1),Cal.MDL)
    sapply(1:n.Island,SIG)
    # Saving the Elite Chromosome
    elitesaving()
#    print(c("Best.Tau",Elite.Tau[Elite.MDL==min(Elite.MDL)][[1]]))
#    print(c("Best.MDL",min(Elite.MDL)[1]))
    if (last.MDL == min(Elite.MDL)) { 
      if (count==10) {break} else {
        count <- count+1
      }
    } else {
      count <- 0
      last.MDL <- min(Elite.MDL)
    }
    sapply(0:(n.Island*n.chor-1),CM, pc=pc, delta=delta)
    if (b%%4==0){
#      print("Do migration every four generations")
      sapply(0:(n.Island*n.chor-1),Cal.MDL)
      sapply(1:n.Island,SIG)
      sapply(1:n.Island,MIG)
    }
    # Close the repeating b loops
  }
  
  # Do one more calculation for finding threshold
  sapply(0:(n.Island*n.chor-1),Cal.MDL)
  sapply(1:n.Island,SIG)
  elitesaving()
  MDL.est<-min(Elite.MDL)
  Tau.est<-Elite.Tau[Elite.MDL==min(Elite.MDL)][[1]]
  Tau.est <- Tau.est[unique(sapply(1:K, function(i) {Choose_Tau(Tau.est,Y[,i])$tau}))]
	# One more fine-tuning step to see whether 0 should be included
  if (0 %in% Tau.est) {
    Tau.temp <- Tau.est[Tau.est!=0] 
  } else {Tau.temp <- c(0,Tau.est)}
  MDL.temp <- MDL_Total(Tau.temp)
  if (MDL.temp < MDL.est) {Tau.est <- Tau.temp; MDL.est <- MDL.temp}
  Tau.est <- sort(Tau.est)
  print(c("result, Tau=",Tau.est))
  return(list(Tau.est = Tau.est, MDL.est = MDL.est))
}

##### Generate the independently distributed sequence of length n with exactly one change point #####
DataGen <- function(Tau,Mu,n){ 
  for (ii in 1:K) {
        tau <- Tau[ii]
        mu1 <- Mu[1,ii]
	mu2 <- Mu[2,ii]
  	if ((tau > n-1) || (tau < 0)) {print("Error, change point out of range!")}
  	y1 <- rnorm(n = tau, mean = mu1, sd = sigma)
 	y2 <- rnorm(n = n-tau, mean = mu2,sd = sigma)
  	Y[,ii] <<- c(y1,y2)
  }
}
 
#set model parameters

n <- 3000 #length of sequences
temp <- read.csv("inputB2.csv",header = T)
Tau0 <- n * temp[,1] #change point of each sequence
Mu0 <- t(temp[,2:3]) #population mean for each subsequence
sigma <- 1 #standard deviation of the sequences

delta <- n * 0.02
K <- length(Tau0)

#set parameters for genetic algorithm 
n.Island <- 20
n.chor <- 40
pc <- 0.8
Y <- matrix(data = 0, nrow = n, ncol = K) # Y store all the sequence

# plot all sequences
#for (i in 1:K){
#  plot(Y[,i],type="l",ylab = paste("y_", i))
#}
#registerDoMC(100)
#rept <- 200

#foreach(i=1:rept) %dopar% {
  
  DataGen(Tau0,Mu0,n)
  t1 <- Sys.time()
  Result <- GA(n.Island,n.chor,pc,delta)
  l <- sapply(1:K,function(i){Choose_Tau(Tau_Set=Result$Tau.est,Y[,i])$tau})
  Tau_Hat <- Result$Tau.est[l]  
  t2 <- Sys.time()
  print(Tau_Hat/n)
  print(t1)
  print(t2)
  print(t2-t1)
#  write.table(t(c(length(Result$Tau.est),(Result$Tau.est),l)),file="modelA2.csv",append=T,quote=F,row.names=F, col.names=F, sep = ",")
 # write.table(t(Tau_Hat/n),file="model1.csv",append=T,quote=F,row.names=F, col.names=F, sep = ",")
#}
