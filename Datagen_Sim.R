setwd("C:/Users/sjoo001/Dropbox/Research/RMSD Item Fit/Item Information/Dichotomous/Sim Results")
library(gdata)

# Number of Replications
nsim <- 200

# DIF proportion (5%, 15%, 30%)
cond1 <- "5%"

# DIF size (small, medium, large)
cond2 <- "small"

# DIF item-level (easy, mod, hard)
cond3 <- "mod"

# Number of Items
J <- 40

# Sample Size
N <- 500

# Number of Groups
G <- 80

# Number of DIF Groups
gdif <- 1

# Number of DIF Itmes 
if(cond1=="5%") ndif <- 2
if(cond1=="15%") ndif <- 6
if(cond1=="30%") ndif <- 12

# Size of DIF
a_dif <- 0
if(cond2=="small") b_dif <- .48
if(cond2=="medium") b_dif <- .77
if(cond2=="large") b_dif <- 1.28

# DIF Type: Uniform (b_dif), Nonuniform (a_dif)
# 1 = Nonuniform DIF (a shift), 2 = Uniform DIF(b shift), 3 = Both
diftype <- 2

# Response Type: Dichotomous, 3-option Polytomous
# 1 = Dichotomous, 2 = 3-option Polytomous
dattype <- 1

# 2PL function
twopl <- function(a,b,theta){
  prob=exp(1.7*(a*(theta-b)))/(1+exp(1.7*(a*(theta-b))))
  return(prob)
}

# GPCM function 
gpcm <- function(a,b,d,theta,cat){
  num <- rep(0,cat)
  for(i in 1:cat){
    if(i==1){
      z <- 1.7*a*(theta-b)
      num[i] <- exp(z)
    }
    if(i>1){ 
      z <- z + 1.7*a*(theta-b+d[i-1])
      num[i] <- exp(z)
    }
  }
  deno <- sum(num)
  p <- num/deno
  return(p)
}

########################
### Start Simulation ###
########################
for(sim in 1:nsim){
  
  # generate item parameters
  if(dattype==1){
     a <- runif(J,min=0.80,max=1.82)
     b.easy <- runif(13,min=-1.91,max=-0.42)
     b.mode <- runif(14,min=-0.42,max=0.54)
     b.hard <- runif(13,min=0.54,max=2.65)
     if(cond3=="easy") b <- c(b.easy,b.mode,b.hard)
     if(cond3=="mod")  b <- c(b.mode,b.easy,b.hard) 
     if(cond3=="hard") b <- c(b.hard,b.easy,b.mode)
     
    # Item parameters are different (some variance across groups) -> more realistic situation
    # This can hold Type I error at .05 with RMSD = .10
     a <- a%*%matrix(1,nc=G)+matrix(1,nr=J)%*%rnorm(G,mean=0,sd=.04)
     b <- b%*%matrix(1,nc=G)+matrix(1,nr=J)%*%rnorm(G,mean=0,sd=.02)
  }
  
  if(dattype==2){
    a <- runif(J,min=0.75,max=2.25)
    b <- runif(J,min=-1.00,max=1.00)
    d <- runif(J,min=-1.00,max=1.00)
    
    # Item parameters are different (some variance across groups) -> more realistic situation
    # This can hold Type I error at .05 with RMSD = .10
    a <- a%*%matrix(1,nc=G)+matrix(1,nr=J)%*%rnorm(G,mean=0,sd=.000004)
    b <- b%*%matrix(1,nc=G)+matrix(1,nr=J)%*%rnorm(G,mean=0,sd=.000002)
    d <- d%*%matrix(1,nc=G)+matrix(1,nr=J)%*%rnorm(G,mean=0,sd=.000002)
  }
  
  # generate group means (randomly generate group means)
  theta <- matrix(0,nr=N,nc=G)
  genmu <- rep(0,G)
  for(g in 2:G){
    if(g==2) genmu[g] <- -1.24 
    else genmu[g] <- runif(1,min=-0.75,max=1)
  }
  for(g in 1:G){
    theta[,g] <- rnorm(N,genmu[g],1)   
  }
  
  # Create DIF
  # Half of groups affected by DIF parameter increased
  # The other half of groups affected by DIF parameter decreased
  if(diftype==1){
    for(j in 1:ndif){
      for(g in 1:gdif){
        a[j,(g+1)] <- a[j,(g+1)]+a_dif 
      }
    }
  }
  if(diftype==2){
    for(j in 1:ndif){
      for(g in 1:gdif){
        b[j,(g+1)] <- b[j,(g+1)]+b_dif 
      }
    }
  }
  if(diftype==3){
    for(j in 1:ndif){
      for(g in 1:gdif){
        a[j,(g+1)] <- a[j,(g+1)]+a_dif 
        b[j,(g+1)] <- b[j,(g+1)]+b_dif 
      }
    }
  }
  
  gcount <- rep(101,N)
  for(g in 1:G){
    if(dattype==1){
      tmpX <- matrix(NA,nr=N,nc=J)
      prob <- matrix(NA,nr=N,nc=J)
      for(i in 1:N){
        for(j in 1:J){
          prob[i,j] <- twopl(a[j,g],b[j,g],theta[i,g])
        }
      }
      
      rand <- matrix(runif(N*J,0,1),nr=N,nc=J)
      for(i in 1:N){
        for(j in 1:J){
          if(prob[i,j]>rand[i,j]){ tmpX[i,j] <- 1 }
          else{tmpX[i,j] <- 0}
        }
      }
      
      tmpX <- cbind(gcount,tmpX)
      gcount <- gcount+1
      
      if(g==1){ X <- tmpX }
      if(g>1){ X <- rbind(X,tmpX)}
    }
    if(dattype==2){
      tmpX <- matrix(NA,nr=N,nc=J)
      prob0 <- matrix(NA,nr=N,nc=J)
      prob1 <- matrix(NA,nr=N,nc=J)
      prob2 <- matrix(NA,nr=N,nc=J)
      for(i in 1:N){
        for(j in 1:J){
          prob0[i,j] <- gpcm(a[j,g],b[j,g],c(d[j,g],-d[j,g]),theta[i,g],3)[1]
          prob1[i,j] <- gpcm(a[j,g],b[j,g],c(d[j,g],-d[j,g]),theta[i,g],3)[2]
          prob2[i,j] <- gpcm(a[j,g],b[j,g],c(d[j,g],-d[j,g]),theta[i,g],3)[3]
        }
      }
      rand <- matrix(runif(N*J,0,1),nr=N,nc=J)
      for(i in 1:N){
        for(j in 1:J){
          if( (prob0[i,j]+prob1[i,j]) < rand[i,j] ) tmpX[i,j] <- 2 
          else if(prob0[i,j]<rand[i,j]) tmpX[i,j] <- 1 
          else tmpX[i,j] <- 0
        }
      }
      
      tmpX <- cbind(gcount,tmpX)
      gcount <- gcount+1
      
      if(g==1) X <- tmpX 
      if(g>1) X <- rbind(X,tmpX)
    }
  }
  write.fwf(X,paste(cond1,"_",cond2,"_",cond3,"/","SimData",sim,".dat",sep=""),colnames=F,sep="")
  
} # End of Replications
