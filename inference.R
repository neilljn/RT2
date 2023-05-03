### bayesian inference from recovery times

# function for the double sum inside other stuff
sigma <- function(i,r,n,N){
  double_sum <- 0
  for(a in 1:n){
    for(b in 1:n){
      double_sum <- double_sum + min(r[a],i[b]) - min(i[a],i[b])
    }
    for(b in (n+1):N){
      double_sum <- double_sum + r[a] - i[a]
    }
  }
  return(double_sum)
}

# function for the number of infectives just before infection j
infectives <- function(j,i,r,n,N){
  current <- i[j]
  num <- length(i[i<current])-length(r[r<current])
  return(num)
}

# function for the log of the conditional pdf of i
log_i <- function(i,r,beta,gamma,n,N){
  sum_log_infectives <- 0
  kappa <- which.min(i)
  for(j in 1:n){
    if(j!=kappa){
      sum_log_infectives <- sum_log_infectives + log(infectives(j,i,r,n,N))
    }
  }
  log_pdf <- sum_log_infectives - (beta*sigma(i,r,n,N)/N) - (gamma*(sum(r-i)))
  return(log_pdf)
}

# function to update the infection j (with log base)
MH_update <- function(j,i,r,beta,gamma,n,N){
  
  # generate new infection period length
  x <- rexp(1,gamma)
  i_new <- r[j] - x
  
  # create i_star
  i_star <- i
  i_star[j] <- i_new
  
  # calculate log alpha
  log_alpha <- log_i(i_star,r,beta,gamma,n,N) + dexp(r[j]-i[j],gamma,log=TRUE) - log_i(i,r,beta,gamma,n,N) - dexp(r[j]-i_new,gamma,log=TRUE)
  log_alpha <- min(0,log_alpha)
  if(is.nan(log_alpha)){
    return(i)
  }
  
  # accept or reject
  u <- runif(1)
  if(log(u)<log_alpha){
    return(i_star)
  }
  else{
    return(i)
  }
  
}

# overall MCMC algorithm
SIR_inference <- function(i_0,r,beta_0,gamma_0,a_beta,b_beta,a_gamma,b_gamma,M,seed){
  
  # setting the seed
  set.seed(seed)
  
  # length of data, with and with infinities
  N <- length(r)
  n <- length(r[r!=Inf])
  
  # removing infinities from i_0 and r
  r <- r[r!=Inf]
  i_0 <- i_0[i_0!=Inf]
  
  # making the current values equal to the initial values
  i <- i_0
  beta <- beta_0
  gamma <- gamma_0
  
  # creating vectors for stored values
  beta_vec <- c()
  gamma_vec <- c()
  sumi_vec <- c()
  i51_vec <- c()
  ll_vec <- c()
  beta_vec[1] <- beta
  gamma_vec[1] <- gamma
  sumi_vec[1] <- sum(i)
  i51_vec[1] <- i[51]
  
  # each MCMC iteration (M times)
  for(k in 1:M){
    
    # update beta
    beta <- rgamma(1,n+a_beta-1,(sigma(i,r,n,N)/N)+b_beta)
    beta_vec[k+1] <- beta
    
    # update gamma
    gamma <- rgamma(1,n+a_gamma,sum(r-i)+b_gamma)
    gamma_vec[k+1] <- gamma
    
    # update i
    for(j in 1:n){
      i <- MH_update(j,i,r,beta,gamma,n,N)
    }
    sumi_vec[k+1] <- sum(i)
    i51_vec[k+1] <- i[51]
    
    # keeping track of the log likelihood
    ll_vec[k] <- log_i(i,r,beta,gamma,n,N)
    
  }
  
  # returning values
  return(list(beta=beta_vec,gamma=gamma_vec,sumi=sumi_vec,i51=i51_vec,loglikelihood=ll_vec))
  
}

# running the inference
inf <- SIR_inference((sim$r)/2,sim$r,0,0,1,0.0001,1,0.0001,1000,1)

# plotting the trace plot with base plot
plot(inf$beta,type="l")
lines(1:1000,rep(2,times=1000),col="red")
plot(inf$gamma,type="l")
lines(1:1000,rep(1,times=1000),col="red")
plot(inf$sumi,type="l")
lines(1:1000,rep(sum(sim$i[1:76]),times=1000),col="red")
plot(inf$i51,type="l")
lines(1:1000,rep(sim$i[51],times=1000),col="red")

# removing burn-in
betas <- inf$beta[-(1:20)]
gammas <- inf$gamma[-(1:20)]
sumis <- inf$sumi[-(1:20)]
i51s <- inf$i51[-(1:20)]

# plotting the histograms with burn-in removed with base plot
hist(betas)
hist(gammas)
hist(sumis)
hist(i51s)

# average parameter estimates with burn-in removed (and comparisons to real values)
mean(betas)
print(2)
mean(gammas)
print(1)
mean(sumis)
print(sum(sim$i[1:76]))
mean(i51s)
print(sim$i[51])
mean(betas/gammas)
print(2)

# plotting the trace plots with ggplot
library(ggplot2)
df_beta <- as.data.frame(cbind(1:1001,inf$beta))
trace_beta <- ggplot() + geom_line(aes(x=V1,y=V2),df_beta) + geom_hline(yintercept=2, color = "red", linewidth=1.2)
trace_beta <- trace_beta + xlab("Iteration") + ylab(expression(beta)) + ggtitle(expression(paste("Trace Plot of ", beta)))
print(trace_beta)
df_gamma <- as.data.frame(cbind(1:1001,inf$gamma))
trace_gamma <- ggplot() + geom_line(aes(x=V1,y=V2),df_gamma) + geom_hline(yintercept=1, color = "red", linewidth=1.2)
trace_gamma <- trace_gamma + xlab("Iteration") + ylab(expression(gamma)) + ggtitle(expression(paste("Trace Plot of ", gamma)))
print(trace_gamma)
df_i51 <- as.data.frame(cbind(1:1001,inf$i51))
trace_i51 <- ggplot() + geom_line(aes(x=V1,y=V2),df_i51) + geom_hline(yintercept=sim$i[51], color = "red", linewidth=1.2)
trace_i51 <- trace_i51 + xlab("Iteration") + ylab("Infection Time of Individual 51") + ggtitle("Trace Plot of Infection Time of Individual 51")
print(trace_i51)
df_sumi <- as.data.frame(cbind(1:1001,inf$sumi))
trace_sumi <- ggplot() + geom_line(aes(x=V1,y=V2),df_sumi) + geom_hline(yintercept=sum(sim$i[1:76]), color = "red", linewidth=1.2)
trace_sumi <- trace_sumi + xlab("Iteration") + ylab("Sum of Infection Times") + ggtitle("Trace Plot of Sum of Infection Times")
print(trace_sumi)

# plotting the histograms with burn-in removed with ggplot
df_beta_hist <- as.data.frame(betas)
hist_beta <- ggplot() + geom_histogram(aes(x=betas),df_beta_hist,bins=10) + geom_vline(xintercept=2, color = "red", linewidth=1.2) + geom_vline(xintercept=mean(betas), color = "blue", linewidth=1.2)
hist_beta <- hist_beta + xlab(expression(beta)) + ylab("Quantity") + ggtitle(expression(paste("Histogram of ", beta)))
print(hist_beta)
df_gamma_hist <- as.data.frame(gammas)
hist_gamma <- ggplot() + geom_histogram(aes(x=gammas),df_gamma_hist,bins=10) + geom_vline(xintercept=1, color = "red", linewidth=1.2) + geom_vline(xintercept=mean(gammas), color = "blue", linewidth=1.2)
hist_gamma <- hist_gamma + xlab(expression(gamma)) + ylab("Quantity") + ggtitle(expression(paste("Histogram of ", gamma)))
print(hist_gamma)
df_i51_hist <- as.data.frame(i51s)
hist_i51 <- ggplot() + geom_histogram(aes(x=i51s),df_i51_hist,bins=10) + geom_vline(xintercept=sim$i[51], color = "red", linewidth=1.2) + geom_vline(xintercept=mean(i51s), color = "blue", linewidth=1.2)
hist_i51 <- hist_i51 + xlab("Infection Time of Individual 51") + ylab("Quantity") + ggtitle("Histogram of Infection Time of Individual 51")
print(hist_i51)
df_sumi_hist <- as.data.frame(sumis)
hist_sumi <- ggplot() + geom_histogram(aes(x=sumis),df_sumi_hist,bins=10) + geom_vline(xintercept=sum(sim$i[1:76]), color = "red", linewidth=1.2) + geom_vline(xintercept=mean(sumis), color = "blue", linewidth=1.2)
hist_sumi <- hist_sumi + xlab("Sum of Infection Times") + ylab("Quantity") + ggtitle("Histogram of Sum of Infection Times")
print(hist_sumi)
df_r0_hist <- as.data.frame(betas/gammas)
hist_r0 <- ggplot() + geom_histogram(aes(x=betas/gammas),df_r0_hist,bins=10) + geom_vline(xintercept=2, color = "red", linewidth=1.2) + geom_vline(xintercept=mean(betas/gammas), color = "blue", linewidth=1.2)
hist_r0 <- hist_r0 + xlab("Basic Reproduction Number") + ylab("Quantity") + ggtitle("Histogram of Basic Reproduction Number")
print(hist_r0)