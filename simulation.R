### simulating an SIR model

# function to run the simulation
SIR_simulate <- function(N,beta,gamma,seed){
  
  # setting the seed
  set.seed(seed)
  
  # vectors for infection and recovery times
  i <- rep(Inf, times=N)
  i[1] <- 0
  i_count <- 2
  r <- rep(Inf, times=N)
  r_count <- 1
  
  # starting values and vectors for the three categories
  S <- N-1
  S_vec <- c()
  S_vec[1] <- S
  I <- 1
  I_vec <- c()
  I_vec[1] <- I
  R <- 0
  R_vec <- c()
  R_vec[1] <- R
  
  # initial time and time vector
  t <- 0
  event_times <- c()
  event_times[1] <- 0
  event_no <- 2
  
  # running the model until I=0
  while(I != 0){
    
    # generate new time
    t_star <- rexp(1, beta*S*I/N + gamma*I)
    t <- t+t_star
    
    # decide which event happens
    u <- runif(1)
    prob <- (beta*S*I/N)/(beta*S*I/N + gamma*I)
    if(u<prob){
      event<-"SI"
    }
    else{
      event<-"IR"
    }
    
    # the event happens
    if(event=="SI"){
      S <- S-1
      I <- I+1
      i[i_count] <- t
      i_count <- i_count+1
    }
    else{
      I <- I-1
      R <- R+1
      r[r_count] <- t
      r_count <- r_count+1
    }
    
    # updating event vectors
    S_vec[event_no] <- S
    I_vec[event_no] <- I
    R_vec[event_no] <- R
    event_times[event_no] <- t
    event_no <- event_no+1
    
  }
  
  # returning values
  return(list(i=i,r=r,S=S_vec,I=I_vec,R=R_vec,t=event_times))
  
}

# running the simulation
sim <- SIR_simulate(100,2,1,11)

# resulting infection and recovery times
sim$i
sim$r

# plotting the results with base plot
plot(sim$t,sim$S,type="l",xlim=c(0,sim$t[length(sim$t)]),ylim=c(0,100))
lines(sim$t,sim$I,col="red")
lines(sim$t,sim$R,col="blue")

# plotting the results with ggplot
library(ggplot2)
df_S <- as.data.frame(cbind(sim$t,sim$S))
df_I <- as.data.frame(cbind(sim$t,sim$I))
df_R <- as.data.frame(cbind(sim$t,sim$R))
sim_plot <- ggplot() + geom_line(aes(x=V1,y=V2),df_S) + geom_line(aes(x=V1,y=V2),df_I,col="red") + geom_line(aes(x=V1,y=V2),df_R,col="blue")
sim_plot <- sim_plot + xlab("Time") + ylab("Individuals") + ggtitle(expression(paste("SIR Model (", beta, "=2, ", gamma, "=1)")))
print(sim_plot)