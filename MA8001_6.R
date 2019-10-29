library(ggplot2)
# 1a)
marg_prob = rep(NA, 50)
marg_prob[1] = 0.01
for(i in 2:50){
  marg_prob[i] = marg_prob[i-1] + 0.05*(1-marg_prob[i-1])
}

# 1c)
set.seed(2606)
tau = 0.3
x = rep(NA,50)
x[1] = rbinom(1, 1, marg_prob[1])
y = rep(NA, 50)
y[1] = rnorm(1, x[1], sd = tau)
for (i in 2:50){
  if(x[i-1] == 1){x[i] =1}
  else{
    if(runif(1) < 0.05){
      x[i] = 1
    }
    else{
      x[i] = 0
    }
  }
  if(i == 20 | i==30){
    y[i] = rnorm(1, x[i], sd = tau)
  }
  else{
    y[i] = rnorm(1, x[i], sd = 100*tau)
  }
}

y = rnorm(50, 0, 100*tau)
y[20] = 0.2
y[30] = 0.7
plot(x)
plot(y)

fb_algo <- function(y_local, tau_local, k){
  y = y_local
  
  # First initiate normality constants and forward probabilities
  norm_const = c(1/(dnorm(y[1],0,100*tau_local)*0.99+dnorm(y[1],1,100*tau_local)*0.01), rep(NA,49))
  forward_prob = matrix(0,nrow=50,ncol = 2)
  if (1%in%k){
    forward_prob[1,1] = norm_const[1]*0.99*dnorm(y[1],0,tau_local)
    forward_prob[1,2] = norm_const[1]*0.01*dnorm(y[1],1,tau_local)
  }
  else{
  forward_prob[1,1] = norm_const[1]*0.99*dnorm(y[1],0,100*tau_local)
  forward_prob[1,2] = norm_const[1]*0.01*dnorm(y[1],1,100*tau_local)
  }
  transition = matrix(NA, nrow = 50, ncol = 4)
  
  # Iterate to find constants and probabilities
  for (t in 2:50){
    if(t %in% k){tau_est = tau_local}
    else{tau_est = tau_local*100}
    p_est = marg_prob[t]
    norm_const[t] = 1/(dnorm(y[t],0,tau_est)*p_est*forward_prob[t-1,1] + (1-p_est)*dnorm(y[t],0,tau_est)*forward_prob[t-1,2] + dnorm(y[t],1,tau_est)*(1-p_est)*forward_prob[t-1,1] + dnorm(y[t],1,tau_est)*p_est*forward_prob[t-1,2])
    transition[t,1] = p_est*dnorm(y[t],0,tau_est)*forward_prob[t-1,1]*norm_const[t]
    transition[t,2] = (1-p_est)*dnorm(y[t],1,tau_est)*forward_prob[t-1,1]*norm_const[t]
    transition[t,3] = (1-p_est)*dnorm(y[t],0,tau_est)*forward_prob[t-1,2]*norm_const[t]
    transition[t,4] = p_est*dnorm(y[t],1,tau_est)*forward_prob[t-1,2]*norm_const[t]
    forward_prob[t,1] = transition[t,1] + transition[t,3]
    forward_prob[t,2] = transition[t,2] + transition[t,4]
  }
  #print(forward_prob)
  # Initiate
  backward_prob = matrix(0,nrow = 50, ncol = 2)
  backward_prob[50,1] = forward_prob[50,1]
  backward_prob[50,2] = forward_prob[50,2]
  back_transition = matrix(NA, ncol = 4, nrow = 50)
  propagation = matrix(NA, ncol = 4, nrow = 50)
  
  
  for (t in 50:2){
    back_transition[t,1] = transition[t,1]/forward_prob[t,1]*backward_prob[t,1]
    back_transition[t,2] = transition[t,2]/forward_prob[t,2]*backward_prob[t,2]
    back_transition[t,3] = transition[t,3]/forward_prob[t,1]*backward_prob[t,1]
    back_transition[t,4] = transition[t,4]/forward_prob[t,2]*backward_prob[t,2]
    backward_prob[t-1,1] = back_transition[t,1] + back_transition[t,2]
    backward_prob[t-1,2] = back_transition[t,3] + back_transition[t,4]
    propagation[t,1] = back_transition[t,1]/backward_prob[t-1,1]
    propagation[t,2] = back_transition[t,2]/backward_prob[t-1,1]
    propagation[t,3] = back_transition[t,3]/backward_prob[t-1,2]
    propagation[t,4] = back_transition[t,4]/backward_prob[t-1,2]
  }
  return(backward_prob)
}
backward_prob = fb_algo(y, tau, k = c(20,30))

df2 = data.frame(prob_0 = backward_prob[,1], prob_1 = backward_prob[,2])
ggplot(df2, aes(x = 1:50, y = prob_1)) + geom_line() + labs(x = "i", y = "prob x_i = 1", title = "Computed marginal probabilities for x_i = 1") + theme_classic(base_size = 19)



# 1d)

monte_carlo_integration <- function(k, tau, n){
  pov = 0
  for (i in 1:n){
  y = rnorm(50, 0, 100*tau)
  y[k] = rnorm(1, x[k], tau)
  probabilities = fb_algo(y, tau, k)
  clean = -100000
  noclean = -5000*sum(probabilities[,2])
  print(noclean)
  pov = pov + max(clean,noclean)/n
  }
  return(pov)
}

pov = rep(NA, 50)
for(i in 1:50){
  print(i)
  pov[i] = monte_carlo_integration(i, 0.3, 10)}
