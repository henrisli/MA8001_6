library(ggplot2)
# 1a)
marg_prob = rep(NA, 50)
marg_prob[1] = 0.01
for(i in 2:50){
  marg_prob[i] = marg_prob[i-1] + 0.05*(1-marg_prob[i-1])
}

# 1c)
set.seed(260619)
tau = 0.3
x = rep(NA,50)
x[1] = rbinom(1, 1, marg_prob[1])
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
}

y = rnorm(50, x, 100*tau)
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
  p_est0 = 0.95
  p_est1 = 1
  
  # Iterate to find constants and probabilities
  for (t in 2:50){
    if(t %in% k){tau_est = tau_local}
    else{tau_est = tau_local*100}
    norm_const[t] = 1/(dnorm(y[t],0,tau_est)*p_est0*forward_prob[t-1,1] + (1-p_est1)*dnorm(y[t],0,tau_est)*forward_prob[t-1,2] + dnorm(y[t],1,tau_est)*(1-p_est0)*forward_prob[t-1,1] + dnorm(y[t],1,tau_est)*p_est1*forward_prob[t-1,2])
    transition[t,1] = p_est0*dnorm(y[t],0,tau_est)*forward_prob[t-1,1]*norm_const[t]
    transition[t,2] = (1-p_est0)*dnorm(y[t],1,tau_est)*forward_prob[t-1,1]*norm_const[t]
    transition[t,3] = (1-p_est1)*dnorm(y[t],0,tau_est)*forward_prob[t-1,2]*norm_const[t]
    transition[t,4] = p_est1*dnorm(y[t],1,tau_est)*forward_prob[t-1,2]*norm_const[t]
    forward_prob[t,1] = transition[t,1] + transition[t,3]
    forward_prob[t,2] = transition[t,2] + transition[t,4]
  }
  
  # Initiate
  backward_prob = matrix(0,nrow = 50, ncol = 2)
  backward_prob[50,1] = forward_prob[50,1]
  backward_prob[50,2] = forward_prob[50,2]
  back_transition = matrix(NA, ncol = 4, nrow = 50)
  #propagation = matrix(NA, ncol = 4, nrow = 50)
  
  
  for (t in 50:2){
    back_transition[t,1] = transition[t,1]/forward_prob[t,1]*backward_prob[t,1]
    back_transition[t,2] = transition[t,2]/forward_prob[t,2]*backward_prob[t,2]
    back_transition[t,3] = transition[t,3]/forward_prob[t,1]*backward_prob[t,1]
    back_transition[t,4] = transition[t,4]/forward_prob[t,2]*backward_prob[t,2]
    backward_prob[t-1,1] = back_transition[t,1] + back_transition[t,2]
    backward_prob[t-1,2] = back_transition[t,3] + back_transition[t,4]
    #propagation[t,1] = back_transition[t,1]/backward_prob[t-1,1]
    #propagation[t,2] = back_transition[t,2]/backward_prob[t-1,1]
    #propagation[t,3] = back_transition[t,3]/backward_prob[t-1,2]
    #propagation[t,4] = back_transition[t,4]/backward_prob[t-1,2]
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
    
    y = rnorm(50, 1/2, 100*tau)
    
    y[k] = rnorm(1, rbinom(1,1,marg_prob[k]), tau)
    
    probabilities = fb_algo(y, tau, k)
    
    clean = -100000
    noclean = -5000*sum(probabilities[,2])
    
    pov = pov + max(clean,noclean)/n
  }
  return(pov)
}

pov = rep(NA, 50)
for(i in 14:50){
  print(i)
  pov[i] = monte_carlo_integration(i, 0.3, 20000)}
write.table(pov, "C://Users//henri//Documents//GitHub//MA8001_6//pov_d.csv")
pov = as.vector(read.table("C://Users//henri//Documents//GitHub//MA8001_6//pov_d.csv"))
voi = pov + 100000
voi = voi[,1]
plot(voi, xlab = "k", ylab = "VoI")



plot(voi, xlab = "k", ylab = "VoI")

# 1e)

monte_carlo_integration_e <- function(k, tau, n){
  pov = 0
  for (i in 1:n){
    
    y = rnorm(50, 1/2, 100*tau)
    x = rep(NA, 50)
    x[k[1]] = rbinom(1, 1, marg_prob[k[1]])
    if(min(k)!=max(k)){
    for (i in (k[1]+1):k[2]){
      if(x[i-1] == 1){x[i] =1}
      else{
        if(runif(1) < 0.05){
          x[i] = 1
        }
        else{
          x[i] = 0
        }
      }
    }
    }
    if(min(k) == max(k)){y[k[1]] = rnorm(1,x[k[1]], tau/sqrt(2))}
    else{
    y[min(k)] = rnorm(1, x[min(k)], tau)
    y[max(k)] = rnorm(1, x[max(k)], tau)}

    probabilities = fb_algo(y, tau, k)
    
    clean = -100000
    noclean = -5000*sum(probabilities[,2])
    
    pov = pov + max(clean,noclean)/n
  }
  return(pov)
}

# positions = matrix(NA, ncol = 2, nrow = 820)
# iter = 1
# for (i in 10:49){
#   for (j in (i+1):50){
#     print(iter)
#     positions[iter,] = c(i,j)
#     iter = iter + 1
#   }}
# pov_e2 = apply(positions, 1, monte_carlo_integration_e, tau = 0.3, n = 100)

pov_e = matrix(NA, ncol = 50, nrow = 50)
# Diagonal elements
for (i in 1:50){
  k = c(i,i)
  pov_e[i,i] = monte_carlo_integration_e(k, 0.3, 100000)
  print(i)
}
# Off diagonal
for (i in 21:50){
  for (j in i:50){
    if (j==i){next}
    k = c(i,j)
    pov_e[i,j] = monte_carlo_integration_e(k, 0.3, 10000)
    print(j)
    pov_e[j,i] = pov_e[i,j]
  }
}

voi_e = pov_e+100000
#for (i in 1:50){voi_e[i,i] = voi[i]}
image.plot(x = 1:50, y = 1:50, z = voi_e, xlab = "i", ylab = "j")
write.table(pov_e, "C://Users//henri//Documents//GitHub//MA8001_6//pov_e.csv")
pov_e = as.matrix(read.table("C://Users//henri//Documents//GitHub//MA8001_6//pov_e.csv"))

max(voi_e)
which(voi_e == max(voi_e), arr.ind = TRUE)

max(voi_e[upper.tri(voi_e, diag = F)])
which(voi_e[upper.tri(voi_e, diag = F)] == max(voi_e[upper.tri(voi_e, diag = F)]), arr.ind = TRUE)

voi_e[20,30]

# 1f)

monte_carlo_integration_f <- function(){
  
  x = rep(NA,50)
  x[1] = rbinom(1, 1, marg_prob[1])
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
  }
  y = rnorm(50, x, 1)
  
  probabilities = fb_algo(y, 1, 1:50)
  
  clean = -100000
  noclean = -5000*sum(probabilities[,2])
  
  pov = max(clean,noclean)
  
  return(pov)
}

pov_f = rep(NA, 50000)
for(i in 1:50000){
  if(i%%1000 == 0){print(i)}
  pov_f[i] = monte_carlo_integration_f()}

plot(x = 1:1000,y=pov_f[1:1000]+100000, ylab = "VoI", xlab = "Monte Carlo Iteration") + abline(a=mean(pov_f+100000), b = 0, col = "red") + axis(2)
