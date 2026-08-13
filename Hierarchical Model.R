#########################################################################
# Parametric bootstrap importance sampling Generalized Bayesian posterior
# Hierarchical model example
#########################################################################

# Check necessary packages

library(mnormt)
library(LaplacesDemon)

############################################################
# Simulated data process
# 
# y_i|theta_i, n_i ~ binomial (n_i, theta_i)
# theta_i|alpha,beta ~ beta(alpha, beta)
# n_i ~ uniform(200, 400)
# alpha ~ exponential(0.01)
# beta ~ exponential(0.01)
#
############################################################

k <- 20
n <- runif(k, min = 200, max = 400)
n <- round(n)
alpha <- rexp(1, rate=0.01)
beta <- rexp(1, rate=0.01)

theta <- rep(NA, k)
for(i in 1:k){
  theta[i] <- rbeta(n[i],alpha,beta)
}

y <- rbinom(20,n,theta)

theta_hat <- y/n

