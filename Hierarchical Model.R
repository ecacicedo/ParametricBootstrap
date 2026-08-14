# Title ####
# Parametric bootstrap importance sampling 
# Generalized Bayesian posterior
# Hierarchical model example

# Check necessary packages ####

library("extraDistr")

i <- 20
j <- 3

# Simulated data process ####
# y_i,j|beta_j ~ normal(beta_j * X_ij, 1)
# beta_j ~ normal(mu_b, 0.15)
# mu_b ~ uniform(0.5,1)

set.seed(123)

X <- round(rnorm(i * j, mean = 50, sd = 10))

X1 <- X[1 : i]
X2 <- X[(i + 1) : (2 * i)]
X3 <- X[(2 * i + 1) : (3 * i)]

mu_b <- runif(1, min = 0.5, max = 1)

beta <- rnorm(j, mean = mu_b, sd = 0.15)

y1 <- X1 * beta[1] + rnorm(i, mean = 0, sd = 1)
y2 <- X2 * beta[2] + rnorm(i, mean = 0, sd = 1)
y3 <- X3 * beta[3] + rnorm(i, mean = 0, sd = 1)
y1
y2
y3

