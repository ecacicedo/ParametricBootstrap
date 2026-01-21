set.seed(1)

# Simulated data
n <- 100
x <- runif(n)
X <- cbind(1, x)
beta_true <- c(1, 2)
sigma_true <- 1

y <- X %*% beta_true + rnorm(n, sd = sigma_true)

# Fit model (MLE)
fit <- lm(y ~ x)
beta_hat <- coef(fit)
sigma_hat <- sqrt(sum(residuals(fit)^2) / n)

# Parametric bootstrap
B <- 10000
beta_boot <- matrix(NA, B, length(beta_hat))
sigma_boot <- numeric(B)

for (b in 1:B) {
  y_star <- X %*% beta_hat + rnorm(n, sd = sigma_hat)
  fit_star <- lm(y_star ~ x)
  beta_boot[b, ] <- coef(fit_star)
  sigma_boot[b] <- sqrt(sum(residuals(fit_star)^2) / n)
}

# Jeffrey's prior
jeff_prior <- function(sigma) sigma^(-2)

# Likelihood ratio term R(theta)

R <- numeric(B)

mu_boot = X %*% beta_boot[b, ]
mu_hat = X %*% beta_hat

for (b in 1:B) {
  R[b] <- dnorm(mu_hat,mean=mu_boot, sd=sigma_boot[b])/dnorm(mu_boot, mean=mu_hat, sd=sigma_hat)
}

w <- jeff_prior(sigma_boot) * R
w <- w / sum(w)

post_mean <- sum(w * beta_boot[, 2])

post_ci <- quantile(beta_boot[, 2], probs = c(0.025, 0.975), weights = w)

# Posterior density
dens_boot <- density(beta_boot[, 2])

bw0 <- density(beta_boot[, 2])$bw
dens_post <- density(beta_boot[, 2], weights = w, bw = bw0)

# Plot
plot(dens_boot, col="grey", main="Raw bootstrap density v LR posterior with Jeffrey's prior")
lines(dens_post, col="black")
legend("topright", legend=c("Bootstrap","Jeffrey's posterior"),
       col=c("grey","black"), lty=1)





# Posterior using MCMC

# Squared Loss
sq_loss <- function(r) {
          0.5 * r^2
}

# Cumulative Loss
loss_cum <- function(beta, X, y) {
  r <- as.vector(y - X %*% beta)
  sum(sq_loss(r))
}

# Log Target distribution
log_target <- function(beta, sigma, eta, X, y) {
  log_prior(beta, sigma) - eta * loss_cum(beta, X, y)
}

# Random Walk Metropolis-Hastings sampler

mh_sampler <- function(
    X, 
    y, 
    eta,
    beta_init,
    cov_init,
    prop_cov,
    n_iter,
    burnin
) {
  p <- length(beta_init)
  samples <- matrix(NA, n_iter, p)
  weights <- matrix(NA, n_iter, 1)
  
  beta_curr <- beta_init
  log_curr  <- log_target(beta_curr, cov_init, eta, X, y)
  
  accept <- 0
  n_total <- n_iter + burnin
  
  for (t in 1:n_total) {
    
    beta_prop <- MASS::mvrnorm(1, beta_curr, prop_cov)
    log_prop  <- log_target(beta_prop, cov_init, eta, X, y)
    
    accepted <- FALSE
    if (log(runif(1)) < log_prop - log_curr) {
      beta_curr <- beta_prop
      log_curr  <- log_prop
      accepted <- TRUE
    }
    
    if (t > burnin) {
      samples[t - burnin, ] <- beta_curr
      weights[t - burnin, ] <- log_curr
      if (accepted) accept <- accept + 1
    }
  }
  
  list(
    weights = weights/sum(weights),
    samples = samples,
    accept_ratio = accept / n_iter,
    mean = colMeans(samples),
    cov = cov(samples)
  )
}