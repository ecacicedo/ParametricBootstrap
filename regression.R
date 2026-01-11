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

# Jeffreys prior
log_prior <- function(beta, sigma) {
  -log(sigma)
}

# Likelihood ratio term R(theta)
logLik_theta <- function(beta, sigma, y) {
  mu <- X %*% beta
  -n * log(sigma) - sum((y - mu)^2) / (2 * sigma^2)
}

logR <- numeric(B)

for (b in 1:B) {
  logR[b] <- logLik_theta(beta_boot[b, ], sigma_boot[b], y) -
    logLik_theta(beta_hat, sigma_hat, y_star = NULL)
}

# Importance weights
logw <- sapply(1:B, function(b) {
  log_prior(beta_boot[b, ], sigma_boot[b]) + logR[b]
})

w <- exp(logw - max(logw))
w <- w / sum(w)

#  Posterior summaries for beta_1
post_mean <- sum(w * beta_boot[, 2])
post_ci <- quantile(beta_boot[, 2], probs = c(0.025, 0.975), weights = w)

list(mean = post_mean, CI = post_ci)





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