# Parametric Bootstrap Importance Sampling ####
# Generalized Bayesian posterior
#
# True DGP:      GARCH(1,1)
# Working model: ARCH(1)
# Parameter vector: theta = (mu, omega, alpha)


library(mnormt)
library(LaplacesDemon)


# True GARCH(1,1) DGP ####

set.seed(123)

n <- 500

mu_true <- 0.0005

omega_garch_true <- 5e-6
alpha_garch_true <- 0.08
beta_garch_true <- 0.90

stopifnot(alpha_garch_true + beta_garch_true < 1)

# Simulate GARCH(1,1) ####
# True DGP: h_t = omega + alpha * (r_{t-1} - mu)^2 + beta * h_{t-1}

simulate_garch11 <- function(n, mu, omega, alpha, beta) {
  
  r <- numeric(n)
  h <- numeric(n)
  
  # Unconditional variance
  
  h[1] <- omega / (1 - alpha - beta)
  r[1] <- mu + sqrt(h[1]) * rnorm(1)
  
  for (t in 2:n) {
    h[t] <- omega + alpha * (r[t - 1] - mu)^2 + beta * h[t - 1]
    r[t] <- mu + sqrt(h[t]) *rnorm(1)
  }
  
  list(returns = r, variance = h)
}

garch_data <- simulate_garch11(
    n = n,
    mu = mu_true,
    omega = omega_garch_true,
    alpha = alpha_garch_true,
    beta = beta_garch_true
    )

r <- garch_data$returns

h_true <- garch_data$variance

# ARCH(1) conditional variance ####
# Working model: h_t = omega + alpha * (r_{t-1} - mu)^2

arch_variance <- function(theta, r_data) {
  
  mu <- theta[1]
  omega <- theta[2]
  alpha <- theta[3]
  
  # Parameter constraints
  
  if (!is.finite(mu) || !is.finite(omega) || !is.finite(alpha) ||
    omega <= 0 || alpha <= 0 || alpha >= 1) {
    return(rep(NA_real_, length(r_data)))
    }
  
  T_obs <- length(r_data)
  h <- numeric(T_obs)
  
  # Stationary ARCH initialization
  
  h[1] <- omega / (1 - alpha)
  
  if (!is.finite(h[1]) || h[1] <= 0) {return(rep(NA_real_, T_obs))}
  
  # ARCH recursion
  
  for (t in 2:T_obs) {h[t] <- omega + alpha *(r_data[t - 1] - mu)^2}
  
  h
}

# Gaussian ARCH loss ####
# L(theta; r) = 1/2 sum_t [log(h_t) + (r_t - mu)^2 / h_t]

arch_loss <- function(theta, r_data) {
  
  mu <- theta[1]
  omega <- theta[2]
  alpha <- theta[3]
  
  # Parameter constraints
  
  if (!is.finite(mu) || !is.finite(omega) || !is.finite(alpha) ||
    omega <= 0 || alpha <= 0 || alpha >= 1) {return(1e100)}
  
  h <- arch_variance(theta = theta, r_data = r_data)
  
  if (any(!is.finite(h)) || any(h <= 0)) {return(1e100)}
  
  residuals <- r_data - mu
  
  loss <- 0.5 * sum(log(h) + residuals^2 / h)
  
  if (!is.finite(loss)) {return(1e100)}
  
  loss
}

# ARCH empirical risk minimizer ####

fit_arch_erm <- function(r_data, start = NULL) {
  
  sample_variance <- var(r_data)
  
  if (is.null(start)) {
    start <- c(mu = mean(r_data), 
               omega = max(0.5 * sample_variance, 1e-8), alpha = 0.20)
    }
  
  # Make sure starting values satisfy constraints
  
  start[2] <- max(start[2], 1e-10)
  
  start[3] <- min(max(start[3], 1e-5), 1 - 1e-5)
  
  # Bounded optimization
  
  fit <- optim(par = start, fn = arch_loss, r_data = r_data,
      method = "L-BFGS-B", lower = c(-Inf, 1e-10, 1e-5),
      upper = c(Inf, Inf, 1 - 1e-5), control = list(maxit = 2000 ))
  
  names(fit$par) <- c("mu", "omega", "alpha")
  
  list(par = fit$par, loss = fit$value, convergence = fit$convergence)}

# Fit working ARCH(1) ####

arch_hat <- fit_arch_erm(r_data = r)

theta_hat <- arch_hat$par

# Priors ####
# mu ~ N(0, 0.02^2)
# omega ~ Lognormal(log(1e-4), 1.5^2)
# alpha ~ Beta(2, 5)

prior_mu_sd <- 0.02

prior_logomega_mean <- log(1e-4)

prior_logomega_sd <- 1.5

prior_alpha_a <- 2
prior_alpha_b <- 5

log_prior <- function(theta) {
  
  mu <- theta[1]
  omega <- theta[2]
  alpha <- theta[3]
  
  # Support restrictions

  if (!is.finite(mu) || !is.finite(omega) || !is.finite(alpha) ||
    omega <= 0 || alpha <= 0 || alpha >= 1) {return(-Inf)}
  
  dnorm(mu, mean = 0, sd = prior_mu_sd, log = TRUE) + 
    dlnorm(omega, meanlog = prior_logomega_mean, 
           sdlog = prior_logomega_sd, log = TRUE) +
    dbeta(alpha,mshape1 = prior_alpha_a, shape2 = prior_alpha_b, log = TRUE)
}

# Generalized Bayesian learning rate ####

eta <- 0.50

# Generalized Bayesian posterior target ####
# pi_G(theta | r) \propto pi(theta) exp{-eta L(theta; r)}


log_target <- function(theta, r_data = r) {
  
  lp <- log_prior(theta)
  
  if (!is.finite(lp)) {return(-Inf)}
  
  loss <- arch_loss(theta = theta, r_data = r_data)
  
  if (!is.finite(loss) || loss >= 1e99) {return(-Inf)}
  
  lp - eta * loss
}

# Simulate from fitted ARCH(1) ####
# Used for parametric bootstrap

simulate_arch1 <- function(n, theta) {
  
  mu <- theta[1]
  omega <- theta[2]
  alpha <- theta[3]
  
  if (omega <= 0 || alpha <= 0 || alpha >= 1) {
    
    stop("Invalid ARCH parameters.")}
  
  r_star <- numeric(n)
  h_star <- numeric(n)
  
  # Stationary initialization

  h_star[1] <- omega / (1 - alpha)
  r_star[1] <- mu + sqrt(h_star[1]) * rnorm(1)
  
  # Simulation

  for (t in 2:n) {
    h_star[t] <- omega + alpha * (r_star[t - 1] - mu)^2
    r_star[t] <- mu + sqrt(h_star[t]) * rnorm(1)
    }
  
  list(returns = r_star, variance = h_star)}

# Positive definite covariance regularization

make_positive_definite <- function(Sigma, min_eigenvalue = 1e-12) {
  
  Sigma <-(Sigma + t(Sigma)) /2
  
  eig <- eigen(Sigma, symmetric = TRUE)
  
  eig$values <-pmax(eig$values, min_eigenvalue)
  
  Sigma_pd <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
  
  (Sigma_pd + t(Sigma_pd)) /2
}

# Weighted covariance

weighted_covariance <- function(draws, weights) {
  
  weights <- weights / sum(weights)
  
  weighted_mean <- colSums(draws * weights)
  
  centered <-sweep(draws, MARGIN = 2, STATS = weighted_mean, FUN = "-")
  
  covariance <-crossprod(centered *sqrt(weights))
  
  covariance / (1 - sum(weights^2))
  }

# Parametric Bootstrap Importance Sampling

bootstrap_is <- function(B_boot = 3000, B_is = 10000, proposal_inflation = 1.25) {
  start_time <- proc.time()[3]
  
  # ERM under misspecified ARCH working model

  observed_fit <- fit_arch_erm(r_data = r)
  
  theta_hat <- observed_fit$par
  
  # Parametric bootstrap distribution of the ARCH ERM

  bootstrap_estimators <- matrix(NA_real_, nrow = B_boot, ncol = 3)
  
  colnames(bootstrap_estimators) <- c("mu", "omega", "alpha")
  successful <- 0
  attempts <-0
  max_attempts <- 2 * B_boot
  
  while (successful < B_boot && attempts < max_attempts) {
    
    attempts <- attempts + 1
    
    # Generate bootstrap data from fitted ARCH model

    simulated_data <- simulate_arch1(n = n, theta = theta_hat)
    
    # Refit ARCH model

    fit_star <- try(fit_arch_erm(r_data = simulated_data$returns, 
                                 start = theta_hat), silent = TRUE)
    
    if (!inherits(fit_star, "try-error") && fit_star$convergence == 0 &&
        all(is.finite(fit_star$par)) && fit_star$par[2] > 0 && 
        fit_star$par[3] > 0 && fit_star$par[3] < 1) {
      
      successful <- successful + 1
      
      bootstrap_estimators[successful,] <- fit_star$par}
    }
  
  if (successful < 100) {stop("Too few successful bootstrap fits.")}
  
  bootstrap_estimators <- bootstrap_estimators[seq_len(successful), ,drop = FALSE]
  
  # Fit multivariate Gaussian proposal
  # q(theta) = N(theta ; mean_boot, Sigma_boot)

  proposal_mean <-colMeans(bootstrap_estimators)
  proposal_covariance <- cov(bootstrap_estimators)
  proposal_covariance <- make_positive_definite(proposal_covariance)
  
  proposal_covariance <- proposal_inflation^2 * proposal_covariance
  
  # Draw importance samples from fitted proposal

  proposal_draws <- mnormt::rmnorm(n = B_is, mean = proposal_mean, 
                                   varcov = proposal_covariance)
  
  colnames(proposal_draws) <- c("mu", "omega", "alpha")
  
  # Evaluate posterior target

  log_target_values <- apply(proposal_draws, MARGIN = 1, FUN = log_target)
  
  # Proposal density

  log_proposal_values <- mnormt::dmnorm(x = proposal_draws, mean = proposal_mean, 
                                        varcov = proposal_covariance, log = TRUE)
  
  # Importance weights

  log_weights <- log_target_values - log_proposal_values
  
  finite_weights <- is.finite(log_weights)
  
  if (!any(finite_weights)) {stop("All importance weights are non-finite.")}
  
  log_weights[!finite_weights] <- -Inf
  
  # Numerically stable normalization

  max_log_weight <- max(log_weights)
  raw_weights <- exp(log_weights - max_log_weight)
  weights <- raw_weights / sum(raw_weights)
  
  # Importance sampling ESS

  ess <- 1 /sum(weights^2)
  
  # Weighted posterior moments

  post_mean <- colSums(proposal_draws * weights)
  
  post_covariance <- weighted_covariance(draws = proposal_draws, weights = weights)
  
  post_summary <- data.frame(parameter = colnames(proposal_draws), mean = post_mean, 
                             sd = sqrt(diag(post_covariance)), row.names = NULL)
  
  # Fraction of proposal draws satisfying ARCH constraints

  valid_proposal <- proposal_draws[, "omega"] > 0 & 
    proposal_draws[, "alpha"] > 0 & 
    proposal_draws[, "alpha"] < 1
  
  valid_proposal_fraction <- mean(valid_proposal)
  
  elapsed_time <- proc.time()[3] - start_time
  
  list(method = "Parametric bootstrap importance sampling", 
       theta_hat = theta_hat, 
       bootstrap_estimators = bootstrap_estimators, 
       proposal_mean = proposal_mean, 
       proposal_covariance = proposal_covariance, 
       draws = proposal_draws, 
       weights = weights, 
       summary = post_summary, 
       covariance = post_covariance, 
       ess = ess, 
       relative_ess = ess / B_is, 
       maximum_weight = max(weights), 
       valid_proposal_fraction = valid_proposal_fraction, 
       elapsed = elapsed_time, 
       ess_per_second = ess / elapsed_time, 
       successful_bootstraps = successful, 
       attempted_bootstraps = attempts)
  }

# Random Walk Metropolis Hastings ####

rwmh <- function(n_iter = 40000, burnin = 10000) {
  start_time <- proc.time()[3]
  p <- 3
  
  # ARCH ERM as starting value

  start_fit <-fit_arch_erm(r_data = r)
  
  start_theta <- start_fit$par
  
  # Negative log posterior

  negative_log_target <- function(theta) {value <- log_target(theta)
      
      if (!is.finite(value)) {return(1e100)}

      - value
  }
  
  
  # Posterior mode using direct constraints

  mode_fit <- optim(par = start_theta, fn = negative_log_target, method = "L-BFGS-B", 
                    lower = c(-Inf, 1e-10, e-5), upper = c(Inf, Inf, 1 - 1e-5), 
                    control = list(maxit = 2000))
  
  posterior_mode <- mode_fit$par
  
  names(posterior_mode) <- c("mu", "omega", "alpha")
  
  # Numerical Hessian around posterior mode

  Hessian <- optimHess(par = posterior_mode, fn = negative_log_target)
  
  Hessian <- (Hessian + t(Hessian)) /2
  
  # Ensure positive definite Hessian

  eig <- eigen(Hessian,symmetric = TRUE)
  
  eig$values <- pmax(eig$values, 1e-8)
  
  Hessian_pd <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
  
  local_covariance <- solve(Hessian_pd)
  
  # Random walk proposal covariance

  proposal_scale <- 2.38^2 / p
  
  proposal_covariance <- proposal_scale * local_covariance
  
  proposal_covariance <- make_positive_definite(proposal_covariance)
  
  # MCMC storage

  draws <- matrix(NA_real_, nrow = n_iter, ncol = p)
  
  colnames(draws) <- c("mu", "omega", "alpha")
  
  current_theta <- posterior_mode
  
  current_log_target <- log_target(current_theta)
  
  accepted <-0
  
  # MCMC

  
  for (i in seq_len(n_iter)) {
    
    proposed_theta <-
      as.vector(mnormt::rmnorm(n = 1, mean = current_theta, 
                               varcov = proposal_covariance))
    
    proposed_log_target <- log_target(proposed_theta)
    
    # Invalid omega or alpha automatically give
    # proposed_log_target = -Inf

    if (is.finite(proposed_log_target)) {
      log_acceptance_ratio <- proposed_log_target - current_log_target
      if (log(runif(1)) < min(0,log_acceptance_ratio)) {
        current_theta <- proposed_theta
        current_log_target <- proposed_log_target
        accepted <- accepted + 1}
      }
    
    draws[i, ] <-current_theta
    }
  
  # Remove burnin
  
  draws_keep <- draws[(burnin + 1):n_iter,, drop = FALSE]
  
  # ESS

  ess_parameter <- apply(draws_keep, MARGIN = 2, FUN = LaplacesDemon::ESS)
  
  overall_ess <- min(ess_parameter)
  
  # Posterior summary
  
  post_summary <- data.frame(parameter = colnames(draws_keep), 
                             mean = colMeans(draws_keep), 
                             sd = apply(draws_keep, 2, sd), 
                             ess = ess_parameter, row.names = NULL)
  
  elapsed_time <- proc.time()[3] - start_time
  
  list(method = "Random walk Metropolis Hastings",
       draws = draws_keep, 
       summary = post_summary, 
       posterior_mode = posterior_mode, 
       proposal_covariance = proposal_covariance, 
       ess_parameter = ess_parameter,
       ess = overall_ess,
       acceptance_rate = accepted / n_iter,
       elapsed = elapsed_time,
       ess_per_second = overall_ess / elapsed_time)
}

# Run PBIS ####

set.seed(456)

bootstrap_result <-
  bootstrap_is(
    B_boot = 3000,
    B_is = 10000,
    proposal_inflation = 1.25
  )

# Run RWMH ####

set.seed(789)

mcmc_result <- rwmh(n_iter = 40000, burnin = 10000)

# Computational cost comparison ####

cost_table <- data.frame(method = c(bootstrap_result$method, mcmc_result$method), 
                         elapsed_seconds = c(bootstrap_result$elapsed, mcmc_result$elapsed), 
                         ESS = c(bootstrap_result$ess, mcmc_result$ess), 
                         ESS_per_second = c(bootstrap_result$ess_per_second, mcmc_result$ess_per_second))

comparison_metrics <- data.frame(
  metric =  c("Time difference: bootstrap minus MCMC", 
              "ESS per second difference: bootstrap minus MCMC", 
              "Relative ESS per second: bootstrap divided by MCMC"),
  value = c(bootstrap_result$elapsed - mcmc_result$elapsed, 
            bootstrap_result$ess_per_second - mcmc_result$ess_per_second, 
            bootstrap_result$ess_per_second / mcmc_result$ess_per_second))

# Posterior comparison ####

post_comparison <- data.frame(parameter = c("mu", "omega", "alpha"),
                              bootstrap_mean = bootstrap_result$summary$mean,
                              mcmc_mean = mcmc_result$summary$mean, 
                              absolute_mean_difference = abs(bootstrap_result$summary$mean - mcmc_result$summary$mean), 
                              bootstrap_sd = bootstrap_result$summary$sd, 
                              mcmc_sd = mcmc_result$summary$sd)

# Importance sampling diagnostics ####

weight_diagnostics <-data.frame(ESS = bootstrap_result$ess, 
                                relative_ESS = bootstrap_result$relative_ess, 
                                maximum_weight = bootstrap_result$maximum_weight, 
                                valid_proposal_fraction = bootstrap_result$valid_proposal_fraction, 
                                successful_bootstraps = bootstrap_result$successful_bootstraps, 
                                attempted_bootstraps = bootstrap_result$attempted_bootstraps)

# Plot simulated returns

plot(r, type = "l", xlab = "Time", ylab = "Return", 
     main = "Returns generated from GARCH(1,1)")

# Plot true GARCH conditional variance ####

plot(h_true, type = "l", xlab = "Time", ylab = expression(h[t]), 
     main = "True GARCH(1,1) conditional variance")

# Posterior density comparison ####

old_par <- par(no.readonly = TRUE)

par(mar = c(5, 5, 4, 10))

parameter_labels <- expression(mu, omega, alpha)

for (j in seq_len(3)) {
  
  # PBIS draws and weights

  bootstrap_draws_j <- bootstrap_result$draws[, j]
  
  # Remove invalid proposal draws before KDE

  valid <- is.finite(bootstrap_draws_j) & bootstrap_result$weights > 0
  
  bootstrap_draws_j <- bootstrap_draws_j[valid]
  bootstrap_weights_j <- bootstrap_result$weights[valid]
  bootstrap_weights_j <- bootstrap_weights_j /sum(bootstrap_weights_j)
  
  # MCMC draws

  mcmc_draws_j <- mcmc_result$draws[, j]
  
  # Common range

  plot_range <- range(bootstrap_draws_j, mcmc_draws_j,finite = TRUE)

  # Weighted PBIS posterior density

  bootstrap_density <- density(x =bootstrap_draws_j, 
                               weights =bootstrap_weights_j,
                               n = 500,
                               from = plot_range[1], 
                               to = plot_range[2])
  
  # RWMH posterior density

  mcmc_density <- density(x = mcmc_draws_j, 
                          n = 500, 
                          from = plot_range[1], 
                          to = plot_range[2])
  
  y_limit <- c(0, 1.05 * max(bootstrap_density$y,mcmc_density$y, na.rm = TRUE))
  
  # Plot

  plot(x = bootstrap_density$x, 
       y = bootstrap_density$y, 
       type = "l", 
       lwd = 2,
       col = "blue",
       xlab = parameter_labels[j],
       ylab = "Posterior density",
       main = bquote("Generalized Bayesian ARCH posterior for " *.(parameter_labels[[j]])),
       xlim = plot_range, ylim = y_limit)
  
  lines(x = mcmc_density$x, y = mcmc_density$y, lwd = 2, col = "red", lty = 2)
  
  # Only mu has a directly comparable true value

  if (j == 1) {abline(v = mu_true, col = "black", lwd = 2, lty =3)
    
    legend("topright", 
           inset = c(-1.15, 0), 
           legend = c("Parametric bootstrap IS", "Random walk MH", "True value"), 
           col = c("blue","red","black"),
           lwd = 2, 
           lty = c(1,2,3),
           bty = "n", 
           xpd = TRUE)} 
  else {legend("topright", 
               inset =c(-1.15, 0),
               legend = c("Parametric bootstrap IS", "Random walk MH"),
               col =c("blue","red"),
               lwd = 2,
               lty = c(1, 2), 
               bty = "n",
               xpd = TRUE)}
  }

par(old_par)

# Conditional variance at posterior means

pbis_mean <-setNames(bootstrap_result$summary$mean, bootstrap_result$summary$parameter)

mcmc_mean <-setNames(mcmc_result$summary$mean, mcmc_result$summary$parameter)

h_arch_pbis <- arch_variance(theta =pbis_mean, r_data =r)

h_arch_mcmc <-arch_variance(theta = mcmc_mean, r_data = r)

# Compare true GARCH variance and fitted ARCH variance

plot(h_true, type = "l", lwd = 2, xlab = "Time", ylab = expression(h[t]),
     main = "True GARCH variance versus fitted ARCH variance")

lines(h_arch_pbis, lwd = 2, col = "blue", lty = 2)

lines(h_arch_mcmc, lwd = 2, col = "red", lty = 3)


legend("topright",
       legend = c("True GARCH variance",
                  "ARCH posterior mean: PBIS",
                  "ARCH posterior mean: RWMH"),
       col = c("black","blue","red"),
       lwd = 2,
       lty = c(1, 2, 3), 
       bty ="n")

# Output ####

cat("\nTrue GARCH(1,1) DGP parameters:\n")

print(c(mu = mu_true, omega = omega_garch_true, 
        alpha = alpha_garch_true, beta = beta_garch_true))

cat("\nARCH(1) empirical risk minimizer:\n")

print(theta_hat)

cat("\nGeneralized Bayesian learning rate:\n")

print(eta)

cat("\nParametric bootstrap IS posterior summary:\n")

print(bootstrap_result$summary)

cat("\nRWMH posterior summary:\n")

print(mcmc_result$summary)

cat("\nRWMH acceptance rate:\n")

print(mcmc_result$acceptance_rate)

cat("\nImportance sampling diagnostics:\n")

print(weight_diagnostics)

cat("\nComputational cost comparison:\n")

print(cost_table)

cat("\nEfficiency comparison:\n")

print(comparison_metrics)

cat("\nPosterior comparison:\n")

print(post_comparison)