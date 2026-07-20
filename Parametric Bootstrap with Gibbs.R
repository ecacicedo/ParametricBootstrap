############################################################
# Parametric bootstrap importance sampling Gibbs posterior 
# Simple linear regression example
############################################################

library(mnormt)
library(LaplacesDemon)
library(topolow)

############################################################
# Simulated data
############################################################

set.seed(123)

n <- 100
x <- runif(n)
X <- cbind(1, x)

colnames(X) <- c("beta0", "beta1")

beta_true <- c(1, 2)
sigma_true <- 1

y <- as.vector(X %*% beta_true + rnorm(n, sd = sigma_true))

p <- ncol(X)

############################################################
# Loss, learning rate and prior
############################################################

loss_total <- function(beta, y_data) {
  residuals <- as.vector(y_data - X %*% beta)
  0.5 * sum(residuals^2)
}

eta <- 1

prior_sd <- 10

log_prior <- function(beta) {-0.5 * sum(beta^2 / prior_sd^2)}

log_target <- function(beta) {log_prior(beta) - eta * loss_total(beta, y)}

############################################################
# Parametric bootstrap importance sampler
############################################################

bootstrap_is <- function(B = 5000) {
  
  start_time <- proc.time()[3]
  
  ##########################################################
  # Empirical risk minimiser
  ##########################################################
  
  beta_hat <- as.vector(solve(t(X) %*% y))
  
  ##########################################################
  # Parametric bootstrap distribution
  #
  # Y_i | beta ~ N(x_i' beta, 1 / eta)
  ##########################################################
  
  sigma_boot <- sqrt(1 / eta)
  
  beta_boot <- matrix(NA_real_, nrow = B, ncol = p)
  
  for (b in seq_len(B)) {
    y_star <- as.vector(X %*% beta_hat + rnorm(n, sd = sigma_boot))
    
    beta_boot[b, ] <- solve(t(X) %*% X,t(X) %*% y_star)}
  
  colnames(beta_boot) <- colnames(X)
  
  ##########################################################
  # Normal approximation to bootstrap proposal
  ##########################################################
  
  Sigma_boot <- cov(beta_boot)
  
  Sigma_boot <- Sigma_boot + 1e-10 * diag(p) # for num stability
  
  ##########################################################
  # Log importance weights
  ##########################################################
  
  log_target_values <- apply(beta_boot, 1, log_target)
  
  log_proposal_values <- mnormt::dmnorm(
    x = beta_boot,
    mean = beta_hat,
    varcov = Sigma_boot,
    log = TRUE
  )
  
  log_weights <- log_target_values - log_proposal_values
  
  # Log sum exp stabilisation
  log_weights <- log_weights - max(log_weights)
  
  raw_weights <- exp(log_weights)
  
  weights <- raw_weights / sum(raw_weights)
  
  ##########################################################
  # Importance sampling ESS
  ##########################################################
  
  ess <- 1 / sum(weights^2)
  
  ##########################################################
  # Weighted posterior summaries
  ##########################################################
  
  posterior_mean <- colSums(beta_boot * weights)
  
  centered_draws <- sweep(beta_boot, 2, posterior_mean, "-")
  
  posterior_covariance <- t(centered_draws) %*% (centered_draws * weights)
  
  posterior_summary <- data.frame(
    parameter = colnames(beta_boot),
    mean = posterior_mean,
    sd = sqrt(diag(posterior_covariance))
  )
  
  elapsed_time <- proc.time()[3] - start_time
  
  list(
    method = "Parametric bootstrap IS",
    beta_hat = beta_hat,
    draws = beta_boot,
    weights = weights,
    summary = posterior_summary,
    ess = ess,
    elapsed = elapsed_time,
    ess_per_second = ess / elapsed_time
  )
}

############################################################
# Random Walk Metropolis Hastings
############################################################

rwmh <- function(n_iter = 30000, burnin = 5000) {
  
  start_time <- proc.time()[3]
  
  beta_hat <- as.vector(solve(t(X) %*% X, t(X) %*% y))
  
  ##########################################################
  # Local Gibbs posterior covariance
  #
  # Neg log posterior Hessian: eta X'X + prior precision
  #
  ##########################################################
  
  prior_precision <- diag(1 / prior_sd^2, p)
  
  local_covariance <- solve(eta * t(X) %*% X + prior_precision)
  
  proposal_scale <- 2.38^2 / p # Gelman's paper
  
  proposal_covariance <- proposal_scale * local_covariance
  
  draws <- matrix(NA_real_, nrow = n_iter, ncol = p)
  
  current_beta <- beta_hat
  current_log_target <- log_target(current_beta)
  
  accepted <- 0
  
  for (iteration in seq_len(n_iter)) {
    
    proposed_beta <- as.vector(
      mnormt::rmnorm(
        n = 1,
        mean = current_beta,
        varcov = proposal_covariance
      )
    )
    
    proposed_log_target <- log_target(proposed_beta)
    
    log_acceptance_ratio <- proposed_log_target - current_log_target
    
    if (log(runif(1)) < min(0, log_acceptance_ratio)) {
      current_beta <- proposed_beta
      current_log_target <- proposed_log_target
      accepted <- accepted + 1
    }
    
    draws[iteration, ] <- current_beta
  }
  
  draws_keep <- draws[(burnin + 1):n_iter, , drop = FALSE]
  
  colnames(draws_keep) <- colnames(X)
  
  ess_by_parameter <- apply(draws_keep, 2, LaplacesDemon::ESS)
  
  posterior_summary <- data.frame(
    parameter = colnames(draws_keep),
    mean = colMeans(draws_keep),
    sd = apply(draws_keep, 2, sd),
    ess = ess_by_parameter,
    row.names = NULL
  )
  
  elapsed_time <- proc.time()[3] - start_time
  
  # Use the smallest parameter ESS as a conservative
  # overall MCMC efficiency measure.
  overall_ess <- min(ess_by_parameter)
  
  list(
    method = "Random walk Metropolis Hastings",
    draws = draws_keep,
    summary = posterior_summary,
    ess_by_parameter = ess_by_parameter,
    ess = overall_ess,
    acceptance_rate = accepted / n_iter,
    elapsed = elapsed_time,
    ess_per_second = overall_ess / elapsed_time
  )
}

############################################################
# Run both methods
############################################################

set.seed(456)

bootstrap_result <- bootstrap_is(B = 5000)

set.seed(789)

mcmc_result <- rwmh(n_iter = 30000, burnin = 5000)

############################################################
# Computational cost comparison
############################################################

cost_table <- data.frame(
  method = c(bootstrap_result$method, 
             mcmc_result$method),
  elapsed_seconds = c(bootstrap_result$elapsed, 
                      mcmc_result$elapsed),
  ESS = c(bootstrap_result$ess, 
          mcmc_result$ess),
  ESS_per_second = c(bootstrap_result$ess_per_second, 
                     mcmc_result$ess_per_second))

comparison_metrics <- data.frame(
  metric = c(
    "Time difference: bootstrap minus MCMC",
    "ESS per second difference: bootstrap minus MCMC",
    "Relative ESS per second: bootstrap divided by MCMC"
  ),
  value = c(bootstrap_result$elapsed - mcmc_result$elapsed,
    bootstrap_result$ess_per_second - mcmc_result$ess_per_second,
    bootstrap_result$ess_per_second / mcmc_result$ess_per_secon)
)

############################################################
# Compare posterior means
############################################################

posterior_comparison <- data.frame(
  parameter = colnames(X),
  bootstrap_mean = bootstrap_result$summary$mean,
  mcmc_mean = mcmc_result$summary$mean,
  absolute_mean_difference = abs(
    bootstrap_result$summary$mean -
      mcmc_result$summary$mean
  ),
  bootstrap_sd = bootstrap_result$summary$sd,
  mcmc_sd = mcmc_result$summary$sd
)

############################################################
# Overlay posterior densities
############################################################

old_par <- par(no.readonly = TRUE)

# Extra right margin for the legend
par(mar = c(5, 5, 4, 10))

for (j in seq_len(p)) {
  
  bootstrap_draws_j <- bootstrap_result$draws[, j]
  mcmc_draws_j <- mcmc_result$draws[, j]
  
  plot_range <- range(bootstrap_draws_j, mcmc_draws_j, finite = TRUE)
  
  ##########################################################
  # Weighted posterior density estimates
  ##########################################################
  
  bootstrap_density <- topolow::weighted_kde(
    x = bootstrap_draws_j,
    weights = bootstrap_result$weights,
    n = 500,
    from = plot_range[1],
    to = plot_range[2]
  )
  
  mcmc_density <- topolow::weighted_kde(
    x = mcmc_draws_j,
    weights = rep(
      1 / length(mcmc_draws_j),
      length(mcmc_draws_j)
    ),
    n = 500,
    from = plot_range[1],
    to = plot_range[2]
  )
  
  ##########################################################
  # Common vertical axis
  ##########################################################
  
  y_limit <- c(0, 1.05 * max(bootstrap_density$y, mcmc_density$y, na.rm = TRUE))
  
  ##########################################################
  # Plot bootstrap IS density
  ##########################################################
  
  plot(
    x = bootstrap_density$x,
    y = bootstrap_density$y,
    type = "l",
    lwd = 2,
    col = "blue",
    xlab = colnames(X)[j],
    ylab = "Posterior density",
    main = paste(
      "Gibbs posterior for",
      colnames(X)[j]
    ),
    xlim = plot_range,
    ylim = y_limit
  )
  
  ##########################################################
  # Add MCMC density
  ##########################################################
  
  lines(x = mcmc_density$x, y = mcmc_density$y, lwd = 2, col = "red", lty = 2)
  
  ##########################################################
  # Add true parameter value
  ##########################################################
  
  abline(v = beta_true[j], col = "black", lwd = 2, lty = 3)
  
  ##########################################################
  # Legend outside the plot
  ##########################################################
  
  legend(
    "topright",
    inset = c(-1.15, 0),
    legend = c(
      "Parametric bootstrap IS",
      "Random walk MH",
      "True value"
    ),
    col = c(
      "blue",
      "red",
      "black"
    ),
    lwd = 2,
    lty = c(1, 2, 3),
    bty = "n",
    xpd = TRUE
  )
}

par(old_par)

############################################################
# Posterior summaries
############################################################

cat("\nParametric bootstrap IS summary:\n")
print(bootstrap_result$summary)

cat("\nMCMC summary:\n")
print(mcmc_result$summary)

cat("\nMCMC acceptance rate:\n")
print(mcmc_result$acceptance_rate)

cat("\nComputational cost comparison:\n")
print(cost_table)

cat("\nEfficiency differences:\n")
print(comparison_metrics)

cat("\nPosterior comparison:\n")
print(posterior_comparison)