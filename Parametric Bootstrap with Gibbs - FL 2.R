############################################################
# Bootstrap importance sampling Gibbs posterior
# Imbalanced binary classification with focal loss
############################################################

library(mnormt)
library(LaplacesDemon)
library(topolow)

set.seed(123)

############################################################
# Simulate imbalanced classification data
############################################################

n <- 500

x1 <- rnorm(n)
x2 <- rbinom(n, size = 1, prob = 0.35)

X <- cbind(intercept = 1, x1 = x1, x2 = x2)

beta_true <- c(intercept = -3.0, x1 = 1.25, x2 = 0.75)

# use 'plogis' function to calculate exp(x)/(1 + exp(x))
# true_prob <- plogis(as.vector(X%*%beta_true)) 
true_prob <- exp(as.vector(X%*%beta_true))/(1+exp(as.vector(X%*%beta_true)))

y <- rbinom(n = n, size = 1, prob = true_prob)

p <- ncol(X)

############################################################
# Focal loss settings
############################################################

# Higher alpha gives greater loss weight to positive cases
alpha_focal <- 0.75

# gamma > 0 downweights confidently classified observations
gamma_focal <- 2

eta <- 1

prior_sd <- 5

############################################################
# Focal loss
############################################################

#Check if really needed
clip_prob <- function(probability, epsilon = 1e-10) {
  pmin(pmax(probability, epsilon),1 - epsilon)
}

focal_loss_vector <- function(beta,X_data,y_data) {
  
  prob_one <- clip_prob(plogis(as.vector(X_data %*% beta)))
  
  prob_true_class <- ifelse(y_data == 1,prob_one,1 - prob_one)
  
  alpha_true_class <- ifelse(y_data == 1,alpha_focal,1 - alpha_focal)
  
  -alpha_true_class * (1 - prob_true_class)^gamma_focal*log(prob_true_class)
}

focal_loss_total <- function(beta,X_data,y_data) {
  sum(focal_loss_vector(beta = beta,X_data = X_data,y_data = y_data))
}

############################################################
# Gradient of focal loss
############################################################

focal_grad <- function(beta,X_data,y_data) {
  
  prob_one <- clip_prob(plogis(as.vector(X_data%*%beta)))
  
  derivative_xbeta <- numeric(length(y_data))
  
  positive <- y_data == 1
  negative <- !positive
  
  ##########################################################
  # Positive class derivative w/ respect to XBeta
  ##########################################################
  
  if (any(positive)) {
    
    prob_positive <- prob_one[positive]
    
    derivative_xbeta[positive] <- alpha_focal*(1 - prob_positive)^gamma_focal*
      (gamma_focal*prob_positive*log(prob_positive) - (1 - prob_positive))
  }
  
  ##########################################################
  # Negative class derivative w/ respect to XBeta
  ##########################################################
  
  if (any(negative)) {
    
    prob_negative <- prob_one[negative]
    
    derivative_xbeta[negative] <- (1 - alpha_focal)*prob_negative^gamma_focal* 
      (prob_negative - gamma_focal*(1 - prob_negative)*log(1 - prob_negative))
  }
  
  as.vector(crossprod(X_data, derivative_xbeta))}

############################################################
# Prior and Gibbs target
############################################################

log_prior <- function(beta) {-sum(beta^2/prior_sd^2)/2}

log_target <- function(beta) {log_prior(beta) - eta * 
    focal_loss_total(beta = beta, X_data = X, y_data = y)}

neg_log_target <- function(beta) {-log_target(beta)}

neg_log_target_grad <- function(beta) {beta/prior_sd^2 + eta * 
    focal_grad(beta = beta, X_data = X, y_data = y)}

############################################################
# Estimate the focal risk minimiser
############################################################

logit_estimate <- function(X_data,y_data) {
  
  fit <- suppressWarnings(glm.fit(x = X_data, y = y_data, family = binomial()))
  
  beta_initial <- fit$coefficients
  
  beta_initial[!is.finite(beta_initial)] <- 0
  
  as.vector(beta_initial)
}

fit_focal_classifier <- function(X_data,y_data,start) {
  
  fit <- optim(par = start, fn = focal_loss_total, gr = focal_grad,
               X_data = X_data, y_data = y_data, method = "BFGS",
               control = list(maxit = 1000, reltol = 1e-10)
  )
  
  if (fit$convergence != 0 || any(!is.finite(fit$par))) {
    stop("Focal loss optimisation failed.")
  }
  
  as.vector(fit$par)
}

beta_initial <- logit_estimate(X_data = X, y_data = y)

beta_hat <- fit_focal_classifier(X_data = X, y_data = y, start = beta_initial)

names(beta_hat) <- colnames(X)

############################################################
# 7. Normalized focal pseudo law
############################################################

focal_pseudo_prob <- function(beta,X_data) {
  
  #score <- clip_prob(plogis(as.vector(X_data%*%beta)))
  score <- exp(as.vector(X%*%beta_true))/(1+exp(as.vector(X%*%beta_true)))
  
  # Loss if the candidate label = 1
  loss_one <- -alpha_focal*(1 - score)^gamma_focal*log(score)
  
  # Loss if the candidate label = 0
  loss_zero <- -(1 - alpha_focal)*score^gamma_focal*log(1 - score)
  
  # log u_1 - log u_0
  log_kernel_ratio <- eta*(loss_zero - loss_one)
  
  #plogis(log_kernel_ratio)
  exp(as.vector(log_kernel_ratio))/(1+exp(as.vector(log_kernel_ratio)))
}

p_hat_positive <- focal_pseudo_prob(beta = beta_hat,X_data = X)

############################################################
# Draw estimators from the fitted bootstrap law
############################################################

draw_bootstrap_estimators <- function(B, beta_center) {
  
  pseudo_prob <- focal_pseudo_prob(beta = beta_center, X_data = X)
  
  beta_boot <- matrix(NA_real_, nrow = B,ncol = p)
  
  for (b in seq_len(B)) {
    
    y_star <- rbinom(n = n,size = 1,prob = pseudo_prob)
    
    beta_boot[b, ] <- fit_focal_classifier(X_data = X, y_data = y_star, 
                                           start = beta_center)
  }
  
  colnames(beta_boot) <- colnames(X)
  
  beta_boot
}

B <- 1500

beta_boot <- draw_bootstrap_estimators(B = B, beta_center = beta_hat)

log_target_values <- apply(beta_boot, 1, log_target)




############################################################
# Multivariate t approximation to the bootstrap estimator law
############################################################
#
# The conditional label law is available exactly, but the
# induced distribution of the focal risk minimiser is not.
#
# We approximate:
#
# q_boot(beta) = Law(beta_hat* | beta_hat)
#
############################################################

log_mean_exp <- function(values) {
  
  maximum <- max(values)
  
  maximum +
    log(
      mean(
        exp(values - maximum)
      )
    )
}

############################################################
# 10. Bootstrap importance sampler
############################################################

bootstrap_is <- function(
    B = 1500
) {
  
  start_time <- proc.time()[3]

  ##########################################################
  # Independent proposal draws
  ##########################################################
  
  beta_boot <- draw_bootstrap_estimators(
    B = B,
    beta_center = beta_hat
  )
  
  ##########################################################
  # Evaluate Gibbs target
  ##########################################################
  
  log_target_values <- apply(
    beta_boot,
    1,
    log_target
  )
  
  ##########################################################
  # Evaluate fitted bootstrap estimator law
  ##########################################################
  
  bootstrap_location <- colMeans(beta_boot)
  bootstrap_covariance <- cov(beta_boot)
  
  degrees_freedom <- 5
  
  bootstrap_scale <-
    bootstrap_covariance *
    (degrees_freedom - 2) /
    degrees_freedom
  
  log_q_boot_values <- LaplacesDemon::dmvt(
    x = beta_boot,
    mu = bootstrap_location,
    S = bootstrap_scale,
    df = degrees_freedom,
    log = TRUE
  )
  
  ##########################################################
  # Importance weights
  ##########################################################
  
  log_weights <-
    log_target_values -
    log_q_boot_values
  
  log_weights <-
    log_weights -
    max(log_weights)
  
  raw_weights <- exp(
    log_weights
  )
  
  weights <- raw_weights /
    sum(raw_weights)
  
  ess <- 1 / sum(
    weights^2
  )
  
  ##########################################################
  # Weighted posterior summaries
  ##########################################################
  
  posterior_mean <- colSums(
    beta_boot * weights
  )
  
  centered_draws <- sweep(
    beta_boot,
    2,
    posterior_mean,
    "-"
  )
  
  posterior_covariance <-
    crossprod(
      centered_draws,
      centered_draws * weights
    )
  
  posterior_summary <- data.frame(
    parameter = colnames(X),
    mean = posterior_mean,
    sd = sqrt(
      diag(posterior_covariance)
    ),
    row.names = NULL
  )
  
  elapsed_time <-
    proc.time()[3] -
    start_time
  
  list(
    method = "Focal bootstrap IS",
    beta_hat = beta_hat,
    draws = beta_boot,
    weights = weights,
    summary = posterior_summary,
    ess = ess,
    maximum_weight = max(weights),
    elapsed = elapsed_time,
    ess_per_second = ess / elapsed_time
  )
}

############################################################
# Random Walk Metropolis Hastings
############################################################

rwmh <- function(n_iter = 30000, burnin = 5000) {
  
  start_time <- proc.time()[3]
  
  ##########################################################
  # Posterior mode and local covariance
  ##########################################################
  
  posterior_mode_fit <- optim(
    par = beta_hat,
    fn = neg_log_target,
    gr = neg_log_target_grad,
    method = "BFGS",
    control = list(maxit = 1000, reltol = 1e-10)
  )
  
  posterior_mode <- as.vector(posterior_mode_fit$par)
  
  posterior_hessian <- optimHess(
    par = posterior_mode,
    fn = neg_log_target,
    gr = neg_log_target_grad)
  
  local_covariance <- solve(posterior_hessian)
  
  proposal_scale <- 2.38^2/p
  
  proposal_covariance <-proposal_scale * local_covariance
  
  ##########################################################
  # Metropolis Hastings iterations
  ##########################################################
  
  draws <- matrix(NA_real_, nrow = n_iter, ncol = p)
  
  current_beta <- posterior_mode
  
  current_log_target <- log_target(current_beta)
  
  accepted <- 0
  
  for (iteration in seq_len(n_iter)) {
    
    proposed_beta <- as.vector(
      mnormt::rmnorm(n = 1, mean = current_beta, varcov = proposal_covariance))
    
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
    parameter = colnames(X),
    mean = colMeans(draws_keep),
    sd = apply(draws_keep, 2, sd),
    ess = ess_by_parameter,
    row.names = NULL
  )
  
  overall_ess <- min(ess_by_parameter)
  
  elapsed_time <-proc.time()[3] - start_time
  
  list(
    method = "Random walk MH",
    draws = draws_keep,
    summary = posterior_summary,
    posterior_mode = posterior_mode,
    local_covariance = local_covariance,
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

bootstrap_result <- bootstrap_is(B = 1500)

set.seed(789)

mcmc_result <- rwmh(n_iter = 30000, burnin = 5000)

############################################################
# Computational comparison
############################################################

cost_table <- data.frame(
  method = c(bootstrap_result$method, mcmc_result$method),
  elapsed_seconds = c(bootstrap_result$elapsed, mcmc_result$elapsed),
  ESS = c(bootstrap_result$ess, mcmc_result$ess),
  ESS_per_second = c(bootstrap_result$ess_per_second, mcmc_result$ess_per_second)
)

comparison_metrics <- data.frame(
  metric = c(
    "Bootstrap minus MCMC time",
    "Bootstrap minus MCMC ESS per second",
    "Bootstrap divided by MCMC ESS per second"
  ),
  value = c(
    bootstrap_result$elapsed - mcmc_result$elapsed,
    bootstrap_result$ess_per_second - mcmc_result$ess_per_second,
    bootstrap_result$ess_per_second/mcmc_result$ess_per_second)
)

posterior_comparison <- data.frame(
  parameter = colnames(X),
  bootstrap_mean = bootstrap_result$summary$mean, 
  mcmc_mean = mcmc_result$summary$mean,
  absolute_mean_difference = abs(bootstrap_result$summary$mean - mcmc_result$summary$mean),
  bootstrap_sd = bootstrap_result$summary$sd,
  mcmc_sd = mcmc_result$summary$sd
)

############################################################
# 16. Overlay posterior densities
############################################################

old_par <- par(no.readonly = TRUE)

par(mfrow = c(1, p), mar = c(5, 4, 4, 1))

for (j in seq_len(p)) {
  
  bootstrap_draws_j <- bootstrap_result$draws[, j]
  
  mcmc_draws_j <- mcmc_result$draws[, j]
  
  plot_range <- range(bootstrap_draws_j, mcmc_draws_j, finite = TRUE)
  
  bootstrap_density <- topolow::weighted_kde(
    x = bootstrap_draws_j, 
    weights = bootstrap_result$weights,
    n = 500,
    from = plot_range[1],
    to = plot_range[2]
    )
  
  mcmc_density <- topolow::weighted_kde(
      x = mcmc_draws_j,
      weights = rep(1 / length(mcmc_draws_j), length(mcmc_draws_j)),
      n = 500,
      from = plot_range[1],
      to = plot_range[2]
    )
  
  y_limit <- c(0, 1.05 * max(bootstrap_density$y, mcmc_density$y, na.rm = TRUE))
  
  plot(
    x = bootstrap_density$x,
    y = bootstrap_density$y,
    type = "l",
    lwd = 2,
    col = "blue",
    xlab = colnames(X)[j],
    ylab = "Posterior density",
    main = colnames(X)[j],
    xlim = plot_range,
    ylim = y_limit
  )
  
  lines(x = mcmc_density$x, y = mcmc_density$y, lwd = 2, col = "red", lty = 2)
  
  legend("topright", 
         legend = c("Bootstrap IS", "Random walk MH", "True value"),
         col = c("blue","red"), lty = c(1, 2, 3), lwd = 2, bty = "n", cex = 0.75)
}

par(old_par)

############################################################
# 17. Print results
############################################################
cat("\nFocal risk minimiser:\n")
print(beta_hat)

cat("\nBootstrap importance sampling summary:\n")
print(bootstrap_result$summary)

cat("\nBootstrap maximum normalized weight:\n")
print(bootstrap_result$maximum_weight)

cat("\nMCMC summary:\n")
print(mcmc_result$summary)

cat("\nMCMC acceptance rate:\n")
print(mcmc_result$acceptance_rate)

cat("\nComputational comparison:\n")
print(cost_table)

cat("\nEfficiency differences:\n")
print(comparison_metrics)

cat("\nPosterior comparison:\n")
print(posterior_comparison)