############################################################
# Generalized Bayesian stochastic volatility example
#
# Correctly specified Kim, Shephard and Chib model
#
# Standard Bayes:
#   eta_R = 1
#   eta_V = 1
#   posterior computed by MCMC
#
# Generalized Bayes:
#   eta_R = 1
#   eta_V allowed to vary
#   posterior computed using:
#
#   1. pseudo likelihood bootstrap
#   2. inner Laplace approximation over h
#   3. outer Laplace proposal over theta
#   4. importance sampling correction
#
# Structural parameters:
#
# theta = (mu, phi, sigma)
#
# where
#
# y_t | h_t ~ N(0, exp(h_t))
#
# h_1 ~ N(mu, sigma^2 / (1 - phi^2))
#
# h_t | h_{t-1}
#   ~ N(
#       mu + phi(h_{t-1} - mu),
#       sigma^2
#     )
#
############################################################


############################################################
# Packages
############################################################

library(mnormt)
library(Matrix)
library(LaplacesDemon)


############################################################
# Simulation settings
############################################################

set.seed(123)

T_obs <- 200


############################################################
# True stochastic volatility parameters
############################################################

mu_true <- -1.02
phi_true <- 0.95
sigma_true <- 0.25

theta_true <-
  c(
    mu = mu_true,
    phi = phi_true,
    sigma = sigma_true
  )


############################################################
# Simulate correctly specified stochastic volatility model
############################################################

simulate_sv <- function(
    T_obs,
    theta
) {
  
  mu <- theta["mu"]
  phi <- theta["phi"]
  sigma <- theta["sigma"]
  
  
  if (
    !is.finite(mu) ||
    !is.finite(phi) ||
    !is.finite(sigma) ||
    abs(phi) >= 1 ||
    sigma <= 0
  ) {
    
    stop(
      "Invalid stochastic volatility parameters."
    )
  }
  
  
  h <- numeric(T_obs)
  y <- numeric(T_obs)
  
  
  ##########################################################
  # Stationary initial volatility
  ##########################################################
  
  h[1] <-
    rnorm(
      1,
      mean = mu,
      sd =
        sigma /
        sqrt(
          1 - phi^2
        )
    )
  
  
  ##########################################################
  # Remaining latent volatility path
  ##########################################################
  
  if (T_obs > 1) {
    
    for (t in 2:T_obs) {
      
      h[t] <-
        rnorm(
          1,
          mean =
            mu +
            phi *
            (
              h[t - 1] -
                mu
            ),
          sd = sigma
        )
    }
  }
  
  
  ##########################################################
  # Mean corrected returns
  ##########################################################
  
  for (t in seq_len(T_obs)) {
    
    y[t] <-
      rnorm(
        1,
        mean = 0,
        sd =
          exp(
            h[t] / 2
          )
      )
  }
  
  
  list(
    y = y,
    h = h
  )
}


simulated_data <-
  simulate_sv(
    T_obs = T_obs,
    theta = theta_true
  )


y <- simulated_data$y

h_true <- simulated_data$h


############################################################
# Priors
#
# Same basic prior structure as the Stan example:
#
# mu    ~ Cauchy(0, 10)
# phi   ~ Uniform(-1, 1)
# sigma ~ half-Cauchy(0, 5)
############################################################

log_prior_theta <- function(theta) {
  
  mu <- theta[1]
  phi <- theta[2]
  sigma <- theta[3]
  
  
  if (
    !is.finite(mu) ||
    !is.finite(phi) ||
    !is.finite(sigma) ||
    abs(phi) >= 1 ||
    sigma <= 0
  ) {
    
    return(-Inf)
  }
  
  
  log_prior_mu <-
    dcauchy(
      mu,
      location = 0,
      scale = 10,
      log = TRUE
    )
  
  
  log_prior_phi <-
    dunif(
      phi,
      min = -1,
      max = 1,
      log = TRUE
    )
  
  
  # Half Cauchy:
  # factor 2 is required because sigma > 0
  
  log_prior_sigma <-
    log(2) +
    dcauchy(
      sigma,
      location = 0,
      scale = 5,
      log = TRUE
    )
  
  
  log_prior_mu +
    log_prior_phi +
    log_prior_sigma
}


############################################################
# Return loss
#
# y_t | h_t ~ N(0, exp(h_t))
#
# Negative log likelihood up to constants:
#
# L_R(h; y)
#
# = 1/2 sum_t [
#
#       h_t
#
#       +
#
#       y_t^2 exp(-h_t)
#
#   ]
############################################################

loss_return <- function(
    h,
    y_data
) {
  
  0.5 *
    sum(
      h +
        y_data^2 *
        exp(-h)
    )
}


############################################################
# Volatility loss
#
# h_1 ~ N(
#          mu,
#          sigma^2 / (1 - phi^2)
#        )
#
# h_t | h_{t-1}
#     ~ N(
#          mu + phi(h_{t-1} - mu),
#          sigma^2
#        )
#
# Constants independent of theta and h are omitted.
############################################################

loss_volatility <- function(
    theta,
    h
) {
  
  mu <- theta[1]
  phi <- theta[2]
  sigma <- theta[3]
  
  
  if (
    !is.finite(mu) ||
    !is.finite(phi) ||
    !is.finite(sigma) ||
    abs(phi) >= 1 ||
    sigma <= 0
  ) {
    
    return(1e100)
  }
  
  
  T_local <- length(h)
  
  
  ##########################################################
  # Initial state contribution
  ##########################################################
  
  initial_variance <-
    sigma^2 /
    (
      1 - phi^2
    )
  
  
  initial_loss <-
    0.5 *
    (
      log(initial_variance) +
        (
          h[1] -
            mu
        )^2 /
        initial_variance
    )
  
  
  ##########################################################
  # Transition contributions
  ##########################################################
  
  if (T_local == 1) {
    
    return(initial_loss)
  }
  
  
  transition_mean <-
    mu +
    phi *
    (
      h[1:(T_local - 1)] -
        mu
    )
  
  
  innovation <-
    h[2:T_local] -
    transition_mean
  
  
  transition_loss <-
    0.5 *
    sum(
      log(sigma^2) +
        innovation^2 /
        sigma^2
    )
  
  
  initial_loss +
    transition_loss
}


############################################################
# Total generalized loss
#
# eta_R = 1
#
# eta_V can vary
############################################################

loss_generalized <- function(
    theta,
    h,
    y_data,
    eta_r = 1,
    eta_v = 1
) {
  
  eta_r *
    loss_return(
      h = h,
      y_data = y_data
    ) +
    
    eta_v *
    loss_volatility(
      theta = theta,
      h = h
    )
}


############################################################
# Gradient of generalized loss with respect to h
#
# Needed for inner Laplace approximation
############################################################

gradient_h <- function(
    h,
    theta,
    y_data,
    eta_r = 1,
    eta_v = 1
) {
  
  mu <- theta[1]
  phi <- theta[2]
  sigma <- theta[3]
  
  T_local <- length(h)
  
  
  ##########################################################
  # Return contribution
  ##########################################################
  
  gradient <-
    0.5 *
    eta_r *
    (
      1 -
        y_data^2 *
        exp(-h)
    )
  
  
  inverse_sigma2 <-
    1 /
    sigma^2
  
  
  ##########################################################
  # Initial volatility contribution
  ##########################################################
  
  gradient[1] <-
    gradient[1] +
    eta_v *
    (
      1 - phi^2
    ) *
    inverse_sigma2 *
    (
      h[1] -
        mu
    )
  
  
  ##########################################################
  # Transition contributions
  ##########################################################
  
  if (T_local > 1) {
    
    innovation <-
      h[2:T_local] -
      mu -
      phi *
      (
        h[1:(T_local - 1)] -
          mu
      )
    
    
    gradient[1] <-
      gradient[1] -
      eta_v *
      phi *
      inverse_sigma2 *
      innovation[1]
    
    
    if (T_local > 2) {
      
      for (t in 2:(T_local - 1)) {
        
        gradient[t] <-
          gradient[t] +
          eta_v *
          inverse_sigma2 *
          (
            innovation[t - 1] -
              phi *
              innovation[t]
          )
      }
    }
    
    
    gradient[T_local] <-
      gradient[T_local] +
      eta_v *
      inverse_sigma2 *
      innovation[T_local - 1]
  }
  
  
  gradient
}


############################################################
# Hessian with respect to latent volatility h
#
# Tridiagonal structure is exploited.
############################################################

hessian_h <- function(
    h,
    theta,
    y_data,
    eta_r = 1,
    eta_v = 1
) {
  
  phi <- theta[2]
  sigma <- theta[3]
  
  T_local <- length(h)
  
  
  ##########################################################
  # Return Hessian contribution
  ##########################################################
  
  diagonal <-
    0.5 *
    eta_r *
    y_data^2 *
    exp(-h)
  
  
  inverse_sigma2 <-
    1 /
    sigma^2
  
  
  ##########################################################
  # Volatility Hessian
  ##########################################################
  
  if (T_local == 1) {
    
    diagonal[1] <-
      diagonal[1] +
      eta_v *
      (
        1 - phi^2
      ) *
      inverse_sigma2
    
    
    return(
      Matrix::Diagonal(
        x = diagonal
      )
    )
  }
  
  
  # First state:
  #
  # initial contribution = 1 - phi^2
  # outgoing transition contribution = phi^2
  #
  # total = 1
  
  diagonal[1] <-
    diagonal[1] +
    eta_v *
    inverse_sigma2
  
  
  # Interior states
  
  if (T_local > 2) {
    
    diagonal[2:(T_local - 1)] <-
      diagonal[2:(T_local - 1)] +
      eta_v *
      inverse_sigma2 *
      (
        1 +
          phi^2
      )
  }
  
  
  # Final state
  
  diagonal[T_local] <-
    diagonal[T_local] +
    eta_v *
    inverse_sigma2
  
  
  off_diagonal <-
    rep(
      -eta_v *
        phi *
        inverse_sigma2,
      T_local - 1
    )
  
  
  Matrix::bandSparse(
    n = T_local,
    k = c(-1, 0, 1),
    diagonals =
      list(
        off_diagonal,
        diagonal,
        off_diagonal
      )
  )
}


############################################################
# Starting value for latent log volatility
############################################################

initial_h_guess <- function(
    theta,
    y_data
) {
  
  mu <- theta[1]
  
  
  small_value <-
    max(
      0.05 *
        var(y_data),
      1e-6
    )
  
  
  empirical_h <-
    log(
      y_data^2 +
        small_value
    )
  
  
  # Shrink noisy log squared returns towards mu
  
  h_start <-
    0.50 *
    empirical_h +
    0.50 *
    mu
  
  
  pmin(
    pmax(
      h_start,
      -12
    ),
    8
  )
}


############################################################
# Conditional mode of h for fixed theta
#
# h_hat(theta)
#
# = argmin_h [
#
#     eta_R L_R
#
#     +
#
#     eta_V L_V
#
#   ]
############################################################

find_h_mode <- function(
    theta,
    y_data,
    eta_r = 1,
    eta_v = 1,
    start_h = NULL
) {
  
  if (
    abs(theta[2]) >= 1 ||
    theta[3] <= 0 ||
    eta_r <= 0 ||
    eta_v <= 0
  ) {
    
    return(
      list(
        convergence = 1,
        h_mode = rep(
          NA_real_,
          length(y_data)
        ),
        value = Inf
      )
    )
  }
  
  
  if (is.null(start_h)) {
    
    start_h <-
      initial_h_guess(
        theta = theta,
        y_data = y_data
      )
  }
  
  
  fit <-
    optim(
      par = start_h,
      
      fn = function(h) {
        
        loss_generalized(
          theta = theta,
          h = h,
          y_data = y_data,
          eta_r = eta_r,
          eta_v = eta_v
        )
      },
      
      gr = function(h) {
        
        gradient_h(
          h = h,
          theta = theta,
          y_data = y_data,
          eta_r = eta_r,
          eta_v = eta_v
        )
      },
      
      method = "BFGS",
      
      control =
        list(
          maxit = 1000,
          reltol = 1e-10
        )
    )
  
  
  list(
    convergence = fit$convergence,
    h_mode = fit$par,
    value = fit$value
  )
}


############################################################
# Inner Laplace approximation
#
# M_eta(theta; y)
#
# = integral
#
# exp{
#   - eta_R L_R(h; y)
#   - eta_V L_V(theta; h)
# }
#
# dh
############################################################

log_marginal_pseudo <- function(
    theta,
    y_data,
    eta_r = 1,
    eta_v = 1
) {
  
  mu <- theta[1]
  phi <- theta[2]
  sigma <- theta[3]
  
  
  if (
    !is.finite(mu) ||
    !is.finite(phi) ||
    !is.finite(sigma) ||
    abs(phi) >= 0.999 ||
    sigma <= 0 ||
    eta_r <= 0 ||
    eta_v <= 0
  ) {
    
    return(-Inf)
  }
  
  
  h_fit <-
    find_h_mode(
      theta = theta,
      y_data = y_data,
      eta_r = eta_r,
      eta_v = eta_v
    )
  
  
  if (
    h_fit$convergence != 0 ||
    any(
      !is.finite(
        h_fit$h_mode
      )
    )
  ) {
    
    return(-Inf)
  }
  
  
  H_h <-
    hessian_h(
      h = h_fit$h_mode,
      theta = theta,
      y_data = y_data,
      eta_r = eta_r,
      eta_v = eta_v
    )
  
  
  determinant_H <-
    determinant(
      H_h,
      logarithm = TRUE
    )
  
  
  if (
    determinant_H$sign <= 0
  ) {
    
    return(-Inf)
  }
  
  
  log_det_H <-
    as.numeric(
      determinant_H$modulus
    )
  
  
  T_local <-
    length(y_data)
  
  
  ##########################################################
  # Laplace approximation
  ##########################################################
  
  -h_fit$value -
    0.5 *
    log_det_H +
    0.5 *
    T_local *
    log(
      2 *
        pi
    )
}


############################################################
# Fit marginal pseudo likelihood
#
# theta_hat_eta =
#
# argmax_theta M_eta(theta; y)
############################################################

fit_pseudo_theta <- function(
    y_data,
    eta_r = 1,
    eta_v = 1,
    start = NULL
) {
  
  if (is.null(start)) {
    
    start <-
      c(
        mu =
          log(
            var(y_data) +
              1e-6
          ),
        
        phi =
          0.90,
        
        sigma =
          0.30
      )
  }
  
  
  objective <- function(theta) {
    
    value <-
      log_marginal_pseudo(
        theta = theta,
        y_data = y_data,
        eta_r = eta_r,
        eta_v = eta_v
      )
    
    
    if (!is.finite(value)) {
      
      return(1e100)
    }
    
    
    -value
  }
  
  
  fit <-
    optim(
      par = start,
      
      fn = objective,
      
      method = "L-BFGS-B",
      
      lower =
        c(
          -8,
          -0.995,
          0.02
        ),
      
      upper =
        c(
          4,
          0.995,
          2
        ),
      
      control =
        list(
          maxit = 300,
          parscale =
            c(
              1,
              0.05,
              0.10
            )
        )
    )
  
  
  names(fit$par) <-
    c(
      "mu",
      "phi",
      "sigma"
    )
  
  
  list(
    par = fit$par,
    log_pseudo =
      -fit$value,
    convergence =
      fit$convergence,
    message =
      fit$message
  )
}


############################################################
# Positive definite covariance regularization
############################################################

make_positive_definite <- function(
    Sigma,
    min_eigenvalue = 1e-8
) {
  
  Sigma <-
    (
      Sigma +
        t(Sigma)
    ) /
    2
  
  
  eig <-
    eigen(
      Sigma,
      symmetric = TRUE
    )
  
  
  eig$values <-
    pmax(
      eig$values,
      min_eigenvalue
    )
  
  
  Sigma_pd <-
    eig$vectors %*%
    diag(
      eig$values
    ) %*%
    t(
      eig$vectors
    )
  
  
  (
    Sigma_pd +
      t(Sigma_pd)
  ) /
    2
}


############################################################
# Outer Laplace approximation
#
# q_eta(theta)
#
# approximately
#
# N(
#   theta_hat_eta,
#   J_eta^{-1}
# )
#
# where
#
# J_eta =
#
# - d^2 / d theta^2
#
# log M_eta(theta; y)
############################################################

laplace_parameter_proposal <- function(
    y_data,
    eta_r = 1,
    eta_v = 1,
    proposal_inflation = 1.15,
    start = NULL
) {
  
  pseudo_fit <-
    fit_pseudo_theta(
      y_data = y_data,
      eta_r = eta_r,
      eta_v = eta_v,
      start = start
    )
  
  
  if (
    pseudo_fit$convergence != 0
  ) {
    
    warning(
      "Pseudo likelihood optimization did not return convergence code 0."
    )
  }
  
  
  theta_hat <-
    pseudo_fit$par
  
  
  negative_log_pseudo <- function(theta) {
    
    value <-
      log_marginal_pseudo(
        theta = theta,
        y_data = y_data,
        eta_r = eta_r,
        eta_v = eta_v
      )
    
    
    if (!is.finite(value)) {
      
      return(1e100)
    }
    
    
    -value
  }
  
  
  J_eta <-
    optimHess(
      par = theta_hat,
      fn = negative_log_pseudo
    )
  
  
  J_eta <-
    (
      J_eta +
        t(J_eta)
    ) /
    2
  
  
  ##########################################################
  # Force positive definite local precision
  ##########################################################
  
  eig <-
    eigen(
      J_eta,
      symmetric = TRUE
    )
  
  
  eig$values <-
    pmax(
      eig$values,
      1e-8
    )
  
  
  J_eta_pd <-
    eig$vectors %*%
    diag(
      eig$values
    ) %*%
    t(
      eig$vectors
    )
  
  
  laplace_covariance <-
    solve(
      J_eta_pd
    )
  
  
  proposal_covariance <-
    proposal_inflation^2 *
    laplace_covariance
  
  
  proposal_covariance <-
    make_positive_definite(
      proposal_covariance
    )
  
  
  dimnames(
    proposal_covariance
  ) <-
    list(
      c(
        "mu",
        "phi",
        "sigma"
      ),
      c(
        "mu",
        "phi",
        "sigma"
      )
    )
  
  
  list(
    theta_hat = theta_hat,
    precision = J_eta_pd,
    laplace_covariance =
      laplace_covariance,
    proposal_covariance =
      proposal_covariance,
    pseudo_fit =
      pseudo_fit
  )
}


############################################################
# Sequential pseudo likelihood bootstrap
#
# Stage 1:
#
# h* generated from normalized
#
# exp{
#   - eta_V L_V
# }
#
# Stage 2:
#
# y* generated from normalized
#
# exp{
#   - eta_R L_R
# }
############################################################

simulate_pseudo_sv <- function(
    T_obs,
    theta,
    eta_r = 1,
    eta_v = 1
) {
  
  mu <- theta[1]
  phi <- theta[2]
  sigma <- theta[3]
  
  
  if (
    abs(phi) >= 1 ||
    sigma <= 0 ||
    eta_r <= 0 ||
    eta_v <= 0
  ) {
    
    stop(
      "Invalid parameter or learning rate."
    )
  }
  
  
  h_star <- numeric(T_obs)
  y_star <- numeric(T_obs)
  
  
  ##########################################################
  # Volatility pseudo likelihood
  #
  # Raising Gaussian loss to eta_V and normalizing gives:
  #
  # variance / eta_V
  ##########################################################
  
  h_star[1] <-
    rnorm(
      1,
      
      mean = mu,
      
      sd =
        sigma /
        sqrt(
          eta_v *
            (
              1 - phi^2
            )
        )
    )
  
  
  if (T_obs > 1) {
    
    for (t in 2:T_obs) {
      
      h_star[t] <-
        rnorm(
          1,
          
          mean =
            mu +
            phi *
            (
              h_star[t - 1] -
                mu
            ),
          
          sd =
            sigma /
            sqrt(
              eta_v
            )
        )
    }
  }
  
  
  ##########################################################
  # Return pseudo likelihood
  #
  # y* | h*
  #
  # ~ N(
  #     0,
  #     exp(h*) / eta_R
  #   )
  #
  # Therefore SD:
  #
  # exp(h*/2) / sqrt(eta_R)
  ##########################################################
  
  for (t in seq_len(T_obs)) {
    
    y_star[t] <-
      rnorm(
        1,
        
        mean = 0,
        
        sd =
          exp(
            h_star[t] /
              2
          ) /
          sqrt(
            eta_r
          )
      )
  }
  
  
  list(
    y = y_star,
    h = h_star
  )
}


############################################################
# Pseudo likelihood bootstrap
#
# The latent bootstrap path is NOT treated as observed when
# theta is re-estimated.
#
# Each bootstrap return series is fitted in exactly the same
# way as the observed return series.
############################################################

pseudo_bootstrap <- function(
    y_data,
    theta_hat,
    eta_r = 1,
    eta_v = 1,
    B_boot = 50
) {
  
  bootstrap_estimators <-
    matrix(
      NA_real_,
      nrow = B_boot,
      ncol = 3
    )
  
  
  colnames(
    bootstrap_estimators
  ) <-
    c(
      "mu",
      "phi",
      "sigma"
    )
  
  
  successful <- 0
  
  attempts <- 0
  
  max_attempts <-
    2 *
    B_boot
  
  
  while (
    successful < B_boot &&
    attempts < max_attempts
  ) {
    
    attempts <-
      attempts +
      1
    
    
    ########################################################
    # Sequential pseudo bootstrap
    ########################################################
    
    pseudo_data <-
      simulate_pseudo_sv(
        T_obs =
          length(
            y_data
          ),
        
        theta =
          theta_hat,
        
        eta_r =
          eta_r,
        
        eta_v =
          eta_v
      )
    
    
    ########################################################
    # Refit marginal pseudo likelihood
    ########################################################
    
    fit_star <-
      try(
        fit_pseudo_theta(
          y_data =
            pseudo_data$y,
          
          eta_r =
            eta_r,
          
          eta_v =
            eta_v,
          
          start =
            theta_hat
        ),
        
        silent = TRUE
      )
    
    
    if (
      !inherits(
        fit_star,
        "try-error"
      ) &&
      fit_star$convergence == 0 &&
      all(
        is.finite(
          fit_star$par
        )
      ) &&
      abs(
        fit_star$par["phi"]
      ) < 1 &&
      fit_star$par["sigma"] > 0
    ) {
      
      successful <-
        successful +
        1
      
      
      bootstrap_estimators[
        successful,
      ] <-
        fit_star$par
    }
  }
  
  
  if (successful == 0) {
    
    stop(
      "No successful pseudo bootstrap fits."
    )
  }
  
  
  bootstrap_estimators <-
    bootstrap_estimators[
      seq_len(
        successful
      ),
      ,
      drop = FALSE
    ]
  
  
  list(
    estimates =
      bootstrap_estimators,
    
    mean =
      colMeans(
        bootstrap_estimators
      ),
    
    covariance =
      cov(
        bootstrap_estimators
      ),
    
    successful =
      successful,
    
    attempts =
      attempts
  )
}


############################################################
# Draw from Laplace proposal subject to model support
#
# The normal proposal is implicitly truncated to
#
# -1 < phi < 1
# sigma > 0
#
# The truncation normalizing constant is constant across
# valid draws and therefore cancels from normalized
# importance weights.
############################################################

draw_valid_proposal <- function(
    n,
    mean,
    covariance
) {
  
  draws <-
    matrix(
      NA_real_,
      nrow = n,
      ncol = 3
    )
  
  
  colnames(draws) <-
    c(
      "mu",
      "phi",
      "sigma"
    )
  
  
  accepted <- 0
  
  
  while (accepted < n) {
    
    remaining <-
      n -
      accepted
    
    
    batch_size <-
      max(
        100,
        2 *
          remaining
      )
    
    
    candidates <-
      mnormt::rmnorm(
        n = batch_size,
        mean = mean,
        varcov = covariance
      )
    
    
    valid <-
      abs(
        candidates[, 2]
      ) < 0.999 &
      candidates[, 3] > 0
    
    
    candidates <-
      candidates[
        valid,
        ,
        drop = FALSE
      ]
    
    
    if (
      nrow(
        candidates
      ) == 0
    ) {
      
      next
    }
    
    
    number_to_take <-
      min(
        remaining,
        nrow(
          candidates
        )
      )
    
    
    index <-
      (
        accepted +
          1
      ):(
        accepted +
          number_to_take
      )
    
    
    draws[
      index,
    ] <-
      candidates[
        seq_len(
          number_to_take
        ),
        ,
        drop = FALSE
      ]
    
    
    accepted <-
      accepted +
      number_to_take
  }
  
  
  draws
}


############################################################
# Weighted covariance
############################################################

weighted_covariance <- function(
    draws,
    weights
) {
  
  weights <-
    weights /
    sum(
      weights
    )
  
  
  weighted_mean <-
    colSums(
      draws *
        weights
    )
  
  
  centered <-
    sweep(
      draws,
      MARGIN = 2,
      STATS = weighted_mean,
      FUN = "-"
    )
  
  
  numerator <-
    crossprod(
      centered *
        sqrt(
          weights
        )
    )
  
  
  correction <-
    1 -
    sum(
      weights^2
    )
  
  
  if (
    correction <= 0
  ) {
    
    return(
      matrix(
        NA_real_,
        ncol(draws),
        ncol(draws)
      )
    )
  }
  
  
  numerator /
    correction
}


############################################################
# Generalized Bayesian PBIS
############################################################

pbis_sv <- function(
    y_data,
    eta_v,
    eta_r = 1,
    B_boot = 50,
    B_is = 1000,
    proposal_inflation = 1.15
) {
  
  start_time <-
    proc.time()[3]
  
  
  ##########################################################
  # Step 1:
  # Laplace proposal based on observed pseudo likelihood
  ##########################################################
  
  proposal <-
    laplace_parameter_proposal(
      y_data = y_data,
      eta_r = eta_r,
      eta_v = eta_v,
      proposal_inflation =
        proposal_inflation
    )
  
  
  theta_hat <-
    proposal$theta_hat
  
  
  ##########################################################
  # Step 2:
  # Pseudo likelihood bootstrap
  ##########################################################
  
  bootstrap_result <-
    pseudo_bootstrap(
      y_data = y_data,
      theta_hat = theta_hat,
      eta_r = eta_r,
      eta_v = eta_v,
      B_boot = B_boot
    )
  
  
  ##########################################################
  # Step 3:
  # Importance proposal draws from outer Laplace
  ##########################################################
  
  proposal_draws <-
    draw_valid_proposal(
      n = B_is,
      mean = theta_hat,
      covariance =
        proposal$proposal_covariance
    )
  
  
  ##########################################################
  # Step 4:
  # Evaluate marginal generalized target
  #
  # target(theta)
  #
  # proportional to
  #
  # prior(theta) *
  #
  # M_eta(theta; y)
  ##########################################################
  
  log_marginal_values <-
    numeric(
      B_is
    )
  
  
  for (s in seq_len(B_is)) {
    
    log_marginal_values[s] <-
      log_marginal_pseudo(
        theta =
          proposal_draws[s, ],
        
        y_data =
          y_data,
        
        eta_r =
          eta_r,
        
        eta_v =
          eta_v
      )
  }
  
  
  log_prior_values <-
    apply(
      proposal_draws,
      MARGIN = 1,
      FUN =
        log_prior_theta
    )
  
  
  log_target_values <-
    log_prior_values +
    log_marginal_values
  
  
  ##########################################################
  # Step 5:
  # Laplace proposal density
  ##########################################################
  
  log_proposal_values <-
    mnormt::dmnorm(
      x =
        proposal_draws,
      
      mean =
        theta_hat,
      
      varcov =
        proposal$proposal_covariance,
      
      log =
        TRUE
    )
  
  
  ##########################################################
  # Step 6:
  # Importance weights
  ##########################################################
  
  log_weights <-
    log_target_values -
    log_proposal_values
  
  
  finite <-
    is.finite(
      log_weights
    )
  
  
  if (
    !any(
      finite
    )
  ) {
    
    stop(
      "All PBIS importance weights are non-finite."
    )
  }
  
  
  log_weights[
    !finite
  ] <-
    -Inf
  
  
  max_log_weight <-
    max(
      log_weights
    )
  
  
  raw_weights <-
    exp(
      log_weights -
        max_log_weight
    )
  
  
  weights <-
    raw_weights /
    sum(
      raw_weights
    )
  
  
  ##########################################################
  # Effective sample size
  ##########################################################
  
  ess <-
    1 /
    sum(
      weights^2
    )
  
  
  ##########################################################
  # Weighted posterior summaries
  ##########################################################
  
  post_mean <-
    colSums(
      proposal_draws *
        weights
    )
  
  
  post_covariance <-
    weighted_covariance(
      draws =
        proposal_draws,
      weights =
        weights
    )
  
  
  post_summary <-
    data.frame(
      parameter =
        c(
          "mu",
          "phi",
          "sigma"
        ),
      
      mean =
        post_mean,
      
      sd =
        sqrt(
          diag(
            post_covariance
          )
        ),
      
      row.names =
        NULL
    )
  
  
  elapsed_time <-
    proc.time()[3] -
    start_time
  
  
  list(
    method =
      paste0(
        "Generalized Bayes PBIS: eta_V = ",
        eta_v
      ),
    
    eta_r =
      eta_r,
    
    eta_v =
      eta_v,
    
    theta_hat =
      theta_hat,
    
    bootstrap_estimators =
      bootstrap_result$estimates,
    
    bootstrap_mean =
      bootstrap_result$mean,
    
    bootstrap_covariance =
      bootstrap_result$covariance,
    
    bootstrap_successful =
      bootstrap_result$successful,
    
    bootstrap_attempts =
      bootstrap_result$attempts,
    
    laplace_covariance =
      proposal$laplace_covariance,
    
    proposal_covariance =
      proposal$proposal_covariance,
    
    draws =
      proposal_draws,
    
    weights =
      weights,
    
    summary =
      post_summary,
    
    covariance =
      post_covariance,
    
    ess =
      ess,
    
    relative_ess =
      ess /
      B_is,
    
    maximum_weight =
      max(
        weights
      ),
    
    elapsed =
      elapsed_time,
    
    ess_per_second =
      ess /
      elapsed_time
  )
}


############################################################
# Standard Bayesian posterior
#
# eta_R = eta_V = 1
#
# Exact joint posterior over theta and h
############################################################

log_joint_standard <- function(
    theta,
    h,
    y_data
) {
  
  mu <- theta[1]
  phi <- theta[2]
  sigma <- theta[3]
  
  
  log_prior <-
    log_prior_theta(
      theta
    )
  
  
  if (
    !is.finite(
      log_prior
    )
  ) {
    
    return(-Inf)
  }
  
  
  ##########################################################
  # Return log likelihood
  ##########################################################
  
  log_return <-
    -0.5 *
    sum(
      h +
        y_data^2 *
        exp(-h)
    )
  
  
  ##########################################################
  # Volatility log likelihood
  ##########################################################
  
  T_local <-
    length(h)
  
  
  initial_quadratic <-
    (
      1 - phi^2
    ) *
    (
      h[1] -
        mu
    )^2
  
  
  if (T_local > 1) {
    
    innovation <-
      h[2:T_local] -
      mu -
      phi *
      (
        h[1:(T_local - 1)] -
          mu
      )
    
    
    transition_quadratic <-
      sum(
        innovation^2
      )
    
  } else {
    
    transition_quadratic <-
      0
  }
  
  
  log_volatility <-
    0.5 *
    log(
      1 - phi^2
    ) -
    
    T_local *
    log(
      sigma
    ) -
    
    0.5 /
    sigma^2 *
    (
      initial_quadratic +
        transition_quadratic
    )
  
  
  log_prior +
    log_return +
    log_volatility
}


############################################################
# Local log conditional contribution for one h_t
#
# Used for efficient single-site MCMC updates
############################################################

log_h_local <- function(
    t,
    h,
    theta,
    y_data
) {
  
  mu <- theta[1]
  phi <- theta[2]
  sigma <- theta[3]
  
  T_local <- length(h)
  
  
  value <-
    -0.5 *
    (
      h[t] +
        y_data[t]^2 *
        exp(
          -h[t]
        )
    )
  
  
  ##########################################################
  # Initial state
  ##########################################################
  
  if (t == 1) {
    
    value <-
      value -
      0.5 *
      (
        1 - phi^2
      ) *
      (
        h[1] -
          mu
      )^2 /
      sigma^2
    
    
    if (T_local > 1) {
      
      innovation_next <-
        h[2] -
        mu -
        phi *
        (
          h[1] -
            mu
        )
      
      
      value <-
        value -
        0.5 *
        innovation_next^2 /
        sigma^2
    }
    
    
    return(value)
  }
  
  
  ##########################################################
  # Incoming transition
  ##########################################################
  
  innovation_current <-
    h[t] -
    mu -
    phi *
    (
      h[t - 1] -
        mu
    )
  
  
  value <-
    value -
    0.5 *
    innovation_current^2 /
    sigma^2
  
  
  ##########################################################
  # Outgoing transition
  ##########################################################
  
  if (t < T_local) {
    
    innovation_next <-
      h[t + 1] -
      mu -
      phi *
      (
        h[t] -
          mu
      )
    
    
    value <-
      value -
      0.5 *
      innovation_next^2 /
      sigma^2
  }
  
  
  value
}


############################################################
# Standard Bayesian MCMC
#
# Metropolis within Gibbs:
#
# 1. componentwise updates for mu, phi, sigma
# 2. single-site latent volatility updates
#
# Proposal scales adapt during burnin.
############################################################

mcmc_standard_bayes <- function(
    y_data,
    n_iter = 20000,
    burnin = 5000,
    adapt_interval = 100
) {
  
  start_time <-
    proc.time()[3]
  
  
  T_local <-
    length(
      y_data
    )
  
  
  ##########################################################
  # Initial structural parameters
  ##########################################################
  
  theta <-
    c(
      mu =
        log(
          var(
            y_data
          ) +
            1e-6
        ),
      
      phi =
        0.90,
      
      sigma =
        0.30
    )
  
  
  ##########################################################
  # Initial latent volatility path
  ##########################################################
  
  initial_h_fit <-
    find_h_mode(
      theta = theta,
      y_data = y_data,
      eta_r = 1,
      eta_v = 1
    )
  
  
  h <-
    initial_h_fit$h_mode
  
  
  ##########################################################
  # Initial proposal scales
  ##########################################################
  
  theta_rw_sd <-
    c(
      mu = 0.08,
      phi = 0.01,
      sigma = 0.02
    )
  
  
  h_rw_sd <- 0.30
  
  
  ##########################################################
  # Storage
  ##########################################################
  
  n_keep <-
    n_iter -
    burnin
  
  
  theta_draws <-
    matrix(
      NA_real_,
      nrow = n_keep,
      ncol = 3
    )
  
  
  colnames(
    theta_draws
  ) <-
    c(
      "mu",
      "phi",
      "sigma"
    )
  
  
  h_posterior_sum <-
    numeric(
      T_local
    )
  
  
  ##########################################################
  # Acceptance counters
  ##########################################################
  
  total_theta_accept <-
    numeric(3)
  
  
  total_h_accept <- 0
  
  
  interval_theta_accept <-
    numeric(3)
  
  
  interval_h_accept <- 0
  
  
  current_log_target <-
    log_joint_standard(
      theta = theta,
      h = h,
      y_data = y_data
    )
  
  
  keep_index <- 0
  
  
  ##########################################################
  # MCMC loop
  ##########################################################
  
  for (iter in seq_len(n_iter)) {
    
    
    ########################################################
    # Structural parameter updates
    ########################################################
    
    for (j in seq_len(3)) {
      
      proposed_theta <-
        theta
      
      
      proposed_theta[j] <-
        rnorm(
          1,
          mean =
            theta[j],
          sd =
            theta_rw_sd[j]
        )
      
      
      proposed_log_target <-
        log_joint_standard(
          theta =
            proposed_theta,
          h =
            h,
          y_data =
            y_data
        )
      
      
      if (
        is.finite(
          proposed_log_target
        )
      ) {
        
        log_acceptance_ratio <-
          proposed_log_target -
          current_log_target
        
        
        if (
          log(
            runif(1)
          ) <
          min(
            0,
            log_acceptance_ratio
          )
        ) {
          
          theta <-
            proposed_theta
          
          
          current_log_target <-
            proposed_log_target
          
          
          total_theta_accept[j] <-
            total_theta_accept[j] +
            1
          
          
          interval_theta_accept[j] <-
            interval_theta_accept[j] +
            1
        }
      }
    }
    
    
    ########################################################
    # Latent volatility updates
    ########################################################
    
    for (t in seq_len(T_local)) {
      
      current_local <-
        log_h_local(
          t = t,
          h = h,
          theta = theta,
          y_data = y_data
        )
      
      
      proposed_h_t <-
        rnorm(
          1,
          mean = h[t],
          sd = h_rw_sd
        )
      
      
      old_h_t <- h[t]
      
      h[t] <-
        proposed_h_t
      
      
      proposed_local <-
        log_h_local(
          t = t,
          h = h,
          theta = theta,
          y_data = y_data
        )
      
      
      log_acceptance_ratio <-
        proposed_local -
        current_local
      
      
      if (
        log(
          runif(1)
        ) <
        min(
          0,
          log_acceptance_ratio
        )
      ) {
        
        total_h_accept <-
          total_h_accept +
          1
        
        
        interval_h_accept <-
          interval_h_accept +
          1
        
      } else {
        
        h[t] <-
          old_h_t
      }
    }
    
    
    ########################################################
    # Recompute exact joint target after h sweep
    ########################################################
    
    current_log_target <-
      log_joint_standard(
        theta = theta,
        h = h,
        y_data = y_data
      )
    
    
    ########################################################
    # Adapt during burnin only
    ########################################################
    
    if (
      iter <= burnin &&
      iter %% adapt_interval == 0
    ) {
      
      theta_acceptance_interval <-
        interval_theta_accept /
        adapt_interval
      
      
      for (j in seq_len(3)) {
        
        theta_rw_sd[j] <-
          theta_rw_sd[j] *
          exp(
            theta_acceptance_interval[j] -
              0.44
          )
      }
      
      
      h_acceptance_interval <-
        interval_h_accept /
        (
          adapt_interval *
            T_local
        )
      
      
      h_rw_sd <-
        h_rw_sd *
        exp(
          h_acceptance_interval -
            0.44
        )
      
      
      ######################################################
      # Prevent unreasonable proposal scales
      ######################################################
      
      theta_rw_sd <-
        pmax(
          theta_rw_sd,
          c(
            1e-4,
            1e-4,
            1e-4
          )
        )
      
      
      theta_rw_sd <-
        pmin(
          theta_rw_sd,
          c(
            2,
            0.20,
            1
          )
        )
      
      
      h_rw_sd <-
        min(
          max(
            h_rw_sd,
            0.02
          ),
          2
        )
      
      
      interval_theta_accept[] <- 0
      
      interval_h_accept <- 0
    }
    
    
    ########################################################
    # Store posterior draws
    ########################################################
    
    if (iter > burnin) {
      
      keep_index <-
        keep_index +
        1
      
      
      theta_draws[
        keep_index,
      ] <-
        theta
      
      
      h_posterior_sum <-
        h_posterior_sum +
        h
    }
  }
  
  
  ##########################################################
  # Posterior latent volatility mean
  ##########################################################
  
  h_posterior_mean <-
    h_posterior_sum /
    n_keep
  
  
  ##########################################################
  # ESS
  ##########################################################
  
  ess_parameter <-
    apply(
      theta_draws,
      2,
      LaplacesDemon::ESS
    )
  
  
  overall_ess <-
    min(
      ess_parameter
    )
  
  
  ##########################################################
  # Posterior summary
  ##########################################################
  
  posterior_summary <-
    data.frame(
      parameter =
        colnames(
          theta_draws
        ),
      
      mean =
        colMeans(
          theta_draws
        ),
      
      sd =
        apply(
          theta_draws,
          2,
          sd
        ),
      
      ess =
        ess_parameter,
      
      row.names =
        NULL
    )
  
  
  elapsed_time <-
    proc.time()[3] -
    start_time
  
  
  list(
    method =
      "Standard Bayesian MCMC",
    
    draws =
      theta_draws,
    
    h_posterior_mean =
      h_posterior_mean,
    
    summary =
      posterior_summary,
    
    acceptance_theta =
      total_theta_accept /
      n_iter,
    
    acceptance_h =
      total_h_accept /
      (
        n_iter *
          T_local
      ),
    
    final_theta_rw_sd =
      theta_rw_sd,
    
    final_h_rw_sd =
      h_rw_sd,
    
    ess =
      overall_ess,
    
    elapsed =
      elapsed_time,
    
    ess_per_second =
      overall_ess /
      elapsed_time
  )
}


############################################################
# Run standard Bayesian MCMC
############################################################

set.seed(456)

mcmc_result <-
  mcmc_standard_bayes(
    y_data = y,
    n_iter = 20000,
    burnin = 5000
  )


############################################################
# Generalized Bayesian eta_V values
############################################################

eta_v_values <-
  c(
    0.50,
    0.75,
    1.00
  )


############################################################
# Run PBIS for each eta_V
#
# For initial testing use:
#
# B_boot = 50
# B_is   = 1000
#
# Increase after validating implementation.
############################################################

pbis_results <-
  vector(
    "list",
    length(
      eta_v_values
    )
  )


names(
  pbis_results
) <-
  paste0(
    "eta_v_",
    eta_v_values
  )


for (
  k in seq_along(
    eta_v_values
  )
) {
  
  set.seed(
    1000 +
      k
  )
  
  
  pbis_results[[k]] <-
    pbis_sv(
      y_data = y,
      
      eta_v =
        eta_v_values[k],
      
      eta_r =
        1,
      
      B_boot =
        50,
      
      B_is =
        1000,
      
      proposal_inflation =
        1.15
    )
}


############################################################
# Posterior comparison table
############################################################

comparison_list <-
  list()


comparison_list[[1]] <-
  data.frame(
    method =
      "Standard Bayesian MCMC",
    
    eta_R =
      1,
    
    eta_V =
      1,
    
    parameter =
      mcmc_result$summary$parameter,
    
    mean =
      mcmc_result$summary$mean,
    
    sd =
      mcmc_result$summary$sd
  )


for (
  k in seq_along(
    pbis_results
  )
) {
  
  comparison_list[[
    k + 1
  ]] <-
    data.frame(
      method =
        pbis_results[[k]]$method,
      
      eta_R =
        1,
      
      eta_V =
        pbis_results[[k]]$eta_v,
      
      parameter =
        pbis_results[[k]]$summary$parameter,
      
      mean =
        pbis_results[[k]]$summary$mean,
      
      sd =
        pbis_results[[k]]$summary$sd
    )
}


posterior_comparison <-
  do.call(
    rbind,
    comparison_list
  )


############################################################
# Computational cost comparison
############################################################

cost_table <-
  data.frame(
    method =
      c(
        mcmc_result$method,
        vapply(
          pbis_results,
          function(x) {
            x$method
          },
          character(1)
        )
      ),
    
    elapsed_seconds =
      c(
        mcmc_result$elapsed,
        vapply(
          pbis_results,
          function(x) {
            x$elapsed
          },
          numeric(1)
        )
      ),
    
    ESS =
      c(
        mcmc_result$ess,
        vapply(
          pbis_results,
          function(x) {
            x$ess
          },
          numeric(1)
        )
      ),
    
    ESS_per_second =
      c(
        mcmc_result$ess_per_second,
        vapply(
          pbis_results,
          function(x) {
            x$ess_per_second
          },
          numeric(1)
        )
      )
  )


############################################################
# PBIS diagnostics
############################################################

pbis_diagnostics <-
  data.frame(
    eta_V =
      eta_v_values,
    
    ESS =
      vapply(
        pbis_results,
        function(x) {
          x$ess
        },
        numeric(1)
      ),
    
    relative_ESS =
      vapply(
        pbis_results,
        function(x) {
          x$relative_ess
        },
        numeric(1)
      ),
    
    maximum_weight =
      vapply(
        pbis_results,
        function(x) {
          x$maximum_weight
        },
        numeric(1)
      ),
    
    bootstrap_successful =
      vapply(
        pbis_results,
        function(x) {
          x$bootstrap_successful
        },
        numeric(1)
      ),
    
    bootstrap_attempts =
      vapply(
        pbis_results,
        function(x) {
          x$bootstrap_attempts
        },
        numeric(1)
      )
  )


############################################################
# Plot simulated returns
############################################################

plot(
  y,
  type = "l",
  xlab = "Time",
  ylab = "Return",
  main =
    "Simulated stochastic volatility returns"
)


############################################################
# Plot true latent log volatility
############################################################

plot(
  h_true,
  type = "l",
  lwd = 2,
  xlab = "Time",
  ylab =
    expression(
      h[t]
    ),
  main =
    "True latent log volatility"
)


############################################################
# Compare true h and standard Bayesian posterior mean h
############################################################

plot(
  h_true,
  type = "l",
  lwd = 2,
  xlab = "Time",
  ylab =
    expression(
      h[t]
    ),
  main =
    "True and MCMC posterior mean log volatility"
)


lines(
  mcmc_result$h_posterior_mean,
  col = "red",
  lwd = 2,
  lty = 2
)


legend(
  "topright",
  legend =
    c(
      "True log volatility",
      "MCMC posterior mean"
    ),
  col =
    c(
      "black",
      "red"
    ),
  lwd =
    2,
  lty =
    c(
      1,
      2
    ),
  bty =
    "n"
)


############################################################
# Posterior density comparison
#
# MCMC versus PBIS for each eta_V
############################################################

parameter_labels <-
  expression(
    mu,
    phi,
    sigma
  )


true_values <-
  c(
    mu_true,
    phi_true,
    sigma_true
  )


plot_colours <-
  c(
    "blue",
    "darkgreen",
    "purple"
  )


old_par <-
  par(
    no.readonly = TRUE
  )


par(
  mar =
    c(
      5,
      5,
      4,
      10
    )
)


for (j in seq_len(3)) {
  
  mcmc_draws_j <-
    mcmc_result$draws[, j]
  
  
  all_values <-
    mcmc_draws_j
  
  
  for (
    k in seq_along(
      pbis_results
    )
  ) {
    
    all_values <-
      c(
        all_values,
        pbis_results[[k]]$draws[, j]
      )
  }
  
  
  plot_range <-
    range(
      all_values,
      finite = TRUE
    )
  
  
  mcmc_density <-
    density(
      mcmc_draws_j,
      n = 500,
      from =
        plot_range[1],
      to =
        plot_range[2]
    )
  
  
  pbis_densities <-
    vector(
      "list",
      length(
        pbis_results
      )
    )
  
  
  max_density <-
    max(
      mcmc_density$y
    )
  
  
  for (
    k in seq_along(
      pbis_results
    )
  ) {
    
    pbis_densities[[k]] <-
      density(
        x =
          pbis_results[[k]]$draws[, j],
        
        weights =
          pbis_results[[k]]$weights,
        
        n =
          500,
        
        from =
          plot_range[1],
        
        to =
          plot_range[2]
      )
    
    
    max_density <-
      max(
        max_density,
        pbis_densities[[k]]$y
      )
  }
  
  
  plot(
    mcmc_density$x,
    mcmc_density$y,
    
    type =
      "l",
    
    lwd =
      2,
    
    col =
      "red",
    
    xlab =
      parameter_labels[j],
    
    ylab =
      "Posterior density",
    
    main =
      bquote(
        "Posterior comparison for " *
          .(
            parameter_labels[[j]]
          )
      ),
    
    xlim =
      plot_range,
    
    ylim =
      c(
        0,
        1.05 *
          max_density
      )
  )
  
  
  for (
    k in seq_along(
      pbis_results
    )
  ) {
    
    lines(
      pbis_densities[[k]]$x,
      pbis_densities[[k]]$y,
      
      col =
        plot_colours[k],
      
      lwd =
        2,
      
      lty =
        k + 1
    )
  }
  
  
  abline(
    v =
      true_values[j],
    
    col =
      "black",
    
    lwd =
      2,
    
    lty =
      3
  )
  
  
  legend(
    "topright",
    
    inset =
      c(
        -1.15,
        0
      ),
    
    legend =
      c(
        "Standard Bayes MCMC",
        
        paste0(
          "PBIS eta_V = ",
          eta_v_values
        ),
        
        "True value"
      ),
    
    col =
      c(
        "red",
        plot_colours,
        "black"
      ),
    
    lwd =
      2,
    
    lty =
      c(
        1,
        seq_along(
          eta_v_values
        ) +
          1,
        3
      ),
    
    bty =
      "n",
    
    xpd =
      TRUE
  )
}


par(
  old_par
)


############################################################
# Output
############################################################

cat(
  "\nTrue stochastic volatility parameters:\n"
)

print(
  theta_true
)


cat(
  "\nStandard Bayesian MCMC posterior:\n"
)

print(
  mcmc_result$summary
)


cat(
  "\nStandard Bayesian MCMC acceptance rates:\n"
)

print(
  mcmc_result$acceptance_theta
)

cat(
  "\nLatent volatility acceptance rate:\n"
)

print(
  mcmc_result$acceptance_h
)


############################################################
# PBIS summaries
############################################################

for (
  k in seq_along(
    pbis_results
  )
) {
  
  cat(
    "\nGeneralized Bayesian PBIS posterior: eta_V = ",
    pbis_results[[k]]$eta_v,
    "\n",
    sep = ""
  )
  
  
  print(
    pbis_results[[k]]$summary
  )
  
  
  cat(
    "\nPseudo likelihood estimator:\n"
  )
  
  
  print(
    pbis_results[[k]]$theta_hat
  )
  
  
  cat(
    "\nBootstrap estimator mean:\n"
  )
  
  
  print(
    pbis_results[[k]]$bootstrap_mean
  )
  
  
  cat(
    "\nLaplace covariance:\n"
  )
  
  
  print(
    pbis_results[[k]]$laplace_covariance
  )
}


cat(
  "\nPosterior comparison:\n"
)

print(
  posterior_comparison
)


cat(
  "\nPBIS diagnostics:\n"
)

print(
  pbis_diagnostics
)


cat(
  "\nComputational cost comparison:\n"
)

print(
  cost_table
)