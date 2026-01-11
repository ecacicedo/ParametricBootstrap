library(MASS)

mech <- c(7, 44, 49, 59, 34, 46, 0, 32, 49, 52, 44, 
          36, 42, 5, 22, 18, 41, 48, 31, 42, 46, 63)
vec <- c(51, 69, 41, 70, 42, 40, 40, 45, 57, 64, 61,
         59, 60, 30, 58, 51, 63, 38, 42, 69, 49, 63)

N <- length(mech)

# Observed correlation
obs_cor <- cor(mech, vec)

# Fit bivariate normal
mu_hat <- colMeans(cbind(mech,vec))
Sigma_hat <- cov(cbind(mech, vec))

# Bootstrap
B <- 10000
boot_cor <- numeric(B)

set.seed(1)
for (i in 1:B) {
  samp <- mvrnorm(n=22, mu=mu_hat, Sigma=Sigma_hat)
  boot_cor[i] <- cor(samp[,1], samp[,2])
}

# Raw bootstrap density
dens_boot <- density(boot_cor)

# Jeffrey's prior weights
jeff_prior <- function(theta) 1/(1-theta^2)

# Fisher transformation z ~ Normal with mean = fisher_z and sd = fisher_sd
fisher_mu <- function(r) 0.5*log((1+r)/(1-r))
fisher_sd <- sd=1/sqrt(N-3)
z_hat <- fisher_mu(obs_cor)
z_i <- fisher_mu(boot_cor)

# Likelihood ratio weight
R_i <- dnorm(z_hat, mean=z_i, sd=fisher_sd) / dnorm(z_i, mean=z_hat, sd=fisher_sd)

w_raw <- jeff_prior(boot_cor) * R_i
w_norm <- w_raw/sum(w_raw)

# Posterior density
bw0 <- density(boot_cor)$bw
dens_post <- density(boot_cor, weights = w_norm, bw = bw0)

# Plot
plot(dens_boot, col="grey", main="Raw bootstrap density v Jeffrey's posterior correlation")
lines(dens_post, col="black")
legend("topright", legend=c("Bootstrap","Jeffrey's posterior"),
       col=c("grey","black"), lty=1)
