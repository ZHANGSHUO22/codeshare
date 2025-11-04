library(ggplot2)
library(coda)

set.seed(123)

# -------------------------------
# 1. Global parameters
# -------------------------------
T           <- 10000
Delta_t     <- 1 / T
mu_true     <- 0.15
sigma_true  <- 0.5  
S0          <- 100
burn        <- 1000
n_iter      <- 10000
n_paths     <- 100
alpha0      <- 2
beta0       <- 0.25

# -------------------------------
# 2. Data generation
# -------------------------------
generate_gbm <- function(T, Delta_t, mu, sigma, S0) {
  S <- numeric(T + 1); S[1] <- S0
  for (t in 1:T) {
    eps <- rnorm(1)
    S[t + 1] <- S[t] * (1 + mu * Delta_t + sigma * sqrt(Delta_t) * eps)
  }
  return(S)
}
all_paths <- lapply(seq_len(n_paths),
                    function(i) generate_gbm(T, Delta_t, mu_true, sigma_true, S0))

# -------------------------------
# 3. Marginal likelihood under EM approximation
# -------------------------------
log_marginal_em <- function(sigma, S, Delta_t) {
  if (sigma <= 0) return(-Inf)
  N    <- length(S) - 1
  St   <- S[1:N]
  Snext <- S[-1]
  x    <- St * Delta_t
  y    <- Snext - St
  w    <- 1 / (St^2 * Delta_t)
  mu_hat <- sum(w * x * y) / sum(w * x^2)
  res    <- y - mu_hat * x
  SSR    <- sum(w * res^2)
  ll <- - (N/2) * log(sigma^2) - SSR / (2 * sigma^2)
  return(ll)
}

# -------------------------------
# 4. Marginal likelihood under exact log-returns
# -------------------------------
log_marginal_exact <- function(sigma, S, Delta_t) {
  if (sigma <= 0) return(-Inf)
  R    <- diff(log(S))
  N    <- length(R)
  var  <- sigma^2 * Delta_t
  SSR  <- sum((R - mean(R))^2)
  ll   <- -(N - 1)/2 * log(var) - SSR/(2 * var)
  return(ll)
}

# -------------------------------
# 5. Posterior functions
# -------------------------------
log_prior <- function(sigma) {
  if (sigma <= 0) return(-Inf)
  sigma2 <- sigma^2
  alpha0 * log(beta0) - lgamma(alpha0) -
    (alpha0 + 1) * log(sigma2) - beta0 / sigma2 + log(2 * sigma)
}

log_posterior <- function(sigma, S_path, Delta_t, method = c("exact", "em")) {
  method <- match.arg(method)
  lp <- log_prior(sigma)
  ll <- switch(method,
               exact = log_marginal_exact(sigma, S_path, Delta_t),
               em    = log_marginal_em   (sigma, S_path, Delta_t)
  )
  lp + ll
}

# -------------------------------
# 6. Adaptive log-normal MH sampler
# -------------------------------
adaptive_lognormal_mh <- function(S_path, n_iter = 10000, burnin = 1000,
                                  init = 1.0, proposal_sd = 0.3,
                                  target_accept = 0.4, adapt_window = 1000,
                                  method = c("exact", "em")) {
  chain <- numeric(n_iter)
  accepted <- logical(n_iter)
  chain[1] <- init
  current_logpost <- log_posterior(init, S_path, Delta_t, method)
  gamma <- 0.05
  
  for (i in 2:n_iter) {
    proposal <- rlnorm(1, log(chain[i - 1]), proposal_sd)
    proposal_logpost <- log_posterior(proposal, S_path, Delta_t, method)
    
    log_q_forward <- dlnorm(proposal, log(chain[i - 1]), proposal_sd, log = TRUE)
    log_q_backward <- dlnorm(chain[i - 1], log(proposal), proposal_sd, log = TRUE)
    log_alpha <- proposal_logpost - current_logpost + log_q_backward - log_q_forward
    
    if (log(runif(1)) < log_alpha) {
      chain[i] <- proposal
      current_logpost <- proposal_logpost
      accepted[i] <- TRUE
    } else {
      chain[i] <- chain[i - 1]
      accepted[i] <- FALSE
    }
    
    if (i <= adapt_window) {
      acc_rate <- mean(accepted[max(1, (i - 100)):i])
      proposal_sd <- proposal_sd * exp(gamma * (acc_rate - target_accept))
    }
  }
  
  samples <- chain[(burnin + 1):n_iter]
  ci <- quantile(samples, c(0.025, 0.975))
  ess <- tryCatch(coda::effectiveSize(samples), error = function(e) NA)
  final_acc   <- mean(accepted[(burnin + 1):n_iter])
  list(mean = mean(samples), lower = ci[1], upper = ci[2],
       covered = sigma_true >= ci[1] & sigma_true <= ci[2],
       ess = ess,
       acc_rate = final_acc,
       chain = chain)
}

plot_diagnostics <- function(chains_list, method_label) {
  chain1 <- chains_list[[1]][(burn+1):n_iter]
  ts.plot(chain1,
          main = paste0("Trace Plot (", method_label, ")"),
          ylab = expression(sigma),
          xlab = "Iteration")
  acf(chain1, main = paste0("ACF (", method_label, ")"))
  plot(density(chain1),
       main = paste0("Posterior Density (", method_label, ")"),
       xlab = expression(sigma),
       ylab = "Density")
}

# -------------------------------
# 7. Simulation over multiple paths
# -------------------------------
run_simulation <- function(method, paths_list) {
  results <- data.frame(mean = numeric(n_paths),
                        lower = numeric(n_paths),
                        upper = numeric(n_paths),
                        covered = logical(n_paths),
                        ess = numeric(n_paths),
                        acc_rate = numeric(n_paths))
  chains_list <- vector("list", n_paths)
  
  for (i in seq_len(n_paths)) {
    S_path <- paths_list[[i]]
    res    <- adaptive_lognormal_mh(S_path,
                                    n_iter       = n_iter,
                                    burnin       = burn,
                                    method       = method)
    results[i, ] <- unlist(res[c("mean","lower","upper","covered","ess","acc_rate")])
    chains_list[[i]] <- res$chain
  }
  
  bias     <- mean(results$mean - sigma_true)
  rmse     <- sqrt(mean((results$mean - sigma_true)^2))
  coverage <- mean(results$covered)
  avg_ess  <- mean(results$ess, na.rm = TRUE)
  avg_acc  <- mean(results$acc_rate)
  
  list(results     = results,
       chains_list = chains_list,
       summary     = list(bias     = bias,
                          rmse     = rmse,
                          coverage = coverage,
                          avg_ess  = avg_ess,
                          avg_acc  = avg_acc))
}

# -------------------------------
# 8. Results
# -------------------------------
res_exact <- run_simulation("exact", all_paths)
res_em    <- run_simulation("em", all_paths)

cat("=== Exact Method ===\n")
print(res_exact$summary)
cat("\n=== Euler–Maruyama Approximation Method ===\n")
print(res_em$summary)

dev.new(width=7, height=7)
plot_diagnostics(res_exact$chains_list, "Exact")

plot_diagnostics(res_em$chains_list, "EM Approx")

df_means <- data.frame(
  method = rep(c("Exact", "EM"), each = n_paths),
  mean   = c(res_exact$results$mean, res_em$results$mean)
)
ggplot(df_means, aes(x = mean, fill = method)) +
  geom_histogram(position = "identity", alpha = 0.6, bins = 20) +
  geom_vline(aes(xintercept = sigma_true), color = "red", linetype = "dashed") +
  facet_wrap(~ method) +
  labs(title = "Posterior Means Distribution",
       x = expression(hat(sigma)), y = "Count") +
  theme_minimal()

df_ess <- data.frame(
  method = rep(c("Exact", "EM"), each = n_paths),
  ess    = c(res_exact$results$ess, res_em$results$ess)
)
ggplot(df_ess, aes(x = ess, fill = method)) +
  geom_histogram(position = "identity", alpha = 0.6, bins = 20) +
  facet_wrap(~ method) +
  labs(title = "Effective Sample Size (ESS) Distribution",
       x = "ESS", y = "Count") +
  theme_minimal()

traj_idx <- c(1, 50, 100)
density_df <- do.call(rbind, lapply(traj_idx, function(i) {
  exact_samps <- res_exact$chains_list[[i]][(burn+1):n_iter]
  em_samps    <- res_em$chains_list[[i]][(burn+1):n_iter]
  data.frame(
    sigma      = c(exact_samps, em_samps),
    method     = rep(c("Exact", "EM Approx"), each = length(exact_samps)),
    trajectory = factor(paste0("Path ", i), levels = paste0("Path ", traj_idx))
  )
}))

ggplot(density_df, aes(x = sigma, color = method)) +
  geom_density(size = 1) +
  scale_color_manual(values = c( "EM Approx" = "green", "Exact" = "blue")) +
  facet_wrap(~ trajectory, nrow = 1) +
  labs(title = "Posterior Density: Exact vs. EM Approximation",
       x = expression(sigma), y = "Density",
       color = "") +
  theme_minimal(base_size = 14)