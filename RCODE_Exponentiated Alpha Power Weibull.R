# ==============================================================================
# EAPW DISTRIBUTION: MULTI-DATASET REAL BENCHMARKS & 5-SET MONTE CARLO SIMULATION
# Datasets: 'faithful$eruptions', 'iris$Petal.Length', 'iris$Petal.Width'
# Simulation Focus: EAPW Parameter Recovery Across 5 Distinct Parameter Vectors
# Optimization: Parallel Processing Engine + Integrated Section 2.6 Profile Likelihood
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. Dependencies Setup
# ------------------------------------------------------------------------------
required_packages <- c("ggplot2", "knitr", "dplyr", "tidyr", "parallel", "gridExtra")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

suppressPackageStartupMessages({
  library(ggplot2)
  library(knitr)
  library(dplyr)
  library(tidyr)
  library(parallel)
  library(gridExtra)
})

set.seed(2026)

# ------------------------------------------------------------------------------
# 1. Distribution Functions (CDF, PDF, Hazard, Quantile Generator)
# ------------------------------------------------------------------------------

# Base Weibull CDF & PDF wrappers
p_wei <- function(x, par) pweibull(x, shape = par[1], scale = par[2])
d_wei <- function(x, par) dweibull(x, shape = par[1], scale = par[2])

# Competitors (for Real Data Benchmark)
p_ew  <- function(x, par) (pweibull(x, par[1], par[2]))^par[3]
d_ew  <- function(x, par) {
  F0 <- pmax(pmin(pweibull(x, par[1], par[2]), 1 - 1e-12), 1e-12)
  par[3] * (F0^(par[3] - 1)) * dweibull(x, par[1], par[2])
}

p_apw <- function(x, par) {
  F0 <- pweibull(x, par[1], par[2])
  if (abs(par[3] - 1) < 1e-5) return(F0)
  (par[3]^F0 - 1) / (par[3] - 1)
}
d_apw <- function(x, par) {
  F0 <- pmax(pmin(pweibull(x, par[1], par[2]), 1 - 1e-12), 1e-12)
  f0 <- dweibull(x, par[1], par[2])
  if (abs(par[3] - 1) < 1e-5) return(f0)
  (log(par[3]) / (par[3] - 1)) * (par[3]^F0) * f0
}

# --- Unified Exponentiated Alpha Power Weibull (EAPW) Core Functions ---
parse_eapw_par <- function(par, k, lam, alp, th) {
  if (!is.null(par)) {
    return(list(k = par[1], lam = par[2], alp = par[3], th = par[4]))
  }
  return(list(k = k, lam = lam, alp = alp, th = th))
}

p_eapw <- function(x, par = NULL, k = NULL, lam = NULL, alp = NULL, th = NULL) {
  p <- parse_eapw_par(par, k, lam, alp, th)
  F0 <- pweibull(x, shape = p$k, scale = p$lam)
  if (abs(p$alp - 1) < 1e-5) return(F0^p$th)
  G <- (p$alp^F0 - 1) / (p$alp - 1)
  return(G^p$th)
}

d_eapw <- function(x, par = NULL, k = NULL, lam = NULL, alp = NULL, th = NULL) {
  p <- parse_eapw_par(par, k, lam, alp, th)
  F0 <- pmax(pmin(pweibull(x, shape = p$k, scale = p$lam), 1 - 1e-12), 1e-12)
  f0 <- dweibull(x, shape = p$k, scale = p$lam)
  if (abs(p$alp - 1) < 1e-5) return(p$th * (F0^(p$th - 1)) * f0)
  G <- pmax(pmin((p$alp^F0 - 1) / (p$alp - 1), 1 - 1e-12), 1e-12)
  return((p$th * log(p$alp) / (p$alp - 1)) * (G^(p$th - 1)) * (p$alp^F0) * f0)
}

h_eapw <- function(x, par = NULL, k = NULL, lam = NULL, alp = NULL, th = NULL) {
  pdf_val  <- d_eapw(x, par = par, k = k, lam = lam, alp = alp, th = th)
  cdf_val  <- p_eapw(x, par = par, k = k, lam = lam, alp = alp, th = th)
  surv_val <- pmax(1 - cdf_val, 1e-12)
  return(pdf_val / surv_val)
}

q_eapw <- function(p, par = NULL, k = NULL, lam = NULL, alp = NULL, th = NULL) {
  pars <- parse_eapw_par(par, k, lam, alp, th)
  p <- pmax(pmin(p, 1 - 1e-12), 1e-12)
  if (abs(pars$alp - 1) < 1e-5) {
    v <- p^(1/pars$th)
  } else {
    v <- log(1 + (p^(1/pars$th)) * (pars$alp - 1)) / log(pars$alp)
  }
  v <- pmax(pmin(v, 1 - 1e-12), 1e-12)
  return(qweibull(v, shape = pars$k, scale = pars$lam))
}

r_eapw <- function(n, par) {
  k <- par[1]; lam <- par[2]; alp <- par[3]; th <- par[4]
  u <- runif(n)
  if (abs(alp - 1) < 1e-5) {
    p <- u^(1/th)
  } else {
    p <- log(1 + (u^(1/th)) * (alp - 1)) / log(alp)
  }
  p <- pmax(pmin(p, 1 - 1e-12), 1e-12)
  return(qweibull(p, shape = k, scale = lam))
}

# Negative Log-Likelihood for Direct 4D EAPW Optimization
nll_eapw_single <- function(par, data) {
  if (any(par <= 1e-4) || abs(par[3] - 1) < 1e-4) return(1e10)
  dens <- d_eapw(data, par = par)
  if (any(dens <= 0) || any(is.na(dens))) return(1e10)
  return(-sum(log(pmax(dens, 1e-300))))
}

# Section 2.6 Compliant Profile Log-Likelihood Function (Analytically Reduced Theta)
profile_nll_eapw <- function(par_3d, data) {
  k   <- par_3d[1]
  lam <- par_3d[2]
  alp <- par_3d[3]
  
  if (k <= 1e-4 || lam <= 1e-4 || alp <= 1e-4 || abs(alp - 1) < 1e-4) return(1e10)
  
  n  <- length(data)
  F0 <- pweibull(data, shape = k, scale = lam)
  G  <- pmax(pmin((alp^F0 - 1) / (alp - 1), 1 - 1e-12), 1e-12)
  
  theta_hat <- -n / sum(log(G))
  if (is.na(theta_hat) || theta_hat <= 0) return(1e10)
  
  f0 <- dweibull(data, shape = k, scale = lam)
  log_dens <- log(theta_hat) + log(log(alp)) - log(alp - 1) + 
    (theta_hat - 1) * log(G) + F0 * log(alp) + log(pmax(f0, 1e-300))
  
  return(-sum(log_dens))
}

# ==============================================================================
# 2. PART 1: REAL DATA EVALUATION FRAMEWORK
# ==============================================================================
fit_real_data <- function(data, data_name = "Dataset") {
  n <- length(data)
  
  nll_wei  <- function(p) { if(any(p <= 1e-4)) return(1e10); -sum(log(pmax(d_wei(data, p), 1e-300))) }
  nll_ew   <- function(p) { if(any(p <= 1e-4)) return(1e10); -sum(log(pmax(d_ew(data, p), 1e-300))) }
  nll_apw  <- function(p) { if(any(p <= 1e-4)) return(1e10); -sum(log(pmax(d_apw(data, p), 1e-300))) }
  nll_eapw <- function(p) { nll_eapw_single(p, data) }
  
  # Weibull baseline
  opt_wei <- optim(c(1.2, mean(data)), nll_wei, method = "Nelder-Mead")
  
  # EW Multi-Start
  inits_ew <- list(c(opt_wei$par, 1.0), c(opt_wei$par, 0.3), c(opt_wei$par, 2.5))
  best_ew_nll <- Inf; best_ew_par <- c(opt_wei$par, 1.0)
  for (init in inits_ew) {
    res <- tryCatch(optim(init, nll_ew, method = "Nelder-Mead"), error = function(e) NULL)
    if (!is.null(res) && res$value < best_ew_nll) { best_ew_nll <- res$value; best_ew_par <- res$par }
  }
  
  # APW Multi-Start
  inits_apw <- list(c(opt_wei$par, 2.0), c(opt_wei$par, 0.1), c(opt_wei$par, 12.0))
  best_apw_nll <- Inf; best_apw_par <- c(opt_wei$par, 2.0)
  for (init in inits_apw) {
    res <- tryCatch(optim(init, nll_apw, method = "Nelder-Mead"), error = function(e) NULL)
    if (!is.null(res) && res$value < best_apw_nll) { best_apw_nll <- res$value; best_apw_par <- res$par }
  }
  
  # EAPW Multi-Start
  inits_eapw <- list(
    c(best_ew_par[1], best_ew_par[2], 1.001, best_ew_par[3]),
    c(best_apw_par[1], best_apw_par[2], best_apw_par[3], 1.0),
    c(opt_wei$par, 10.0, 2.0),
    c(opt_wei$par, 30.0, 3.5)
  )
  best_eapw_nll <- Inf; best_eapw_par <- c(opt_wei$par, 10.0, 2.0)
  for (init in inits_eapw) {
    res <- tryCatch(optim(init, nll_eapw, method = "Nelder-Mead", control = list(maxit = 4000)), error = function(e) NULL)
    if (!is.null(res) && res$value < best_eapw_nll) { best_eapw_nll <- res$value; best_eapw_par <- res$par }
  }
  
  models <- list(
    "Weibull" = list(par = opt_wei$par, nll = opt_wei$value, n_par = 2, pfun = p_wei, dfun = d_wei),
    "EW"      = list(par = best_ew_par, nll = best_ew_nll, n_par = 3, pfun = p_ew, dfun = d_ew),
    "APW"     = list(par = best_apw_par, nll = best_apw_nll, n_par = 3, pfun = p_apw, dfun = d_apw),
    "EAPW"    = list(par = best_eapw_par, nll = best_eapw_nll, n_par = 4, pfun = p_eapw, dfun = d_eapw)
  )
  
  summary_df <- data.frame()
  s_data <- sort(data)
  i <- 1:n
  
  for (m_name in names(models)) {
    m <- models[[m_name]]
    aic <- 2 * m$n_par + 2 * m$nll
    bic <- log(n) * m$n_par + 2 * m$nll
    
    z <- pmax(pmin(m$pfun(s_data, m$par), 1 - 1e-12), 1e-12)
    ks_stat <- suppressWarnings(as.numeric(ks.test(data, function(x) m$pfun(x, m$par))$statistic))
    w_star  <- sum((z - (2*i - 1)/(2*n))^2) + 1/(12*n)
    a_star  <- -n - (1/n) * sum((2*i - 1)*log(z) + (2*n + 1 - 2*i)*log(1 - z))
    
    summary_df <- rbind(summary_df, data.frame(
      Dataset = data_name, Model = m_name, NegLogLik = m$nll,
      AIC = aic, BIC = bic, KS = ks_stat, W_star = w_star, A_star = a_star,
      MLE_Estimates = paste(round(m$par, 4), collapse = ", ")
    ))
  }
  
  return(list(summary = summary_df, models = models))
}

# --- Density Plot Function ---
plot_real_densities <- function(data, fit_obj, title_text = "Fitted Densities") {
  x_seq <- seq(min(data) * 0.8, max(data) * 1.15, length.out = 500)
  m <- fit_obj$models
  
  df_densities <- data.frame(
    x = x_seq,
    Weibull = d_wei(x_seq, m$Weibull$par),
    EW      = d_ew(x_seq, m$EW$par),
    APW     = d_apw(x_seq, m$APW$par),
    EAPW    = d_eapw(x_seq, par = m$EAPW$par)
  )
  
  df_long <- pivot_longer(df_densities, cols = -x, names_to = "Model", values_to = "Density")
  
  p <- ggplot() +
    geom_histogram(data = data.frame(x = data), aes(x = x, y = after_stat(density)),
                   bins = 25, fill = "slategray2", color = "white", alpha = 0.7) +
    geom_line(data = df_long, aes(x = x, y = Density, color = Model, linetype = Model), linewidth = 1.1) +
    labs(title = title_text, x = "Observed Value", y = "Density") +
    theme_minimal(base_size = 12) +
    scale_color_manual(values = c("Weibull" = "#E41A1C", "EW" = "#377EB8", "APW" = "#4DAF4A", "EAPW" = "#984EA3")) +
    theme(legend.position = "top", plot.title = element_text(face = "bold", hjust = 0.5))
  
  return(p)
}

# ------------------------------------------------------------------------------
# 2a. BENCHMARK DATASET 1: faithful$eruptions (n = 272)
# ------------------------------------------------------------------------------
data("faithful", package = "datasets")
real_fit_faithful <- fit_real_data(faithful$eruptions, "faithful$eruptions (n=272)")

cat("\n======================================================================\n")
cat("        REAL DATASET BENCHMARK 1: faithful$eruptions                  \n")
cat("======================================================================\n\n")
print(knitr::kable(real_fit_faithful$summary, digits = 4, format = "simple"))
print(plot_real_densities(faithful$eruptions, real_fit_faithful, "Fitted Densities on faithful$eruptions Data"))

# ------------------------------------------------------------------------------
# 2b. BENCHMARK DATASET 2: iris$Petal.Length (n = 150)
# ------------------------------------------------------------------------------
data("iris", package = "datasets")
real_fit_petal_length <- fit_real_data(iris$Petal.Length, "iris$Petal.Length (n=150)")

cat("\n======================================================================\n")
cat("        REAL DATASET BENCHMARK 2: iris$Petal.Length (n=150)           \n")
cat("======================================================================\n\n")
print(knitr::kable(real_fit_petal_length$summary, digits = 4, format = "simple"))
print(plot_real_densities(iris$Petal.Length, real_fit_petal_length, "Fitted Densities on iris$Petal.Length Data"))

# ------------------------------------------------------------------------------
# 2c. BENCHMARK DATASET 3: iris$Petal.Width (n = 150)
# ------------------------------------------------------------------------------
real_fit_petal_width <- fit_real_data(iris$Petal.Width, "iris$Petal.Width (n=150)")

cat("\n======================================================================\n")
cat("        REAL DATASET BENCHMARK 3: iris$Petal.Width (n=150)            \n")
cat("======================================================================\n\n")
print(knitr::kable(real_fit_petal_width$summary, digits = 4, format = "simple"))
print(plot_real_densities(iris$Petal.Width, real_fit_petal_width, "Fitted Densities on iris$Petal.Width Data"))

# ==============================================================================
# 3. PART 2: PARALLELIZED HIGH-SPEED MONTE CARLO SIMULATION (5 PARAMETER SETS)
# Evaluate MLE Performance: Mean, Bias, Relative Bias %, and RMSE
# Across Sample Sizes (n), Replications (M), and 5 Parameter Vectors
# ==============================================================================

param_sets <- list(
  "Set_1" = c(k = 1.5, lam = 2.0, alp = 5.0,  th = 1.5),
  "Set_2" = c(k = 2.0, lam = 1.5, alp = 15.0, th = 2.5),
  "Set_3" = c(k = 1.0, lam = 3.0, alp = 2.0,  th = 0.8),
  "Set_4" = c(k = 2.5, lam = 1.0, alp = 8.0,  th = 3.0),
  "Set_5" = c(k = 0.8, lam = 4.0, alp = 25.0, th = 1.2)
)

sample_sizes <- c(30, 50, 100, 200, 500) 
M_replications <- 500                    

run_eapw_simulation <- function(p_sets, n_vec, M) {
  num_cores <- max(1, detectCores() - 1)
  cl <- makeCluster(num_cores)
  
  # Export essential functions to worker nodes
  clusterExport(cl, c("r_eapw", "d_eapw", "p_eapw", "parse_eapw_par", "nll_eapw_single"))
  clusterEvalQ(cl, { set.seed(2026) })
  
  cat("\n======================================================================\n")
  cat(sprintf("   STARTING PARALLEL EAPW MONTE CARLO SIMULATION (%d CPU Cores)   \n", num_cores))
  cat("======================================================================\n\n")
  
  results_accumulator <- list()
  param_names <- c("k (shape)", "lambda (scale)", "alpha (power)", "theta (exp)")
  
  for (setName in names(p_sets)) {
    true_par <- p_sets[[setName]]
    cat(sprintf("-> Processing %s: [k=%.2f, lambda=%.2f, alpha=%.2f, theta=%.2f]\n", 
                setName, true_par[1], true_par[2], true_par[3], true_par[4]))
    
    for (n in n_vec) {
      cat(sprintf("   * Sample Size n = %d (%d replications in parallel)... ", n, M))
      
      # Parallelized replication mapping
      estimates_list <- parLapply(cl, 1:M, function(r) {
        sim_data <- r_eapw(n, true_par)
        
        init_guesses <- list(
          true_par,
          true_par * c(1.1, 0.9, 1.2, 0.8),
          true_par * c(0.8, 1.1, 0.7, 1.3)
        )
        
        best_nll <- Inf
        best_par <- NULL
        
        for (init in init_guesses) {
          fit <- tryCatch({
            optim(par = init, 
                  fn = nll_eapw_single, 
                  data = sim_data, 
                  method = "L-BFGS-B",
                  lower = c(0.001, 0.001, 0.001, 0.001),
                  upper = c(100, 100, 200, 100))
          }, error = function(e) NULL)
          
          if (!is.null(fit) && fit$convergence == 0 && fit$value < best_nll) {
            best_nll <- fit$value
            best_par <- fit$par
          }
        }
        return(best_par)
      })
      
      # Aggregate non-NULL parameter fits
      estimates_mat <- do.call(rbind, estimates_list[!sapply(estimates_list, is.null)])
      valid_count <- ifelse(!is.null(estimates_mat), nrow(estimates_mat), 0)
      
      if (!is.null(estimates_mat) && valid_count > 0) {
        for (j in 1:4) {
          par_true <- true_par[j]
          par_sim  <- estimates_mat[, j]
          
          mean_est <- mean(par_sim)
          bias     <- mean_est - par_true
          rbias    <- (bias / par_true) * 100
          rmse     <- sqrt(mean((par_sim - par_true)^2))
          
          results_accumulator[[length(results_accumulator) + 1]] <- data.frame(
            Parameter_Set = setName,
            Sample_Size   = n,
            Parameter     = param_names[j],
            True_Value    = par_true,
            Mean_MLE      = round(mean_est, 4),
            Bias          = round(bias, 4),
            RBias_Pct     = round(rbias, 2),
            RMSE          = round(rmse, 4),
            Converged_Reps= valid_count
          )
        }
      }
      cat("Done.\n")
    }
  }
  
  stopCluster(cl)
  return(do.call(rbind, results_accumulator))
}

# Run optimized parallel Monte Carlo Simulation across all 5 parameter sets
sim_summary <- run_eapw_simulation(param_sets, sample_sizes, M_replications)

cat("\n======================================================================\n")
cat("              EAPW MONTE CARLO SIMULATION RESULTS TABLE               \n")
cat("======================================================================\n\n")
print(knitr::kable(sim_summary, digits = 4, format = "simple"))

# ==============================================================================
# VISUALIZING THEORETICAL DENSITY CURVES FOR 5 EAPW PARAMETER SETS
# ==============================================================================

# Define Named Parameter Vectors for Density Comparison
density_param_sets <- list(
  "Set 1: Baseline (k=1.5, λ=2.0, α=5, θ=1.5)"   = c(k = 1.5, lam = 2.0, alp = 5.0,  th = 1.5),
  "Set 2: High Alpha (k=2.0, λ=1.5, α=15, θ=2.5)" = c(k = 2.0, lam = 1.5, alp = 15.0, th = 2.5),
  "Set 3: Exponential Base (k=1.0, λ=3.0, α=2, θ=0.8)" = c(k = 1.0, lam = 3.0, alp = 2.0,  th = 0.8),
  "Set 4: High Shape/Theta (k=2.5, λ=1.0, α=8, θ=3.0)" = c(k = 2.5, lam = 1.0, alp = 8.0,  th = 3.0),
  "Set 5: Heavy-Tailed (k=0.8, λ=4.0, α=25, θ=1.2)"   = c(k = 0.8, lam = 4.0, alp = 25.0, th = 1.2)
)

x_seq <- seq(0.01, 8.0, length.out = 1000)

plot_data <- bind_rows(lapply(names(density_param_sets), function(setName) {
  p <- density_param_sets[[setName]]
  dens <- d_eapw(x_seq, par = p)
  data.frame(x = x_seq, density = dens, Configuration = setName)
}))

# Plot 1: Overlaid Density Comparison
p_overlay <- ggplot(plot_data, aes(x = x, y = density, color = Configuration, linetype = Configuration)) +
  geom_line(linewidth = 1.1) +
  coord_cartesian(ylim = c(0, 1.25), xlim = c(0, 7.5)) +
  labs(
    title = "EAPW Theoretical Density Curves (Overlaid Comparison)",
    subtitle = "Visualizing structural variations across parameter spaces",
    x = "x",
    y = "Density f(x)"
  ) +
  scale_color_brewer(palette = "Set1") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.direction = "vertical",
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

# Plot 2: Side-by-Side Faceted Grid
p_facet <- ggplot(plot_data, aes(x = x, y = density, fill = Configuration)) +
  geom_area(alpha = 0.4, color = "black", linewidth = 0.4) +
  facet_wrap(~ Configuration, scales = "free_y", ncol = 3) +
  coord_cartesian(xlim = c(0, 7.5)) +
  labs(
    title = "EAPW Individual Density Profiles Side-by-Side",
    x = "x",
    y = "Density f(x)"
  ) +
  scale_fill_brewer(palette = "Set1") +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5),
    strip.background = element_rect(fill = "grey92"),
    strip.text = element_text(face = "bold", size = 8.5)
  )

# Render Plots
print(p_overlay)
print(p_facet)

# ==============================================================================
# EXTENDED DIAGNOSTIC FRAMEWORK: ADVANCED TABLES & COLORFUL PLOTS FOR EAPW
# Features:
#   1. Theoretical Hazard Rate Function (HRF) Visualizer (Bathtub, Unimodal, etc.)
#   2. Standard Errors & 95% Confidence Interval Table via Hessian Matrix
#   3. Comprehensive Model Ranking & Fitness Metric Summary Table
#   4. Empirical CDF vs. Fitted Model CDF Comparison Plots
#   5. Diagnostic P-P (Probability-Probability) & Q-Q Plots
#   6. Monte Carlo RMSE & Bias Convergence Trend Line Plots
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Theoretical Hazard Rate Function (HRF) Visualizer
# ------------------------------------------------------------------------------

hrf_param_sets <- list(
  "Increasing"   = c(k = 2.5, lam = 1.5, alp = 5.0,  th = 2.0),
  "Decreasing"   = c(k = 0.6, lam = 2.0, alp = 0.5,  th = 0.8),
  "Bathtub"      = c(k = 0.5, lam = 3.0, alp = 40.0, th = 0.4),
  "Unimodal"     = c(k = 1.8, lam = 1.0, alp = 15.0, th = 0.5),
  "Upside-Down" = c(k = 1.2, lam = 2.5, alp = 8.0,  th = 0.3)
)

x_hrf_seq <- seq(0.05, 5.0, length.out = 800)
hrf_plot_df <- bind_rows(lapply(names(hrf_param_sets), function(shape_name) {
  p <- hrf_param_sets[[shape_name]]
  hrf_vals <- h_eapw(x_hrf_seq, par = p)
  data.frame(x = x_hrf_seq, Hazard = hrf_vals, Shape = shape_name)
}))

p_hrf <- ggplot(hrf_plot_df, aes(x = x, y = Hazard, color = Shape)) +
  geom_line(linewidth = 1.2) +
  coord_cartesian(ylim = c(0, 4.0), xlim = c(0, 5.0)) +
  scale_color_viridis_d(option = "plasma", end = 0.9) +
  labs(
    title = "EAPW Theoretical Hazard Rate Function (HRF) Flexibility",
    subtitle = "Demonstrating Monotonic, Bathtub, Unimodal, and Upside-Down Shapes",
    x = "Time / Value (x)",
    y = "Hazard Rate h(x)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

print(p_hrf)

# ------------------------------------------------------------------------------
# 2. Table: Real Data Parameter Estimates, Standard Errors, & 95% CIs
# ------------------------------------------------------------------------------

compute_eapw_se <- function(data, data_name, warm_start) {
  nll <- function(p) -sum(log(d_eapw(data, par = p)))
  
  opt <- tryCatch({
    optim(warm_start, nll, hessian = TRUE, method = "L-BFGS-B",
          lower = c(1e-3, 1e-3, 1e-3, 1e-3))
  }, error = function(e) NULL)
  
  if (is.null(opt) || any(is.na(opt$hessian))) return(NULL)
  
  cov_mat <- tryCatch(solve(opt$hessian), error = function(e) matrix(NA, 4, 4))
  se <- sqrt(diag(cov_mat))
  
  par_names <- c("k (shape)", "lambda (scale)", "alpha (power)", "theta (exp)")
  estimates <- opt$par
  lower_ci  <- estimates - 1.96 * se
  upper_ci  <- estimates + 1.96 * se
  
  data.frame(
    Dataset     = data_name,
    Parameter   = par_names,
    MLE         = round(estimates, 4),
    Std_Error   = round(se, 4),
    CI_95_Lower = round(lower_ci, 4),
    CI_95_Upper = round(upper_ci, 4)
  )
}

se_table_faithful  <- compute_eapw_se(faithful$eruptions, "faithful$eruptions", real_fit_faithful$models$EAPW$par)
se_table_petal_len <- compute_eapw_se(iris$Petal.Length, "iris$Petal.Length", real_fit_petal_length$models$EAPW$par)
se_table_petal_wid <- compute_eapw_se(iris$Petal.Width, "iris$Petal.Width", real_fit_petal_width$models$EAPW$par)

all_se_table <- rbind(se_table_faithful, se_table_petal_len, se_table_petal_wid)

cat("\n======================================================================\n")
cat("      EAPW PARAMETER ESTIMATES & 95% CONFIDENCE INTERVALS TABLE       \n")
cat("======================================================================\n\n")
print(knitr::kable(all_se_table, digits = 4, format = "simple"))

# ------------------------------------------------------------------------------
# 3. Table: Comprehensive Cross-Dataset Model Ranking Summary
# ------------------------------------------------------------------------------

combined_benchmark_df <- rbind(
  real_fit_faithful$summary,
  real_fit_petal_length$summary,
  real_fit_petal_width$summary
)

ranked_summary_df <- combined_benchmark_df %>%
  group_by(Dataset) %>%
  mutate(
    Rank_AIC    = min_rank(AIC),
    Rank_BIC    = min_rank(BIC),
    Rank_KS     = min_rank(KS),
    Rank_A_Star = min_rank(A_star)
  ) %>%
  ungroup() %>%
  select(Dataset, Model, NegLogLik, AIC, Rank_AIC, BIC, Rank_BIC, KS, Rank_KS, A_star, Rank_A_Star)

cat("\n======================================================================\n")
cat("         CROSS-DATASET MODEL PERFORMANCE & RANKING SUMMARY            \n")
cat("======================================================================\n\n")
print(knitr::kable(ranked_summary_df, digits = 4, format = "simple"))

# ------------------------------------------------------------------------------
# 4. Plots: Empirical CDF vs. Fitted Model CDF Comparison
# ------------------------------------------------------------------------------

plot_cdf_comparison <- function(data, fit_obj, title_text) {
  x_seq <- seq(min(data) * 0.8, max(data) * 1.15, length.out = 400)
  m     <- fit_obj$models
  
  cdf_df <- data.frame(
    x       = x_seq,
    Weibull = p_wei(x_seq, m$Weibull$par),
    EW      = p_ew(x_seq, m$EW$par),
    APW     = p_apw(x_seq, m$APW$par),
    EAPW    = p_eapw(x_seq, par = m$EAPW$par)
  ) %>%
    pivot_longer(cols = -x, names_to = "Model", values_to = "CDF")
  
  ggplot() +
    stat_ecdf(data = data.frame(x = data), aes(x = x), color = "black", linewidth = 0.8, linetype = "dashed") +
    geom_line(data = cdf_df, aes(x = x, y = CDF, color = Model), linewidth = 1.1) +
    scale_color_manual(values = c("Weibull" = "#E41A1C", "EW" = "#377EB8", "APW" = "#4DAF4A", "EAPW" = "#984EA3")) +
    labs(title = title_text, subtitle = "Black dashed line indicates Empirical CDF", x = "Value", y = "Cumulative Probability") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), legend.position = "bottom")
}

p_cdf_1 <- plot_cdf_comparison(faithful$eruptions, real_fit_faithful, "CDF Comparison: faithful$eruptions")
p_cdf_2 <- plot_cdf_comparison(iris$Petal.Length, real_fit_petal_length, "CDF Comparison: iris$Petal.Length")
p_cdf_3 <- plot_cdf_comparison(iris$Petal.Width, real_fit_petal_width, "CDF Comparison: iris$Petal.Width")

grid.arrange(p_cdf_1, p_cdf_2, p_cdf_3, ncol = 3)

# ------------------------------------------------------------------------------
# 5. Plots: Diagnostic P-P and Q-Q Plots for EAPW Model
# ------------------------------------------------------------------------------

generate_pp_qq_plots <- function(data, par, title_suffix) {
  s_data <- sort(data)
  n      <- length(s_data)
  p_emp  <- (1:n - 0.5) / n
  
  p_theo <- p_eapw(s_data, par = par)
  q_theo <- q_eapw(p_emp, par = par)
  
  df_diag <- data.frame(
    Observed_Q    = s_data, 
    Theoretical_Q = q_theo,
    Empirical_P   = p_emp, 
    Theoretical_P = p_theo
  )
  
  p_pp <- ggplot(df_diag, aes(x = Theoretical_P, y = Empirical_P)) +
    geom_point(color = "#984EA3", alpha = 0.7, size = 1.8) +
    geom_abline(intercept = 0, slope = 1, color = "darkred", linetype = "dashed", linewidth = 0.8) +
    labs(title = paste("P-P Plot:", title_suffix), x = "Theoretical Probability", y = "Empirical Probability") +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
  
  p_qq <- ggplot(df_diag, aes(x = Theoretical_Q, y = Observed_Q)) +
    geom_point(color = "#377EB8", alpha = 0.7, size = 1.8) +
    geom_abline(intercept = 0, slope = 1, color = "darkred", linetype = "dashed", linewidth = 0.8) +
    labs(title = paste("Q-Q Plot:", title_suffix), x = "Theoretical Quantile", y = "Observed Quantile") +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
  
  return(list(pp = p_pp, qq = p_qq))
}

# Generate diagnostics for all 3 datasets
diag_1 <- generate_pp_qq_plots(faithful$eruptions, real_fit_faithful$models$EAPW$par, "faithful$eruptions")
diag_2 <- generate_pp_qq_plots(iris$Petal.Length, real_fit_petal_length$models$EAPW$par, "iris$Petal.Length")
diag_3 <- generate_pp_qq_plots(iris$Petal.Width, real_fit_petal_width$models$EAPW$par, "iris$Petal.Width")

# Arrange in a 3x2 grid (P-P on left, Q-Q on right)
grid.arrange(
  diag_1$pp, diag_1$qq,
  diag_2$pp, diag_2$qq,
  diag_3$pp, diag_3$qq,
  ncol = 2, nrow = 3
)

# ------------------------------------------------------------------------------
# 6. Plots: Monte Carlo Simulation Convergence Trends
# ------------------------------------------------------------------------------

p_sim_rmse <- ggplot(sim_summary, aes(x = Sample_Size, y = RMSE, color = Parameter, group = Parameter)) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 2.2) +
  facet_wrap(~ Parameter_Set, scales = "free_y", ncol = 3) +
  scale_color_brewer(palette = "Set2") +
  scale_x_continuous(breaks = c(30, 50, 100, 200, 500)) +
  labs(
    title = "Monte Carlo RMSE Convergence Profiles Across Parameter Sets",
    subtitle = "Demonstrating Asymptotic Consistency as Sample Size (n) Increases",
    x = "Sample Size (n)",
    y = "Root Mean Squared Error (RMSE)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold")
  )

print(p_sim_rmse)

p_sim_rbias <- ggplot(sim_summary, aes(x = Sample_Size, y = RBias_Pct, color = Parameter, group = Parameter)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  geom_line(linewidth = 1.0) +
  geom_point(size = 2.0) +
  facet_wrap(~ Parameter_Set, scales = "free_y", ncol = 3) +
  scale_color_brewer(palette = "Dark2") +
  scale_x_continuous(breaks = c(30, 50, 100, 200, 500)) +
  labs(
    title = "Monte Carlo Relative Bias (%) Trajectories",
    subtitle = "Convergence of Parameter Estimates Toward Unbiased Zero Line",
    x = "Sample Size (n)",
    y = "Relative Bias (%)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold")
  )

print(p_sim_rbias)
