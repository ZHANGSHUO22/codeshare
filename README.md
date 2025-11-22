# Bayesian Estimation of GBM Volatility

This repository contains the **R** implementation of the simulation study presented in the Master's Thesis: **"Bayesian Statistics for Stochastic Processes in Finance"** (Charles University, 2025).

## 📖 Overview

This project implements a Bayesian inference framework to estimate the volatility parameter ($\sigma$) of a **Geometric Brownian Motion (GBM)** model using high-frequency financial data.

The core contribution is a systematic comparison between two likelihood construction methods:

1.  **Exact Method:** Based on the exact log-return density.
2.  **Euler-Maruyama (EM) Approximation:** Based on the discretized SDE solution.

The repository features an **Adaptive Log-Normal Random-Walk Metropolis-Hastings (MH) sampler** that analytically marginalizes the drift parameter ($\mu$) to focus solely on volatility estimation.

## 🚀 Key Features

  * **Data Generation:** Simulates GBM price paths using the Euler-Maruyama scheme.
  * **Likelihood Comparison:** Implements both exact and approximate marginal likelihood functions.
  * **Analytic Marginalization:** Treats the drift ($\mu$) as a nuisance parameter and integrates it out to improve sampling efficiency.
  * **Adaptive MCMC:** Uses a Metropolis-Hastings sampler with adaptive tuning of the proposal distribution during the burn-in phase.
  * **Diagnostics:** Includes calculation of Bias, RMSE, Coverage probability, Effective Sample Size (ESS), and visualizations (Trace plots, ACF, Posterior Density).

## 🛠️ Dependencies

The code is written in **R** and requires the following packages:

  * `ggplot2` (for visualization)
  * `coda` (for MCMC diagnostics)

You can install them via CRAN:

```r
install.packages(c("ggplot2", "coda"))
```

## 📂 Usage

1.  Clone this repository.
2.  Open the R script (e.g., `main.R`).
3.  Run the script. It is self-contained and performs the following steps:
      * Sets global parameters (Time horizon $T=10000$, $\sigma_{true}=0.5$, etc.).
      * Generates 100 synthetic GBM paths.
      * Runs the Adaptive MH sampler for both Exact and EM methods.
      * Outputs statistical summaries (Bias, RMSE, Coverage).
      * Generates diagnostic plots.

## 📊 Methodology

### Model

The asset price $S_t$ follows a Geometric Brownian Motion:
$$dS_t = \mu S_t dt + \sigma S_t dW_t$$

### Bayesian Inference

  * **Prior:** Inverse-Gamma prior for variance $\sigma^2 \sim IG(2, 0.25)$ and a flat prior for drift $\pi(\mu) \propto 1$.
  * **Sampling:** The algorithm targets the marginal posterior $p(\sigma^2 | Data)$ after integrating out $\mu$.
  * **Adaptation:** The proposal scale $\tau$ is updated using a Robbins-Monro scheme during burn-in to achieve a target acceptance rate of $\approx 0.4$.

## 📈 Results Summary

As detailed in the thesis, the simulation results demonstrate that:

  * Both methods provide well-calibrated 95% credible intervals.
  * For high-frequency data (large $T$), the **Euler-Maruyama approximation yields virtually identical posterior precision** to the exact method while being computationally efficient.
  * The bias vanishes as the time step $\Delta t$ decreases.

## 📝 Citation

If you use this code or the methodology in your research, please cite the thesis:

```bibtex
@mastersthesis{ZhangShuo2025,
  author  = {Shuo Zhang},
  title   = {Bayesian Statistics for Stochastic Processes in Finance},
  school  = {Faculty of Mathematics and Physics, Charles University},
  year    = {2025},
  address = {Prague},
  type    = {Master's Thesis}
}
```


## 📄 License

[ MIT License]

-----

