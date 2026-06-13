
### clear memory / fresh start
rm(list = ls())

### Table of Contents

# start -------------------------------------------------------------------
# Libraries ---------------------------------------------------------------

library(tidyverse)
library(dplyr)
library(tidyr)
library(cluster)
library(factoextra)
library(broom)
library(ggplot2)
library(ggcorrplot)
library(scales)

library(MASS)
library(quadprog)   # apply for portfolio opt.

select <- dplyr::select
filter <- dplyr::filter

if (!require("conflicted")) install.packages("conflicted")
library(conflicted)

# Force conflict resolution
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("lag", "dplyr")
conflict_prefer("map", "purrr")

# Load Data ---------------------------------------------------------------
setwd("/Users/chrisbyrialsen/Desktop/SDU/ds.econ/thesis ideas")

portfolio_dataset <- read.csv("data/Portfolio_Dataset.csv")
macro_regime_clean <- read.csv("data/Macro_Regime_Dataset.csv")

# PCA Loadings ---------------------------------------------------------------------

macro_pca_data <- macro_regime_clean %>%
  select(country, date, cpi, policy_rate, sovereign_yield) %>%
  na.omit()

macro_scaled <- macro_pca_data %>%
  select(cpi, policy_rate, sovereign_yield) %>%
  scale()

pca <- prcomp(macro_scaled, center = FALSE, scale. = FALSE)

summary(pca)
pca$rotation
plot(pca, type = "l")

macro_regime_pca <- macro_pca_data %>%
  mutate(
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    date = as.Date(paste0(date, "-01"))   # parse dates once, here
  )


# Global Monthly Regimes --------------------------------------------------

# Average PCs across countries per month → one row per date
macro_global <- macro_regime_pca %>%
  group_by(date) %>%
  summarise(
    PC1 = mean(PC1),
    PC2 = mean(PC2),
    .groups = "drop"
  )

# K-Selection: Elbow + Silhouette 
pca_data <- macro_global[, c("PC1", "PC2")]

# Elbow plot (WSS)
fviz_nbclust(
  pca_data,
  kmeans,
  method  = "wss",
  k.max   = 8,
  linecolor = "#2c3e50"
) +
  labs(
    title    = "Elbow Method: Optimal Number of Regimes",
    subtitle = "Look for the 'elbow' where WSS reduction flattens",
    x        = "Number of Clusters (k)",
    y        = "Total Within-Cluster Sum of Squares"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

ggsave("elbow_plot.pdf", width = 7, height = 5, device = "pdf")

# Silhouette plot
fviz_nbclust(
  pca_data,
  kmeans,
  method    = "silhouette",
  k.max     = 8,
  linecolor = "#e74c3c"
) +
  labs(
    title    = "Silhouette Method: Optimal Number of Regimes",
    subtitle = "Higher average silhouette = better-defined clusters",
    x        = "Number of Clusters (k)",
    y        = "Average Silhouette Width"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

ggsave("silhouette_plot.pdf", width = 7, height = 5, device = "pdf")

# BIC via GMM (formal model selection)
library(mclust)

gmm_bic <- mclustBIC(pca_data, G = 1:8)
plot(gmm_bic)
summary(gmm_bic)

set.seed(123)

kmeans_global <- kmeans(
  macro_global[, c("PC1", "PC2")],
  centers  = 3,
  iter.max = 100,
  nstart   = 50
)

macro_global$regime <- kmeans_global$cluster

# Inspect
table(macro_global$regime)

fviz_cluster(
  kmeans_global,
  data = macro_global[, c("PC1", "PC2")],
  ellipse.type = "norm"
)

# Label regimes based on average macro conditions
macro_global %>%
  left_join(
    macro_regime_pca %>%
      group_by(date) %>%
      summarise(across(c(cpi, policy_rate, sovereign_yield), mean)),
    by = "date"
  ) %>%
  group_by(regime) %>%
  summarise(across(c(cpi, policy_rate, sovereign_yield), mean), n = n())

# Assign names after inspecting the table above
macro_global <- macro_global %>%
  mutate(
    regime_name = case_when(
      regime == 1 ~ "Eased Monetary Conditions",
      regime == 2 ~ "Tight Monetary Conditions",
      regime == 3 ~ "Mixed Conditions"
    )
  )

save(macro_global, file = "Macro_Regime_Clusters.rdata")
write_csv(macro_global, "data/Macro_Regime_Clusters.csv")


# ETF Returns -------------------------------------------------------------

portfolio_returns <- portfolio_dataset %>%
  arrange(date) %>%
  mutate(
    date     = as.Date(paste0(date, "-01")),   # parse dates once, here
    IEF_ret  = IEF  / lag(IEF)  - 1,
    LQD_ret  = LQD  / lag(LQD)  - 1,
    HYG_ret  = HYG  / lag(HYG)  - 1,
    AIGI_ret = AIGI / lag(AIGI) - 1,
    AHYG_ret = AHYG / lag(AHYG) - 1
  ) %>%
  na.omit()


# Merge Returns with Regimes ----------------------------------------------

regime_returns <- portfolio_returns %>%
  left_join(
    macro_global %>% select(date, regime, regime_name),
    by = "date"
  ) %>%
  na.omit()   # drops dates with no regime match

nrow(regime_returns)
table(regime_returns$regime)


# Mean Returns by Regime --------------------------------------------------

regime_returns %>%
  group_by(regime, regime_name) %>%
  summarise(across(ends_with("_ret"), mean, na.rm = TRUE), .groups = "drop")


# Volatility by Regime ----------------------------------------------------

regime_returns %>%
  group_by(regime, regime_name) %>%
  summarise(across(ends_with("_ret"), sd, na.rm = TRUE), .groups = "drop")


# Sharpe Ratios by Regime -------------------------------------------------

regime_returns %>%
  group_by(regime, regime_name) %>%
  summarise(
    across(ends_with("_ret"),
           ~ mean(.x, na.rm = TRUE) / sd(.x, na.rm = TRUE)),
    .groups = "drop"
  )


# ANOVA Significance Tests ------------------------------------------------

etf_cols <- c("IEF_ret", "LQD_ret", "HYG_ret", "AIGI_ret", "AHYG_ret")

anova_results <- map(etf_cols, ~ {
  aov(reformulate("factor(regime)", response = .x), data = regime_returns)
})

names(anova_results) <- etf_cols

map(anova_results, summary)

# Post-hoc: Tukey HSD (run for any ETF with significant ANOVA)
map(anova_results, TukeyHSD)












# Phase 1: Credit Spreads -------------------------------------------------

# Compute implied spreads (excess return over IEF)
portfolio_returns <- portfolio_dataset %>%
  arrange(date) %>%
  mutate(
    date     = as.Date(paste0(date, "-01")),
    IEF_ret  = IEF  / lag(IEF)  - 1,
    LQD_ret  = LQD  / lag(LQD)  - 1,
    HYG_ret  = HYG  / lag(HYG)  - 1,
    AIGI_ret = AIGI / lag(AIGI) - 1,
    AHYG_ret = AHYG / lag(AHYG) - 1,
    BRJP_ret = BRJP / lag(BRJP) - 1
  ) %>%
  na.omit()

portfolio_returns <- portfolio_returns %>%
  mutate(
    spread_LQD  = LQD_ret  - IEF_ret,
    spread_HYG  = HYG_ret  - IEF_ret,
    spread_AIGI = AIGI_ret - IEF_ret,
    spread_AHYG = AHYG_ret - IEF_ret,
    spread_BRJP = BRJP_ret - IEF_ret
  )

# Attach global regime
portfolio_returns <- portfolio_returns %>%
  left_join(
    macro_global %>% select(date, regime, regime_name),
    by = "date"
  ) %>%
  na.omit()

# Inspect spreads by regime
portfolio_returns %>%
  group_by(regime, regime_name) %>%
  summarise(
    across(starts_with("spread_"), \(x) mean(x, na.rm = TRUE)),
    .groups = "drop"
  )

# Spread forecasting models (AR(1) + regime
spread_cols <- c("spread_LQD", "spread_HYG", "spread_AIGI", "spread_AHYG", "spread_BRJP")
spread_models <- map(spread_cols, ~ {
  df <- portfolio_returns %>%
    arrange(date) %>%
    mutate(
      spread_lag = lag(.data[[.x]]),
      regime_f   = factor(regime)
    ) %>%
    na.omit()
  
  lm(reformulate(c("spread_lag", "regime_f"), response = .x), data = df)
})

names(spread_models) <- spread_cols

# Results
map(spread_models, tidy)
map(spread_models, glance)

# Phase 2 & 3: Returns/MPT -------

# 1. Define the MVO function (Requires library(quadprog) active!)
mvo <- function(mu, sigma, gamma = 3) {
  n <- length(mu)
  Dmat <- gamma * sigma
  dvec <- mu
  Amat <- cbind(rep(1, n), diag(n))
  bvec <- c(1, rep(0, n))
  solve.QP(Dmat, dvec, Amat, bvec, meq = 1)$solution
}

# 2. Forecast Next Month's Expected Returns using CURRENT Month's Data
expected_returns <- portfolio_returns %>%
  arrange(date) %>%
  mutate(
    exp_LQD  = predict(spread_models$spread_LQD,  newdata = tibble(spread_lag = spread_LQD,  regime_f = factor(regime))),
    exp_HYG  = predict(spread_models$spread_HYG,  newdata = tibble(spread_lag = spread_HYG,  regime_f = factor(regime))),
    exp_AIGI = predict(spread_models$spread_AIGI, newdata = tibble(spread_lag = spread_AIGI, regime_f = factor(regime))),
    exp_AHYG = predict(spread_models$spread_AHYG, newdata = tibble(spread_lag = spread_AHYG, regime_f = factor(regime))),
    exp_BRJP = predict(spread_models$spread_BRJP, newdata = tibble(spread_lag = spread_BRJP, regime_f = factor(regime))),
    exp_IEF  = mean(IEF_ret, na.rm = TRUE) # Static baseline anchor
  ) %>%
  # Push the forecasts forward by 1 month so they apply to the trading month
  mutate(across(starts_with("exp_"), lag)) %>%
  na.omit() %>%
  dplyr::select(date, regime, regime_name, starts_with("exp_"))

# 3. Attach actual returns for the trading month
opt_data <- expected_returns %>%
  left_join(
    portfolio_returns %>% dplyr::select(date, IEF_ret, LQD_ret, HYG_ret, AIGI_ret, AHYG_ret, BRJP_ret),
    by = "date"
  )

etf_cols <- c("IEF_ret", "LQD_ret", "HYG_ret", "AIGI_ret", "AHYG_ret", "BRJP_ret")
exp_cols  <- c("exp_IEF", "exp_LQD", "exp_HYG", "exp_AIGI", "exp_AHYG", "exp_BRJP")

# 4. Calculate Historical Regime Covariance Matrices
regime_cov <- opt_data %>%
  group_by(regime) %>%
  summarise(cov_matrix = list(cov(across(all_of(etf_cols)))), .groups = "drop")

# 5. Generate Weights for the Trading Month 
weights <- opt_data %>%
  mutate(trading_regime = lag(regime)) %>% 
  filter(!is.na(trading_regime)) %>% 
  left_join(regime_cov, by = c("trading_regime" = "regime")) %>%
  rowwise() %>%
  mutate(
    mu     = list(c_across(all_of(exp_cols))),
    w      = list(mvo(unlist(mu), cov_matrix)),
    w_IEF  = w[[1]],
    w_LQD  = w[[2]],
    w_HYG  = w[[3]],
    w_AIGI = w[[4]],
    w_AHYG = w[[5]],
    w_BRJP = w[[6]]
  ) %>%
  ungroup() %>%
  dplyr::select(date, regime, regime_name, starts_with("w_"))

# Phase 4: Back-testing ----------------------------------------------------

# Attach actual returns to weights
backtest <- weights %>%
  left_join(
    portfolio_returns %>% dplyr::select(date, IEF_ret, LQD_ret, HYG_ret, AIGI_ret, AHYG_ret, BRJP_ret),
    by = "date"
  )

# Dynamic strategy return (regime-based MVO)
backtest <- backtest %>%
  mutate(
    ret_dynamic = w_IEF  * IEF_ret  +
      w_LQD  * LQD_ret  +
      w_HYG  * HYG_ret  +
      w_AIGI * AIGI_ret +
      w_AHYG * AHYG_ret +
      w_BRJP * BRJP_ret
  )

### TransactionCost and TurnoverConstraint

backtest <- backtest %>%
  arrange(date) %>%
  mutate(
    # 1. Calculate the turnover (absolute weight change) per asset
    # Using rowSums on the differences of columns starting with "w_"
    turnover = rowSums(abs(across(starts_with("w_"), ~ .x - coalesce(lag(.x), .x)))),
    
    # 2. Apply transaction cost penalty (assumed avg 0.2% or 20bps per trade)
    # This reflects slippage and brokerage fees
    tx_cost = turnover * 0.0020,
    
    # 3. Calculate Net Strategy Return
    ret_dynamic_net = ret_dynamic - tx_cost
  )

# Add a constraint check: 
# If turnover is > 0.5 (50% shift), it's a high-friction regime
backtest <- backtest %>%
  mutate(friction_regime = ifelse(turnover > 0.5, "High Friction", "Low Friction"))

# Benchmark 1: Equal-weight (1/6 each)
backtest <- backtest %>%
  mutate(
    ret_ew = (IEF_ret + LQD_ret + HYG_ret + AIGI_ret + AHYG_ret + BRJP_ret) / 6
  )

# Benchmark 2: Static MVO (unconditional moments, no regime)
mu_static    <- colMeans(portfolio_returns[, paste0(c("IEF","LQD","HYG","AIGI","AHYG","BRJP"), "_ret")])
sigma_static <- cov(portfolio_returns[, paste0(c("IEF","LQD","HYG","AIGI","AHYG","BRJP"), "_ret")])
w_static     <- mvo(mu_static, sigma_static)

backtest <- backtest %>%
  mutate(
    ret_static = w_static[1] * IEF_ret  +
      w_static[2] * LQD_ret  +
      w_static[3] * HYG_ret  +
      w_static[4] * AIGI_ret +
      w_static[5] * AHYG_ret +
      w_static[6] * BRJP_ret
  )

# Performance summary
backtest %>%
  summarise(
    across(
      c(ret_dynamic, ret_ew, ret_static),
      list(
        mean   = \(x) mean(x, na.rm = TRUE),
        sd     = \(x) sd(x, na.rm = TRUE),
        sharpe = \(x) mean(x, na.rm = TRUE) / sd(x, na.rm = TRUE)
      )
    )
  ) %>%
  pivot_longer(everything()) %>%
  separate(name, into = c("strategy", "metric"), sep = "_(?=[^_]+$)") %>%
  pivot_wider(names_from = metric, values_from = value)

# Phase 4.1: CVaR and Max Drawdown =============================

cvar_95 <- function(x, alpha = 0.05) {
  threshold <- quantile(x, alpha, na.rm = TRUE)
  mean(x[x <= threshold], na.rm = TRUE)
}

max_drawdown <- function(x) {
  cum <- cumprod(1 + x)
  running_max <- cummax(cum)
  drawdown <- (cum - running_max) / running_max
  min(drawdown)
}

# Extended performance summary
backtest %>%
  summarise(
    across(
      c(ret_dynamic, ret_ew, ret_static),
      list(
        mean    = \(x) mean(x, na.rm = TRUE),
        sd      = \(x) sd(x, na.rm = TRUE),
        sharpe  = \(x) mean(x, na.rm = TRUE) / sd(x, na.rm = TRUE),
        cvar95  = \(x) cvar_95(x),
        max_dd  = \(x) max_drawdown(x)
      )
    )
  ) %>%
  pivot_longer(everything()) %>%
  separate(name, into = c("strategy", "metric"), sep = "_(?=[^_]+$)") %>%
  pivot_wider(names_from = metric, values_from = value)


# Phase 4.2: Portfolio Diagnostics =================================================

ret_cols    <- c("IEF_ret","LQD_ret","HYG_ret","AIGI_ret","AHYG_ret","BRJP_ret")
asset_names <- c("IEF","LQD","HYG","AIGI","AHYG","BRJP")
regime_list <- unique(portfolio_returns$regime_name)

# --- 1. Regime-conditional correlation matrices ---------------------------

cor_by_regime <- map_dfr(regime_list, function(reg) {
  df <- portfolio_returns %>% filter(regime_name == reg)
  cor_mat <- cor(df[, ret_cols], use = "complete.obs")
  
  as.data.frame(cor_mat) %>%
    rownames_to_column("asset1") %>%
    pivot_longer(-asset1, names_to = "asset2", values_to = "correlation") %>%
    mutate(
      asset1 = recode(asset1, !!!setNames(asset_names, ret_cols)),
      asset2 = recode(asset2, !!!setNames(asset_names, ret_cols)),
      regime_name = reg
    )
})

p_corr <- ggplot(cor_by_regime, aes(x = asset1, y = asset2, fill = correlation)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", correlation)), size = 3) +
  scale_fill_gradient2(
    low = "#2980b9", mid = "white", high = "#e74c3c",
    midpoint = 0, limits = c(-1, 1), name = "Correlation"
  ) +
  facet_wrap(~regime_name, ncol = 1) +
  labs(
    title    = "Regime-Conditional Correlation Matrices",
    subtitle = "Pairwise return correlations by macroeconomic regime",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title  = element_text(face = "bold"),
    strip.text  = element_text(face = "bold")
  )

p_corr

ggsave("correlation_matrices.pdf", plot = p_corr, width = 7, height = 14, device = "pdf")


# --- 2. Risk Contribution table --------------------------------------------

risk_contribution <- function(w, sigma) {
  w <- as.numeric(w)
  port_var <- as.numeric(t(w) %*% sigma %*% w)
  marginal <- as.numeric(sigma %*% w)
  rc <- w * marginal
  rc / port_var   # normalized to sum to 1
}

avg_weights <- weights %>%
  group_by(regime_name) %>%
  summarise(
    across(starts_with("w_"), mean),
    .groups = "drop"
  )

risk_contrib_table <- map_dfr(regime_list, function(reg) {
  w_row <- avg_weights %>% filter(regime_name == reg)
  w_vec <- w_row %>% select(starts_with("w_")) %>% as.numeric()
  
  regime_id <- portfolio_returns %>%
    filter(regime_name == reg) %>%
    pull(regime) %>%
    first()
  
  sigma <- regime_cov %>%
    filter(regime == regime_id) %>%
    pull(cov_matrix) %>%
    .[[1]]
  
  rc <- risk_contribution(w_vec, sigma)
  
  tibble(
    regime_name = reg,
    asset       = asset_names,
    weight      = w_vec,
    risk_contribution = rc
  )
})

risk_contrib_table


# --- 3. Risk Premium vs Risk Contribution table -----------------------------

avg_premiums <- expected_returns %>%
  group_by(regime_name) %>%
  summarise(
    across(starts_with("exp_"), mean),
    .groups = "drop"
  ) %>%
  pivot_longer(-regime_name, names_to = "asset", values_to = "risk_premium") %>%
  mutate(asset = gsub("exp_", "", asset))

premium_vs_risk <- risk_contrib_table %>%
  left_join(avg_premiums, by = c("regime_name", "asset")) %>%
  mutate(premium_per_risk = risk_premium / risk_contribution)

premium_vs_risk
# Phase 5: Jensen's Alpha -------------------------------------------------

# Use IEF as the market factor (bond market benchmark)
backtest <- backtest %>%
  mutate(
    excess_dynamic = ret_dynamic - IEF_ret,
    excess_static  = ret_static  - IEF_ret,
    excess_ew      = ret_ew      - IEF_ret
  )

alpha_dynamic <- lm(excess_dynamic ~ IEF_ret, data = backtest)
alpha_static  <- lm(excess_static  ~ IEF_ret, data = backtest)
alpha_ew      <- lm(excess_ew      ~ IEF_ret, data = backtest)

map(
  list(dynamic = alpha_dynamic, static = alpha_static, ew = alpha_ew),
  tidy
)

# Phase 6: Monte Carlo Robustness -----------------------------------------

set.seed(123)
n_sim    <- 1000
n_months <- nrow(backtest)

# --- 6a: Bootstrap (primary) 

boot_sharpe <- replicate(n_sim, {
  idx <- sample(1:n_months, n_months, replace = TRUE)
  
  r_dynamic <- backtest$ret_dynamic[idx]
  r_static  <- backtest$ret_static[idx]
  r_ew      <- backtest$ret_ew[idx]
  
  c(
    sharpe_dynamic = mean(r_dynamic) / sd(r_dynamic),
    sharpe_static  = mean(r_static)  / sd(r_static),
    sharpe_ew      = mean(r_ew)      / sd(r_ew)
  )
}) %>%
  t() %>%
  as.data.frame()

cat("--- Bootstrap Results ---\n")
summary(boot_sharpe)
cat("P(dynamic > static): ",  mean(boot_sharpe$sharpe_dynamic > boot_sharpe$sharpe_static), "\n")
cat("P(dynamic > ew):     ",  mean(boot_sharpe$sharpe_dynamic > boot_sharpe$sharpe_ew), "\n")










# Forecast ---------------------------------------------------------------------

ret_cols <- c(
  "IEF_ret",
  "LQD_ret",
  "HYG_ret",
  "AIGI_ret",
  "AHYG_ret",
  "BRJP_ret"
)

forecast_models <- lapply(ret_cols, function(asset){
  
  formula <- as.formula(
    paste0(
      asset,
      " ~ factor(regime_name)"
    )
  )
  
  lm(formula,
     data = portfolio_returns)
  
})

names(forecast_models) <- ret_cols


# Forecast Expected Returns (today) -----------------------------------------------

current_regime <- "Eased Monetary Conditions"

expected_returns <- sapply(
  seq_along(ret_cols),
  function(i){
    
    predict(
      forecast_models[[i]],
      newdata = data.frame(
        regime_name = current_regime
      )
    )
    
  }
)

names(expected_returns) <- ret_cols

expected_returns


# Estimate Covariance Matrix ----------------------------------------------

regime_data <- portfolio_returns %>%
  filter(
    regime_name == current_regime
  )

Sigma <- cov(
  regime_data[, ret_cols],
  use = "complete.obs"
)

Sigma


# Portfolio Optimization --------------------------------------------------

n <- length(expected_returns)

Amat <- cbind(
  rep(1,n),
  diag(n)
)

bvec <- c(
  1,
  rep(0,n)
)

sol <- solve.QP(
  Dmat = 2*Sigma,
  dvec = expected_returns,
  Amat = Amat,
  bvec = bvec,
  meq = 1
)

weights <- sol$solution

weights


# Time-Series visualization ------------------------------------------------------

# Regime shading helper
regime_bands <- backtest %>%
  arrange(date) %>%
  mutate(
    regime_change = regime != lag(regime, default = first(regime)),
    band_id       = cumsum(regime_change)
  ) %>%
  group_by(band_id, regime, regime_name) %>%
  summarise(
    xmin = min(date),
    xmax = max(date),
    .groups = "drop"
  )

regime_colors <- c(
  "Eased Monetary Conditions" = "#d4e6f1",
  "Tight Monetary Conditions" = "#fadbd8",
  "Mixed Conditions"          = "#d5f5e3"
)

# Cumulative returns
cumret <- backtest %>%
  arrange(date) %>%
  mutate(
    cum_dynamic = cumprod(1 + ret_dynamic),
    cum_static  = cumprod(1 + ret_static),
    cum_ew      = cumprod(1 + ret_ew)
  ) %>%
  select(date, cum_dynamic, cum_static, cum_ew) %>%
  pivot_longer(-date, names_to = "strategy", values_to = "value") %>%
  mutate(strategy = recode(strategy,
                           cum_dynamic = "Dynamic (Regime MVO)",
                           cum_static  = "Static MVO",
                           cum_ew      = "Equal-Weight"
  ))

p_cumret <- ggplot() +
  geom_rect(
    data = regime_bands,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = regime_name),
    alpha = 0.3
  ) +
  scale_fill_manual(values = regime_colors, name = "Regime") +
  geom_line(
    data = cumret,
    aes(x = date, y = value, color = strategy, linetype = strategy),
    linewidth = 0.8
  ) +
  scale_color_manual(
    values = c(
      "Dynamic (Regime MVO)" = "#2c3e50",
      "Static MVO"           = "#e74c3c",
      "Equal-Weight"         = "#3498db"
    ),
    name = "Strategy"
  ) +
  scale_linetype_manual(
    values = c(
      "Dynamic (Regime MVO)" = "solid",
      "Static MVO"           = "dashed",
      "Equal-Weight"         = "dotted"
    ),
    name = "Strategy"
  ) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  labs(
    title    = "Cumulative Portfolio Performance by Strategy",
    subtitle = "Shaded regions indicate macroeconomic regime",
    x        = NULL,
    y        = "Cumulative Return (USD, base = 1)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position   = "bottom",
    panel.grid.minor  = element_blank(),
    plot.title        = element_text(face = "bold"),
    legend.box        = "vertical"
  )

p_cumret

ggsave(
  "cumulative_returns.pdf",
  plot   = p_cumret,
  width  = 10,
  height = 6,
  device = "pdf"
)





# Credit Spread Visualzation ----------------------------------------------

spread_long <- portfolio_returns %>%
  select(date, regime_name, spread_LQD, spread_HYG, spread_AIGI, spread_AHYG, spread_BRJP) %>%
  pivot_longer(-c(date, regime_name), names_to = "asset", values_to = "spread") %>%
  mutate(asset = recode(asset,
                        spread_LQD  = "LQD (US IG)",
                        spread_HYG  = "HYG (US HY)",
                        spread_AIGI = "AIGI (ASEAN IG)",
                        spread_AHYG = "AHYG (ASEAN HY)",
                        spread_BRJP = "BRJP (ASEAN Broad)"
  ))

p_spreads <- ggplot() +
  geom_rect(
    data = regime_bands,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = regime_name),
    alpha = 0.3
  ) +
  scale_fill_manual(values = regime_colors, name = "Regime") +
  geom_line(
    data = spread_long,
    aes(x = date, y = spread, color = asset),
    linewidth = 0.6
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  facet_wrap(~asset, ncol = 1, scales = "free_y") +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  scale_color_brewer(palette = "Set1", guide = "none") +
  labs(
    title    = "Implied Credit Spreads by Asset and Macroeconomic Regime",
    subtitle = "Spread = ETF return minus IEF return; shaded regions indicate regime",
    x        = NULL,
    y        = "Monthly Implied Spread"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position  = "bottom",
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold"),
    strip.text       = element_text(face = "bold")
  )

p_spreads

ggsave(
  "credit_spreads.pdf",
  plot   = p_spreads,
  width  = 10,
  height = 12,
  device = "pdf"
)

# Efficient Frontier Visuals ------------------------------------------------------

ret_cols <- c("IEF_ret","LQD_ret","HYG_ret","AIGI_ret","AHYG_ret","BRJP_ret")
asset_names <- c("IEF","LQD","HYG","AIGI","AHYG","BRJP")

regime_list <- unique(portfolio_returns$regime_name)

frontier_by_regime <- map_dfr(regime_list, function(reg) {
  df <- portfolio_returns %>% filter(regime_name == reg)
  
  mu_r    <- colMeans(df[, ret_cols])
  sigma_r <- cov(df[, ret_cols])
  n       <- length(mu_r)
  
  targets <- seq(min(mu_r), max(mu_r) * 1.1, length.out = 200)
  
  map_dfr(targets, function(target) {
    tryCatch({
      Amat <- cbind(rep(1, n), mu_r, diag(n))
      bvec <- c(1, target, rep(0, n))
      sol  <- solve.QP(2 * sigma_r, rep(0, n), Amat, bvec, meq = 2)
      tibble(ret = target, vol = sqrt(sol$value), regime_name = reg)
    }, error = function(e) NULL)
  })
})

assets_by_regime <- map_dfr(regime_list, function(reg) {
  df <- portfolio_returns %>% filter(regime_name == reg)
  mu_r    <- colMeans(df[, ret_cols])
  sigma_r <- cov(df[, ret_cols])
  tibble(
    asset       = asset_names,
    ret         = mu_r,
    vol         = sqrt(diag(sigma_r)),
    regime_name = reg
  )
})

strategies_by_regime <- backtest %>%
  group_by(regime_name) %>%
  summarise(
    dynamic_ret = mean(ret_dynamic),
    dynamic_vol = sd(ret_dynamic),
    static_ret  = mean(ret_static),
    static_vol  = sd(ret_static),
    ew_ret      = mean(ret_ew),
    ew_vol      = sd(ret_ew),
    .groups = "drop"
  )

p_frontier_regime <- ggplot() +
  geom_path(
    data = frontier_by_regime,
    aes(x = vol, y = ret),
    color = "#2c3e50", linewidth = 1.2
  ) +
  geom_point(
    data = assets_by_regime,
    aes(x = vol, y = ret),
    shape = 21, size = 3, fill = "#bdc3c7", color = "grey30"
  ) +
  geom_text(
    data = assets_by_regime,
    aes(x = vol, y = ret, label = asset),
    vjust = -1, size = 3, color = "grey30"
  ) +
  geom_point(
    data = strategies_by_regime,
    aes(x = dynamic_vol, y = dynamic_ret),
    shape = 17, size = 4, color = "#e74c3c"
  ) +
  geom_text(
    data = strategies_by_regime,
    aes(x = dynamic_vol, y = dynamic_ret, label = "Dynamic"),
    vjust = 2, size = 3, color = "#e74c3c"
  ) +
  geom_point(
    data = strategies_by_regime,
    aes(x = static_vol, y = static_ret),
    shape = 15, size = 4, color = "#27ae60"
  ) +
  geom_text(
    data = strategies_by_regime,
    aes(x = static_vol, y = static_ret, label = "Static MVO"),
    vjust = 2, size = 3, color = "#27ae60"
  ) +
  geom_point(
    data = strategies_by_regime,
    aes(x = ew_vol, y = ew_ret),
    shape = 16, size = 4, color = "#2980b9"
  ) +
  geom_text(
    data = strategies_by_regime,
    aes(x = ew_vol, y = ew_ret, label = "Equal-Weight"),
    vjust = 2, size = 3, color = "#2980b9"
  ) +
  facet_wrap(~regime_name, ncol = 1, scales = "free") +
  labs(
    title    = "Regime-Conditional Efficient Frontiers",
    subtitle = "Frontier and strategy positions vary across macroeconomic regimes",
    x        = "Monthly Volatility (Std. Dev.)",
    y        = "Monthly Mean Return"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold"),
    strip.text       = element_text(face = "bold", size = 11)
  )

p_frontier_regime

ggsave(
  "efficient_frontier_by_regime.pdf",
  plot   = p_frontier_regime,
  width  = 8,
  height = 14,
  device = "pdf"
)












# SML / Jensen's Alpha visual -----------------------------------

# Extract alpha and beta for each strategy
sml_data <- map_dfr(
  list(
    "Dynamic (Regime MVO)" = alpha_dynamic,
    "Static MVO"           = alpha_static,
    "Equal-Weight"         = alpha_ew
  ),
  ~ tibble(
    alpha = coef(.x)[1],
    beta  = coef(.x)[2]
  ),
  .id = "strategy"
)

# Market (IEF) reference point
market_point <- tibble(
  strategy = "Market (IEF)",
  beta     = 1,
  ret      = mean(backtest$IEF_ret)
)

# SML line
sml_line <- tibble(
  beta = seq(-0.5, 1.5, length.out = 100),
  ret  = mean(backtest$IEF_ret) + beta * mean(backtest$IEF_ret)
)

# Strategy expected returns (alpha + beta * market)
sml_data <- sml_data %>%
  mutate(ret = alpha + beta * mean(backtest$IEF_ret))

p_sml <- ggplot() +
  geom_line(
    data = sml_line,
    aes(x = beta, y = ret),
    color = "grey40", linewidth = 1, linetype = "dashed"
  ) +
  geom_point(
    data = sml_data,
    aes(x = beta, y = ret, color = strategy, shape = strategy),
    size = 5
  ) +
  geom_text(
    data = sml_data,
    aes(x = beta, y = ret, label = strategy, color = strategy),
    vjust = -1, size = 3.5
  ) +
  geom_point(
    data = market_point,
    aes(x = beta, y = ret),
    shape = 21, size = 4, fill = "white", color = "grey40"
  ) +
  geom_text(
    data = market_point,
    aes(x = beta, y = ret, label = strategy),
    vjust = -1, size = 3.5, color = "grey40"
  ) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey60") +
  scale_color_manual(
    values = c(
      "Dynamic (Regime MVO)" = "#e74c3c",
      "Static MVO"           = "#27ae60",
      "Equal-Weight"         = "#2980b9"
    ),
    name = NULL
  ) +
  scale_shape_manual(
    values = c(
      "Dynamic (Regime MVO)" = 17,
      "Static MVO"           = 15,
      "Equal-Weight"         = 16
    ),
    name = NULL
  ) +
  labs(
    title    = "Security Market Line and Jensen's Alpha",
    subtitle = "Strategies plotted relative to bond market factor (IEF)",
    x        = "Beta (IEF)",
    y        = "Monthly Expected Return"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position  = "none",
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold")
  )

p_sml

ggsave(
  "security_market_line.pdf",
  plot   = p_sml,
  width  = 8,
  height = 6,
  device = "pdf"
)

# Jensen's Alpha Bar Chart ------------------------------------------------

alpha_bars <- sml_data %>%
  select(strategy, alpha) %>%
  mutate(positive = alpha > 0)

p_alpha <- ggplot(alpha_bars, aes(x = strategy, y = alpha, fill = positive)) +
  geom_col(width = 0.5) +
  geom_hline(yintercept = 0, color = "grey30", linewidth = 0.8) +
  geom_text(
    aes(
      label = sprintf("%.4f", alpha),
      vjust = ifelse(positive, -0.5, 1.5)
    ),
    size = 4, color = "grey20"
  ) +
  scale_fill_manual(
    values = c("TRUE" = "#1abc9c", "FALSE" = "#e74c3c"),
    labels = c("TRUE" = "Alpha > 0", "FALSE" = "Alpha < 0"),
    name   = NULL
  ) +
  scale_x_discrete(limits = c("Equal-Weight", "Static MVO", "Dynamic (Regime MVO)")) +
  labs(
    title    = "Jensen's Alpha by Strategy",
    subtitle = "Alpha estimated relative to bond market factor (IEF)",
    x        = NULL,
    y        = "Monthly Alpha"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position  = "bottom",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title       = element_text(face = "bold")
  )

p_alpha

ggsave(
  "jensens_alpha.pdf",
  plot   = p_alpha,
  width  = 7,
  height = 5,
  device = "pdf"
)

# Bootstrap Sharpe Distributions ------------------------------------------

boot_long <- boot_sharpe %>%
  pivot_longer(everything(), names_to = "strategy", values_to = "sharpe") %>%
  mutate(strategy = recode(strategy,
                           sharpe_dynamic = "Dynamic (Regime MVO)",
                           sharpe_static  = "Static MVO",
                           sharpe_ew      = "Equal-Weight"
  ))

p_boot <- ggplot(boot_long, aes(x = sharpe, fill = strategy, color = strategy)) +
  geom_density(alpha = 0.3, linewidth = 0.8) +
  geom_vline(
    data = boot_long %>%
      group_by(strategy) %>%
      summarise(median = median(sharpe), .groups = "drop"),
    aes(xintercept = median, color = strategy),
    linetype = "dashed", linewidth = 0.8
  ) +
  scale_fill_manual(
    values = c(
      "Dynamic (Regime MVO)" = "#e74c3c",
      "Static MVO"           = "#27ae60",
      "Equal-Weight"         = "#2980b9"
    ),
    name = NULL
  ) +
  scale_color_manual(
    values = c(
      "Dynamic (Regime MVO)" = "#e74c3c",
      "Static MVO"           = "#27ae60",
      "Equal-Weight"         = "#2980b9"
    ),
    name = NULL
  ) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey40") +
  labs(
    title    = "Bootstrap Distribution of Sharpe Ratios",
    subtitle = "1,000 bootstrap resamples; dashed lines indicate median",
    x        = "Sharpe Ratio",
    y        = "Density"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position  = "bottom",
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold")
  )

p_boot

ggsave(
  "bootstrap_sharpe.pdf",
  plot   = p_boot,
  width  = 8,
  height = 5,
  device = "pdf"
)




# prelim: for forecast Vis... ----------------------------------------------

# 1. Rebuild weights (got overwritten earlier)
weights <- opt_data %>%
  mutate(trading_regime = lag(regime)) %>% 
  filter(!is.na(trading_regime)) %>% 
  left_join(regime_cov, by = c("trading_regime" = "regime")) %>%
  rowwise() %>%
  mutate(
    mu     = list(c_across(all_of(exp_cols))),
    w      = list(mvo(unlist(mu), cov_matrix)),
    w_IEF  = w[[1]],
    w_LQD  = w[[2]],
    w_HYG  = w[[3]],
    w_AIGI = w[[4]],
    w_AHYG = w[[5]],
    w_BRJP = w[[6]]
  ) %>%
  ungroup() %>%
  dplyr::select(date, regime, regime_name, starts_with("w_"))

# 2. Build implementation table
implementation_table <- weights %>%
  group_by(regime_name) %>%
  summarise(
    Avg_IEF  = mean(w_IEF),
    Avg_LQD  = mean(w_LQD),
    Avg_HYG  = mean(w_HYG),
    Avg_AIGI = mean(w_AIGI),
    Avg_AHYG = mean(w_AHYG),
    Avg_BRJP = mean(w_BRJP),
    .groups = "drop"
  )

# Forecast Visualization --------------------------------------------------

# 1. Take the summary table you already generated successfully
as_tibble(implementation_table) %>%
  # 2. Reshape all the Avg_ asset columns into long format
  pivot_longer(
    cols = starts_with("Avg_"), 
    names_to = "Asset", 
    values_to = "Weight"
  ) %>%
  # 3. Clean up the "Avg_" text prefix so the legend looks pristine
  mutate(Asset = gsub("Avg_", "", Asset)) %>%
  # 4. Build the visualization
  ggplot(aes(x = regime_name, y = Weight, fill = Asset)) +
  geom_bar(stat = "identity", position = "stack", width = 0.55) +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Dynamic Credit Allocation Across Macro Regimes",
    subtitle = "Empirical asset mix optimized via predictive walk-forward MVO",
    x = "Identified Macro Environment", 
    y = "Portfolio Allocation %", 
    fill = "Asset Class"
  ) +
  scale_y_continuous(labels = scales::percent)



# Friction Gap (cost basis) -----------------------------------------------

# 1. Prepare the data including the new net returns
cumret_gap <- backtest %>%
  arrange(date) %>%
  mutate(
    cum_gross = cumprod(1 + ret_dynamic),
    cum_net   = cumprod(1 + ret_dynamic_net), # Your new cost-adjusted data
    cum_ew    = cumprod(1 + ret_ew)
  ) %>%
  select(date, cum_gross, cum_net, cum_ew) %>%
  pivot_longer(-date, names_to = "strategy", values_to = "value") %>%
  mutate(strategy = recode(strategy,
                           cum_gross = "Dynamic (Gross)",
                           cum_net   = "Dynamic (Net of Fees)",
                           cum_ew    = "Equal-Weight"))

# 2. Build the chart
ggplot() +
  geom_line(data = cumret_gap, 
            aes(x = date, y = value, color = strategy, linetype = strategy), 
            linewidth = 1) +
  # Use shading to fill the "Friction Gap"
  geom_ribbon(data = cumret_gap %>% pivot_wider(names_from = strategy, values_from = value),
              aes(x = date, ymin = `Dynamic (Net of Fees)`, ymax = `Dynamic (Gross)`),
              fill = "red", alpha = 0.1) +
  scale_color_manual(values = c("Dynamic (Gross)" = "gray50", 
                                "Dynamic (Net of Fees)" = "#2c3e50", 
                                "Equal-Weight" = "#3498db")) +
  scale_linetype_manual(values = c("Dynamic (Gross)" = "dashed", 
                                   "Dynamic (Net of Fees)" = "solid", 
                                   "Equal-Weight" = "dotted")) +
  theme_minimal() +
  labs(title = "The Friction Gap: Gross vs. Net Dynamic Strategy Performance",
       subtitle = "Red area represents cumulative transaction cost impact (20bps per trade)",
       y = "Cumulative Wealth (Base = 1)", x = NULL)

# Risk Contribution Bar Chart ----------------------------------------------
p_riskcontrib <- risk_contrib_table %>%
  mutate(asset = factor(asset, levels = asset_names)) %>%
  ggplot(aes(x = regime_name, y = risk_contribution, fill = asset)) +
  geom_bar(stat = "identity", position = "stack", width = 0.55) +
  scale_fill_brewer(palette = "Set2") +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title    = "Risk Contribution by Asset and Regime",
    subtitle = "Share of total portfolio volatility attributable to each asset",
    x        = "Macroeconomic Regime",
    y        = "Risk Contribution (%)",
    fill     = "Asset"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position  = "bottom",
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold")
  )

p_riskcontrib

ggsave(
  "risk_contribution.pdf",
  plot   = p_riskcontrib,
  width  = 8,
  height = 6,
  device = "pdf"
)
# Risk Premium v. Risk Cont. Scatter ---------------------------------
p_premium_risk <- premium_vs_risk %>%
  mutate(asset = factor(asset, levels = asset_names)) %>%
  ggplot(aes(x = risk_contribution, y = risk_premium, color = asset)) +
  geom_point(size = 4) +
  geom_text(aes(label = asset), vjust = -1, size = 3.5, show.legend = FALSE) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey50") +
  facet_wrap(~regime_name, ncol = 1, scales = "free") +
  scale_color_brewer(palette = "Set2") +
  scale_x_continuous(labels = scales::percent) +
  labs(
    title    = "Risk Premium vs. Risk Contribution by Regime",
    subtitle = "Assets above/left are efficient: high premium relative to risk contributed",
    x        = "Risk Contribution (% of portfolio volatility)",
    y        = "Risk Premium (expected excess return)",
    color    = "Asset"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position  = "none",
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold"),
    strip.text       = element_text(face = "bold")
  )

p_premium_risk

ggsave(
  "premium_vs_risk_contribution.pdf",
  plot   = p_premium_risk,
  width  = 8,
  height = 14,
  device = "pdf"
)
# end ---------------------------------------------------------------------


