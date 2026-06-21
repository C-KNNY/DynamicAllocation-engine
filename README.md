# DynamicAllocation-engine
(RStudio: a regime-aware multi-asset credit allocation framework)

# Volatility

Realized volatility:

$$
\sigma_{real}
=
\sqrt{
252\cdot Var
\left(
\ln\frac{S_t}{S_{t-1}}
\right)
}
$$

---

# Portfolio Optimization

$$
w^*
=
\arg\max_w
\left(
\mu^Tw-\lambda w^T\Sigma w
\right)
$$
