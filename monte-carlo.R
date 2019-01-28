library(econet)
library(statnet)
library(tidyr)
library(igraph)
library(ggplot2)

set.seed(1234)

net_mc <- function(params_true, psi, sigma_e, n, terms, repl = 100) {

  # Extract parameters
  alpha = params_true["alpha"]
  beta = params_true[c("beta_X")]
  phi = params_true["phi"]

  # Set NLLS formula
  formula <- "y ~ X"

  # Create a seed for network
  G <- erdos.renyi.game(n, 3*n, type="gnm", directed = FALSE)
  G <- as.network(get.adjacency(G))

  # Generate ERGM distribution using seed
  net_formation_obs <- ergm(as.formula(paste("G ~ sender + receiver +", terms)))

  net_formation_true <- net_formation_obs
  net_formation_true$coef <- net_formation_true$coef[terms]
  net_formation_true$formula <-  as.formula(paste("G ~", terms))

  # Simulate a network using the true model distribution
  G_true <- simulate(net_formation_true)
  G_obs <- simulate(net_formation_obs)
  plot(G)

  # Compute xi from ERGM distribution
  xi <- net_formation_obs$coef[1:(2 * (n - 1))]
  xi <- (xi[1:n - 1] + xi[n:(2 * (n - 1))]) / 2
  xi <- c(-sum(xi), xi)

  # Create y and X
  X <- matrix(rnorm(n * length(beta)), n)
  epsilon <- rnorm(n, 0, sigma_e)
  y <- solve((diag(n) - phi * as.matrix(G_true)),  alpha + X %*% beta + psi * xi + epsilon)

  # Gather y, X in a data frame
  df <-  data.frame(y, X)

  # Estimate ERGM on G
  net_formation_est <- ergm(as.formula(paste("G_obs ~ sender + receiver +", terms)))

  # Drop the xi from the coefficient of the model
  net_formation_est$coef <- net_formation_est$coef[terms]
  net_formation_est$formula <-  as.formula(paste("G ~", terms))


  fit_nlls <- net_dep(formula = formula, df, G = as.matrix(G_obs), model = "model_B", estimation = "NLLS",  hypothesis = "lim",start.val = params_true, endogeneity = FALSE)
  params_nlls <- coef(fit_nlls[[1]])
  params_ergm <- data.frame(t(params_true))

  for (i in 1:repl) {
    G_n <- simulate(net_formation_est)
    G_n <- as.matrix(G_n)
    tryCatch({
      fit <- net_dep(formula = formula, df, G = G_n, model = "model_B", estimation = "NLLS", hypothesis = "lim", start.val = params_true, endogeneity = FALSE)
      params_ergm[i, ] <- coef(fit[[1]])
    })
  }
  return(list(params_ergm, params_nlls, G = G))
}

## Example 1
# Monte Carlo simulation
params_true <- c(alpha = -1, beta_X = 1, phi = .9)
result <- net_mc(params_true, psi = 1, sigma_e = 10, n = 50, terms = c("triangle"), repl = 250)

# Plot coefficients
fig_1 <- ggplot(gather(result[[1]]), aes(value, fill=key)) +
  geom_density(kernel = "gaussian", alpha = 0.5) +
  geom_vline(aes(xintercept=params_true), data=gather(as.data.frame(t(params_true))))+
  geom_vline(aes(xintercept=value), data=gather(as.data.frame(t(result[[2]]))),linetype = 2)+
  facet_wrap( ~ key) +
  xlim(-5,5)

fig_1
