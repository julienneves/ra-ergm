library(tidyr)
library(ggplot2)
library(parallel)

library(econet)
library(statnet)


GenerateNetwork <- function(params_net){
  n <- params_net$n
  terms <- params_net$terms
  target_stats <- params_net$target_stats
  xi <- rnorm(n)
  
  formula_obs <- as.formula(paste("G ~ nodecov('xi') + ", paste(terms, collapse = " + ")))
  formula_true <- as.formula(paste("G ~ ", paste(terms, collapse = " + ")))
  
  G <- network.initialize(n , directed = FALSE)
  G %v% 'xi' <- xi
  
  # Generate ERGM distribution using seed
  net_formation_obs <- ergm(formula_obs, target.stats = c(min(target_stats), target_stats))

  net_formation_true<- ergm(formula_true, target.stats = target_stats)
  net_formation_true$coef <- net_formation_obs$coef[terms]
  
  G_true <- simulate(net_formation_true)
  G_obs <- simulate(net_formation_obs)

  return(list(G_true = G_true, G_obs = G_obs, xi = xi, net_formation_obs = net_formation_obs, net_formation_true = net_formation_true))
}

EstimateModel <- function(dgp_net, params_net, params_model, repl = 100) {
  # Extract parameters
  alpha <- params_model$alpha
  beta <-  params_model$beta_X
  phi <-  params_model$phi
  psi <-  params_model$psi
  sigma_e <-  params_model$sigma_e
  
  coef_true<- list(alpha = alpha, 
                   beta_X = beta, 
                   phi = phi)
  
  n <-  params_net$n
  terms <- params_net$terms
  
  G_true <- dgp_net$G_true
  G_obs <- dgp_net$G_true
  
  xi <- dgp_net$xi
  
  # Create y and X
  X <- matrix(rnorm(n * length(beta)), n)
  epsilon <- rnorm(n, 0, sigma_e)
  y <- solve((diag(n) - phi * as.matrix(G_true)),  alpha + X %*% beta + psi * xi + epsilon)

  # Gather y, X in a data frame
  df <-  data.frame(y, X)

  # Estimate ERGM on G
  net_formation_est <- ergm(as.formula(paste("G_obs ~ sociality +", paste(terms, collapse = " + "))))

  # Drop the sendder/receiver (xi) from the coefficients of the model
  net_formation_est$coef <- net_formation_est$coef[terms]
  net_formation_est$formula <-  as.formula(paste("G_obs ~", paste(terms, collapse = " + ")))


  fit_nlls <- net_dep(formula = "y ~ X", df, G = as.matrix(G_obs), 
                      model = "model_B", estimation = "NLLS",  hypothesis = "lim", 
                      start.val = coef_true, endogeneity = FALSE)
  coef_nlls <- t(coef(fit_nlls[[1]]))
  coef_ergm <- data.frame(coef_nlls)

  for (i in 1:repl) {
    tryCatch({
      G_n <- simulate(net_formation_est)
      G_n <- as.matrix(G_n)
      fit <- net_dep(formula = "y ~ X", df, G = G_n, 
                     model = "model_B", estimation = "NLLS", hypothesis = "lim", 
                     start.val = coef_true, endogeneity = FALSE)
      coef_ergm[i, ] <- t(coef(fit[[1]]))
    })
  }
  return(list(coef_ergm = coef_ergm, coef_nlls = coef_nlls, coef_true = as.data.frame(coef_true)))
}


## Monte Carlo simulation
# Generate network
params_net = list(n = 50,
                  target_stats = c(10,100), 
                  terms = c("triangle", "edges"))

# Generate network
dgp_net <- GenerateNetwork(params_net)

# Plot true and observed network
par(mfrow=c(1,2)) 
plot(dgp_net$G_true, coord = expand.grid(1:ceiling(sqrt(params_net$n)),1:ceiling(sqrt(params_net$n)))[1:params_net$n,])
plot(dgp_net$G_obs,  coord = expand.grid(1:ceiling(sqrt(params_net$n)),1:ceiling(sqrt(params_net$n)))[1:params_net$n,])

# Run simulation
# Set parameters
params_model_1 <- list(alpha = .5, 
                     beta_X = .5, 
                     phi = 0,
                     psi = 1,
                     sigma_e = 1)
result_1 <- mclapply(1:100, function(x, ...) EstimateModel(...), dgp_net, params_net, params_model_1)

params_model_2 <- list(alpha = .5, 
                       beta_X = .5, 
                       phi = 0.8,
                       psi = 1,
                       sigma_e = 1)
result_2 <- mclapply(1:100, function(x, ...) EstimateModel(...), dgp_net, params_net, params_model_2)


coef_nlls <- t(sapply(result, function(x) x$coef_nlls))
coef_ergm <- t(sapply(result, function(x) mean(x$coef_ergm)))

# Plot coefficients
fig_1 <- ggplot2::ggplot(gather(coef_nlls), aes(value, fill=key)) +
  geom_density(kernel = "gaussian", alpha = 0.5) +
  geom_vline(aes(xintercept=value), data=gather(as.data.frame(t(result[[2]]))),linetype = 2)+
  facet_wrap( ~ key) +
  xlim(-1,1)

fig_1
