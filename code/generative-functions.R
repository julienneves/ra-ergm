GenerateData <- function(params_model, dgp_net){
  # Extract parameters
  alpha <- params_model$alpha
  beta <-  params_model$beta_X
  phi <-  params_model$phi
  psi <-  params_model$psi
  sigma_e <-  params_model$sigma_e
  
  G_true <- dgp_net$G_true
  
  n <-  G_true$gal$n
  xi <- dgp_net$xi
  epsilon <- rnorm(n, 0, sigma_e)
  
  # Create y and X
  X <- dgp_net$X
  
  y <- solve((diag(n) - phi * as.matrix(G_true)),  alpha + X %*% beta + psi * xi + epsilon)
  
  # Gather y, X in a data frame
  df <-  data.frame(y, X)
  
  return(df)
}

GenerateNetwork <- function(params_net){
  n <- params_net$n
  terms <- params_net$terms
  target_stats <- params_net$target_stats
  xi <- rnorm(n)
  
  # Create y and X
  X <- matrix(rnorm(n), n)
  
  formula_obs <- as.formula(paste("G ~ nodecov('xi') + ", paste(terms, collapse = " + ")))
  formula_true <- as.formula(paste("G ~ ", paste(terms, collapse = " + ")))
  
  G <- network.initialize(n , directed = FALSE)
  G %v% 'xi' <- xi
  
  # Generate ERGM distribution using seed
  net_formation_obs <- ergm(formula_obs, target.stats = c(mean(target_stats), target_stats))
  
  net_formation_true<- ergm(formula_true, target.stats = target_stats)
  net_formation_true$coef <- net_formation_obs$coef[terms]
  
  G_true <- simulate(net_formation_true)
  G_obs <- simulate(net_formation_obs)
  
  # Estimate ERGM on G
  net_formation_est <- ergm(as.formula(paste("G_obs ~ sociality +", paste(terms, collapse = " + "))))
  
  # Drop the sender/receiver (xi) from the coefficients of the model
  net_formation_est$coef <- net_formation_est$coef[terms]
  net_formation_est$formula <-  as.formula(paste("G_obs ~", paste(terms, collapse = " + ")))
  
  return(list(G_true = G_true, 
              G_obs = G_obs,
              xi = xi, 
              X = X,
              net_formation_est = net_formation_est, 
              net_formation_obs = net_formation_obs, 
              net_formation_true = net_formation_true))
}