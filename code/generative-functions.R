GenerateData <- function(params_model, dgp_net){
  # Extract parameters
  alpha <- params_model$alpha
  beta <-  params_model$beta_X
  phi <-  params_model$phi
  psi <- params_model$psi
  sigma_e <-  params_model$sigma_e
  
  G_true <- dgp_net$G_true
  X <- dgp_net$X
  xi <- dgp_net$xi
  
  n <-  G_true$gal$n
  epsilon <- rnorm(n, 0, sigma_e)
  
  # Create y
  y <- solve((diag(n) - phi * as.matrix(G_true)),  alpha + X %*% beta + psi * xi + epsilon)
  
  # Gather y, X in a data frame
  df <-  data.frame(y, X)
  
  return(df)
}

GenerateNetwork <- function(params_net){
  n <- params_net$n
  terms_true <- params_net$terms_true
  terms_alumni <- params_net$terms_alumni
  target_stats_alumni <- params_net$target_stats_alumni
  target_stats_true <- params_net$target_stats_true
  target_stats_xi <- params_net$target_stats_xi
  xi <- params_net$xi

  # Create y and X
  X <- matrix(rnorm(n), n)
  
  formula_alumni <- as.formula(paste("G ~ ", paste(terms_alumni, collapse = " + ")))
  formula_obs <- as.formula(paste("G ~ nodecov('xi') + ", paste(terms_true, collapse = " + ")))
  formula_true <- as.formula(paste("G ~ ", paste(terms_true, collapse = " + ")))
  formula_est <- as.formula(paste("G_obs ~ sociality + ", paste(terms_true, collapse = " + ")))

  G <- network.initialize(n , directed = FALSE)

  # Generate ERGM distribution using seed
  
  if (target_stats_alumni == "density"){
    G_alumni<-network(n, directed=FALSE, density=0.1)
  } else{
    net_formation_alumni <- ergm(formula_alumni, target.stats = target_stats_alumni)
    G_alumni <- simulate(net_formation_alumni)   
  }
  
  G %v% 'xi' <- xi
  
  net_formation_obs <- ergm(formula_obs, target.stats = c(target_stats_xi, target_stats_true))
  net_formation_true <- ergm(formula_true, target.stats = target_stats_true)
  net_formation_true$coef <- net_formation_obs$coef[-1]
  net_formation_true$MCMCtheta <- net_formation_obs$MCMCtheta[-1]
  
  G_obs <- simulate(net_formation_obs)
  G_true <- simulate(net_formation_true)

  # Estimate ERGM on G
  net_formation_est <- ergm(formula_est)
  
  # Drop the sender/receiver (xi) from the coefficients of the model
  xi_est <- c(1-sum(net_formation_est$coef[1:(n-1)]), net_formation_est$coef[1:(n-1)])
  net_formation_est$coef <- net_formation_est$coef[-(1:(n-1))]
  net_formation_est$formula <-  as.formula(paste("G_obs ~", paste(terms, collapse = " + ")))
  
  return(list(G_true = G_true, 
              G_alumni = G_alumni, 
              G_obs = G_obs,
              X = X,
              xi = xi,
              xi_est = xi_est,
              net_formation_est = net_formation_est, 
              net_formation_obs = net_formation_obs, 
              net_formation_alumni = net_formation_alumni, 
              net_formation_true = net_formation_true))
}



GenerateNetworkEndo <- function(params_net){
  n <- params_net$n
  terms_true <- params_net$terms_true
  terms_alumni <- params_net$terms_alumni
  terms_est <- params_net$terms_est
  target_stats_alumni <- params_net$target_stats_alumni
  target_stats_true <- params_net$target_stats_true
  target_stats_xi <- params_net$target_stats_xi
  xi <- params_net$xi
  
  # Create y and X
  X <- matrix(rnorm(n), n)
  
  formula_alumni <- as.formula(paste("G ~ ", paste(terms_alumni, collapse = " + ")))
  formula_obs <- as.formula(paste("G ~ offset(nodecov('xi')) + ", paste(terms_true, collapse = " + ")))
  formula_true <- as.formula(paste("G ~ ", paste(terms_true, collapse = " + ")))
  formula_est <- as.formula(paste("G_obs ~  sociality + ", paste(terms_est, collapse = " + ")))
  
  # Generate ERGM distribution using seed
  G <- network.initialize(n , directed = FALSE)
  G %v% 'xi' <- xi
  
  
  if (terms_alumni == "density"){
    G_alumni<-network(n, directed=FALSE , density= target_stats_alumni)
  } else{
    net_formation_alumni <- ergm(formula_alumni, offset.coef = target_stats_alumni)
    G_alumni <- simulate(net_formation_alumni)   
  }
  
  net_formation_obs <- ergm(formula_obs, offset.coef = c(target_stats_xi, target_stats_true))
  net_formation_true <- ergm(formula_true, offset.coef = target_stats_true)
  
  G_obs <- simulate(net_formation_obs)
  G_true <- simulate(net_formation_true)
  
  # Estimate ERGM on G
  net_formation_est <- ergm(formula_est)
  
  # Drop the sender/receiver (xi) from the coefficients of the model
  xi_est <- c(1-sum(net_formation_est$coef[1:(n-1)]), net_formation_est$coef[1:(n-1)])
  net_formation_est$coef <- net_formation_est$coef[-(1:(n-1))]
  net_formation_est$formula <-  as.formula(paste("G_obs ~", paste(terms_est, collapse = " + ")))
  
  return(list(G_true = G_true, 
              G_alumni = G_alumni, 
              G_obs = G_obs,
              X = X,
              xi = xi,
              xi_est = xi_est,
              net_formation_est = net_formation_est, 
              net_formation_obs = net_formation_obs, 
              net_formation_alumni = net_formation_alumni, 
              net_formation_true = net_formation_true))
}
