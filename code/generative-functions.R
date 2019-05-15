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
  xi <- params_net$xi
  X <- params_net$X
  n <- length(xi)
  
  target_stats_alumni <- params_net$target_stats_alumni
  target_stats_true <- params_net$target_stats_true
  target_stats_xi <- params_net$target_stats_xi
  
  terms_alumni <- params_net$terms_alumni
  terms_est <- params_net$terms
  terms_true <- paste("offset(",paste(terms_est, ")", sep = ""), sep = "")

  
  formula_alumni <- as.formula(paste("G ~ ", paste(terms_alumni, collapse = " + ")))
  formula_true <- as.formula(paste("G ~ ", paste(terms_true, collapse = " + ")))
  formula_obs <- as.formula(paste("G ~ offset(nodecov('xi')) + ", paste(terms_true, collapse = " + ")))
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
  net_formation_est <- ergm(formula_est, estimate = "MPLE")
  
  # Drop the sender/receiver (xi) from the coefficients of the model
  xi_hat <- c(1-sum(net_formation_est$coef[1:(n-1)]), net_formation_est$coef[1:(n-1)])
  net_formation_est$coef <- net_formation_est$coef[-(1:(n-1))]
  net_formation_est$MCMCtheta <- net_formation_est$MCMCtheta[-(1:(n-1))]
  net_formation_est$formula <-  as.formula(paste("G_obs ~", paste(terms_est, collapse = " + ")))
  
  return(list(G_true = G_true, 
              G_alumni = G_alumni, 
              G_obs = G_obs,
              X = X,
              xi = xi,
              xi_hat = xi_hat,
              net_formation_est = net_formation_est, 
              net_formation_obs = net_formation_obs,  
              net_formation_true = net_formation_true))
}
