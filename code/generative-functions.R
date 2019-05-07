GenerateData <- function(params_model, dgp_net){
  # Extract parameters
  alpha <- params_model$alpha
  beta <-  params_model$beta_X
  phi <-  params_model$phi
  sigma_e <-  params_model$sigma_e
  
  G_true <- dgp_net$G_true
  X <- dgp_net$X
  
  n <-  G_true$gal$n
  epsilon <- rnorm(n, 0, sigma_e)
  
  # Create y
  y <- solve((diag(n) - phi * as.matrix(G_true)),  alpha + X %*% beta + epsilon)
  
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

  # Create y and X
  X <- matrix(rnorm(n), n)
  
  formula_alumni <- as.formula(paste("G ~ ", paste(terms_alumni, collapse = " + ")))
  formula_true <- as.formula(paste("G ~ ", paste(terms_true, collapse = " + ")))
  formula_est <- as.formula(paste("G_obs ~ ", paste(terms_true, collapse = " + ")))

  G <- network.initialize(n , directed = FALSE)

  # Generate ERGM distribution using seed
  
  if (target_stats_alumni == "density"){
    G_alumni<-network(n, directed=FALSE, density=0.1)
  } else{
    net_formation_alumni <- ergm(formula_alumni, target.stats = target_stats_alumni)
    G_alumni <- simulate(net_formation_alumni)   
  }
  
  net_formation_true <- ergm(formula_true, target.stats = target_stats_true)
  G_obs <- simulate(net_formation_true)
  G_true <- simulate(net_formation_true)

  # Estimate ERGM on G
  net_formation_est <- ergm(formula_est)
  
  return(list(G_true = G_true, 
              G_alumni = G_alumni, 
              G_obs = G_obs,
              X = X,
              net_formation_est = net_formation_est, 
              net_formation_alumni = net_formation_alumni, 
              net_formation_true = net_formation_true))
}