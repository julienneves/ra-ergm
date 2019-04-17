EstimateModel <- function(dgp_net, params_model) {
  
  data <- GenerateData(params_model, dgp_net)
  
  coef_true <- params_model[c("alpha", "beta_X", "phi")]
  coef_nlls <- EstimateNLLS(data, dgp_net, start_val = coef_true)
  coef_ergm <- tryCatch({EstimateERGM(data, dgp_net, start_val = coef_true)})
  
  return(list(coef_ergm = coef_ergm, coef_nlls = coef_nlls, coef_true = as.data.frame(coef_true)))
}

EstimateNLLS <- function(data, dgp_net, start_val) {
  
  G <- as.matrix(dgp_net$G_obs)
  
  fit<- net_dep(formula = "y ~ X", data = data, G = G, 
                      model = "model_B", estimation = "NLLS",  hypothesis = "lim", 
                      start.val = start_val, endogeneity = FALSE)
  coef_nlls <- t(coef(fit[[1]]))
  return(coef_nlls)
}

EstimateERGM <- function(data, dgp_net, start_val, repl = 99) {

  G_obs <- dgp_net$G_obs
  
  df <- bind_rows(replicate(repl, data, simplify = FALSE)) 
  
  G <- simulate(dgp_net$net_formation_est, nsim = repl)
  G <- as.matrix(bdiag(lapply(G, as.matrix)))
  
  fit <- net_dep(formula = "y ~ X", data = df, G = G, 
                      model = "model_B", estimation = "NLLS",  hypothesis = "lim", 
                      start.val = start_val, endogeneity = FALSE)
  coef_ergm <- t(coef(fit[[1]]))
  return(coef_ergm)
}

