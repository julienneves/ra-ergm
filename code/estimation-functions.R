EstimateModel <- function(dgp_net, params_model) {
  
  data <- GenerateData(params_model, dgp_net)
  
  coef_true <- params_model[c("alpha", "beta_X", "phi")]
  
  cat(sprintf("True network\n"))
  fit_true <- EstimateNLLS(data, dgp_net, start_val = coef_true, type = "true")
  
  cat(sprintf("Alumni network\n"))
  fit_alumni <- EstimateNLLS(data, dgp_net, start_val = coef_true, type = "alumni")
  
  cat(sprintf("Observed network\n"))
  fit_obs <- EstimateNLLS(data, dgp_net, start_val = coef_true, type = "obs")
  
  cat(sprintf("ERGM network\n"))
  fit_ergm <- tryCatch({EstimateERGM(data, dgp_net, start_val = coef_true, repl = params_model$B)})
  
  return(list(fit_ergm = summary(fit_ergm), 
              fit_true = summary(fit_true), 
              fit_obs = summary(fit_obs), 
              fit_alumni = summary(fit_alumni), 
              coef_true = as.data.frame(coef_true)))
}

EstimateNLLS <- function(data, dgp_net, start_val, type) {
  
  if (type == "true"){
    G <- as.matrix(dgp_net$G_true)
  } else if (type == "obs") {
    G <- as.matrix(dgp_net$G_obs)
  } else if (type == "alumni"){
    G <- as.matrix(dgp_net$G_alumni)
  } else {
    warning("No G specified")
  }
  
  fit<- net_dep(formula = "y ~ X", data = data, G = G, 
                      model = "model_B", estimation = "NLLS",  hypothesis = "lim", 
                      start.val = start_val, endogeneity = FALSE)
  return(fit)
}


EstimateERGM <- function(data, dgp_net, start_val, repl = 99) {

  G_obs <- dgp_net$G_obs
  G_alumni <- dgp_net$G_alumni
  
  df <- bind_rows(replicate(repl, data, simplify = FALSE)) 
  
  G <- simulate(dgp_net$net_formation_est, nsim = repl)
  G <- as.matrix(bdiag(lapply(G, as.matrix)))
  
  fit <- net_dep(formula = "y ~ X", data = df, G = G, 
                      model = "model_B", estimation = "NLLS",  hypothesis = "lim", 
                      start.val = start_val, endogeneity = FALSE)
  return(fit)
}

