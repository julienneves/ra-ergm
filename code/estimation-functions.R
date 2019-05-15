EstimateModel <- function(params_model,params_net, dgp_net = NULL) {
  
  if (is.null(dgp_net)){
    dgp_net <- GenerateNetwork(params_net)
    data <- GenerateData(params_model, dgp_net)
  } else {
    data <- GenerateData(params_model, dgp_net)
  }
  
  
  coef_true <- params_model[c("alpha", "beta_X", "phi")]
  start_val <- coef_true
  
  cat(sprintf("ERGM network\n"))
  fit_ergm <- EstimateERGM(data, dgp_net, start_val = coef_true, repl = params_model$B)
  
  cat(sprintf("True network\n"))
  fit_true <- EstimateNLLS(data, dgp_net, start_val = coef_true, type = "true")
  
  cat(sprintf("Alumni network\n"))
  fit_alumni <- EstimateNLLS(data, dgp_net, start_val = coef_true, type = "alumni")
  
  
  cat(sprintf("Observed network\n"))
  fit_obs <- EstimateNLLS(data, dgp_net, start_val = coef_true, type = "obs")

  
  cat(sprintf("Observed with correction\n"))
  start_val$beta_xi <- params_model$psi
  fit_cor <- EstimateNLLS(data, dgp_net, start_val = coef_true, type = "obs", correction = TRUE)
  
  est <- data.frame(ergm = summary(fit_ergm)$coefficients[,"Estimate"],
                    true = summary(fit_true)$coefficients[,"Estimate"],
                    obs = summary(fit_obs)$coefficients[,"Estimate"],
                    alumni = summary(fit_alumni)$coefficients[,"Estimate"],
                    cor = summary(fit_cor)$coefficients[1:3,"Estimate"])
  std <-  data.frame(ergm = summary(fit_ergm)$coefficients[,"Std. Error"],
                     true = summary(fit_true)$coefficients[,"Std. Error"],
                     obs = summary(fit_obs)$coefficients[,"Std. Error"],
                     alumni = summary(fit_alumni)$coefficients[,"Std. Error"],
                     cor = summary(fit_cor)$coefficients[1:3,"Std. Error"])
  
  
  coef_true <- data.frame(params = unlist(coef_true))
  
  tval <- (est-rep(coef_true,ncol(est)))/std
  
  return(list(est = est,
              std = std,
              tval = tval))
}

EstimateNLLS <- function(data, dgp_net, start_val, type, correction = NULL) {
  
  if (type == "true"){
    G <- as.matrix(dgp_net$G_true)
  } else if (type == "obs") {
    G <- as.matrix(dgp_net$G_obs)
  } else if (type == "alumni"){
    G <- as.matrix(dgp_net$G_alumni)
  } else {
    warning("No G specified")
  }
  
  
  if (is.null(correction)){
    fit<- net_dep(formula = "y ~ X", data = data, G = G, 
                  model = "model_B", estimation = "NLLS",  hypothesis = "lim", 
                  start.val = start_val, endogeneity = FALSE)
  } else {
    data <- bind_cols(data, xi = dgp_net$xi_hat)
    
    fit<- net_dep(formula = "y ~ X", data = data, G = G, 
                  model = "model_B", estimation = "NLLS",  hypothesis = "lim", 
                  start.val = start_val, endogeneity = FALSE)   
  }
  return(fit)
}


EstimateERGM <- function(data, dgp_net, start_val, repl = 99, correction = NULL) {

  G_obs <- dgp_net$G_obs
  G_alumni <- dgp_net$G_alumni
  
  G <- simulate(dgp_net$net_formation_est, nsim = repl)
  G <- as.matrix(bdiag(lapply(G, as.matrix)))
  
  if (is.null(correction)){
    df <- bind_rows(replicate(repl, data, simplify = FALSE)) 
    
    
    fit <- net_dep(formula = "y ~ X", data = df, G = G, 
                   model = "model_B", estimation = "NLLS",  hypothesis = "lim", 
                   start.val = start_val, endogeneity = FALSE)
  } else {
    data <- bind_cols(data, xi = dgp_net$xi_hat)
    df <- bind_rows(replicate(repl, data, simplify = FALSE)) 
    
    fit <- net_dep(formula = "y ~ X + xi", data = df, G = G, 
                   model = "model_B", estimation = "NLLS",  hypothesis = "lim", 
                   start.val = start_val, endogeneity = FALSE)    
  }
  
  return(fit)
}

