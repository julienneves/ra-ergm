PlotStudent <- function(result, params, folder) {
  tval_phi = as.data.frame(result[[1]][[2]]$tval["phi",])[FALSE,]
  tval_alpha = as.data.frame(result[[1]][[2]]$tval["alpha",])[FALSE,]
  tval_beta = as.data.frame(result[[1]][[2]]$tval["beta_X",])[FALSE,]
  
  for (i in 1:length(result)) {
    
    a <- result[[i]]
    
    a <- a[-which(sapply(a, is.null))]

    est = a[[1]]$est
    std = a[[1]]$std
    tval = abs(a[[1]]$tval) <= 1.645
    
    for (j in 2:length(a)) {
      est = est + a[[j]]$est
      std = std + a[[j]]$std
      tval = tval + (abs(a[[j]]$tval) <= 1.645)
    }
    
    est = est / length(a)
    std = std / length(a)
    tval = tval / length(a)
    
    tval_alpha[i, ] = tval["alpha",]
    tval_beta[i, ] = tval["beta_X",]
    tval_phi[i,] = tval["phi",]}
  
  df_phi  <- bind_cols(params, tval_phi)
  df_alpha <- bind_cols(params, tval_alpha)
  df_beta <- bind_cols(params, tval_beta)
  
  fig1 <- ggplot(data = df_phi ) +
    geom_smooth(aes(phi, ergm, colour = "ergm")) +    
    geom_smooth(aes(phi, cor, colour = "cor")) +
    geom_smooth(aes(phi, true, colour = "true")) +
    geom_smooth(aes(phi, obs, colour = "obs")) +
    geom_hline(yintercept = .9) +
    facet_wrap(~ psi + sigma_e)  +
    ylab("Coverage") +
    ggtitle("phi") +
    theme_minimal()
  fig1
  
  ggsave(paste(folder, "/phi.pdf", sep = ""))
  
  
  fig2 <- ggplot(data = df_beta) +
    geom_smooth(aes(phi, ergm, colour = "ergm")) +
    geom_smooth(aes(phi, cor, colour = "cor")) +
    geom_smooth(aes(phi, true, colour = "true")) +
    geom_smooth(aes(phi, obs, colour = "obs")) +
    geom_hline(yintercept = .9) +
    facet_wrap(~ psi + sigma_e)  +
    ylab("Coverage") +
    ggtitle("beta") +
    theme_minimal()
  fig2
  
  ggsave(paste(folder, "/beta.pdf", sep = ""))
  
  fig3 <- ggplot(data = df_alpha) +
    geom_smooth(aes(phi, ergm, colour = "ergm")) +  
    geom_smooth(aes(phi, cor, colour = "cor")) +
    geom_smooth(aes(phi, true, colour = "true")) +
    geom_smooth(aes(phi, obs, colour = "obs")) +
    geom_hline(yintercept = .9) +
    facet_wrap(~ psi + sigma_e)  +
    ylab("Coverage") +
    ggtitle("alpha") +
    theme_minimal()
  fig3
  
  ggsave(paste(folder, "/alpha.pdf", sep = ""))
}
