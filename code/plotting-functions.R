PlotStudent <- function(result, folder) {
  tval_phi = as.data.frame(result[[1]][[1]]$tval["phi",])[FALSE,]
  tval_alpha = as.data.frame(result[[1]][[1]]$tval["alpha",])[FALSE,]
  tval_beta = as.data.frame(result[[1]][[1]]$tval["beta_X",])[FALSE,]
  
  for (i in 1:length(result)) {
    a <- result[[i]]
    
    
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
    tval_phi[i,] = tval["phi",]
  }
  
  df_phi <- bind_cols(params, tval_phi)
  df_alpha <- bind_cols(params, tval_alpha)
  df_beta <- bind_cols(params, tval_beta)
  
  fig1 <- ggplot(data = df_phi) +
    geom_line(aes(phi, ergm, colour = "ergm")) +    
    geom_line(aes(phi, cor, colour = "cor")) +
    geom_line(aes(phi, true, colour = "true")) +
    geom_line(aes(phi, alumni, colour = "alumni")) +
    geom_line(aes(phi, obs, colour = "obs")) +
    geom_hline(yintercept = .9) +
    facet_wrap(c("beta_X", "psi")) +
    ylab("Coverage") +
    ggtitle("phi") +
    theme_minimal()
  fig1
  
  ggsave(paste(folder, "/phi.pdf", sep = ""))
  
  
  fig2 <- ggplot(data = df_beta) +
    geom_line(aes(phi, ergm, colour = "ergm")) +
    geom_line(aes(phi, cor, colour = "cor")) +
    geom_line(aes(phi, true, colour = "true")) +
    geom_line(aes(phi, alumni, colour = "alumni")) +
    geom_line(aes(phi, obs, colour = "obs")) +
    geom_hline(yintercept = .9) +
    facet_wrap(c("beta_X", "psi")) +
    ylab("Coverage") +
    ggtitle("beta") +
    theme_minimal()
  fig2
  
  ggsave(paste(folder, "/beta.pdf", sep = ""))
  
  fig3 <- ggplot(data = df_alpha) +
    geom_line(aes(phi, ergm, colour = "ergm")) +  
    geom_line(aes(phi, cor, colour = "cor")) +
    geom_line(aes(phi, true, colour = "true")) +
    geom_line(aes(phi, alumni, colour = "alumni")) +
    geom_line(aes(phi, obs, colour = "obs")) +
    geom_hline(yintercept = .9) +
    facet_wrap(c("beta_X", "psi")) +
    ylab("Coverage") +
    ggtitle("alpha") +
    theme_minimal()
  fig3
  
  ggsave(paste(folder, "/alpha.pdf", sep = ""))
}
