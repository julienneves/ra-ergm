library(tidyr)
library(dplyr)
library(ggplot2)
library(parallel)

library(miscTools)
library(econet)
library(Matrix)
library(statnet)

coef_ergm <- as.data.frame(result[[1]][[1]]$coef_ergm)[FALSE,]
coef_nlls <- coef_ergm

mean_ergm <- coef_ergm
sd_ergm <- coef_ergm

mean_nlls <- coef_ergm
sd_nlls <- coef_ergm

coef_true <- coef_ergm
mean_true <- coef_ergm

for (i in 1:length(result)) {
  a <- result[[i]]
  
  
  for (j in 1:length(a)) {
    coef_ergm[j, ] <- a[[j]]$coef_ergm
    coef_nlls[j, ] <- a[[j]]$coef_nlls
    coef_true[j, ] <- a[[j]]$coef_true
  }
  
  mean_ergm[i, ] <- as.data.frame(t(sapply(coef_ergm, mean, na.rm = TRUE)))
  sd_ergm[i, ] <- as.data.frame(t(sapply(coef_ergm, sd, na.rm = TRUE)))
  
  mean_nlls[i, ] <- as.data.frame(t(sapply(coef_nlls, mean, na.rm = TRUE)))
  sd_nlls[i, ] <- as.data.frame(t(sapply(coef_nlls, sd, na.rm = TRUE)))
  
  mean_true[i, ] <- as.data.frame(t(sapply(coef_true, mean, na.rm = TRUE)))
  
}  

names(params)[2] <- "beta_X"

names(mean_ergm) <- c("alpha_ergm", "beta_X_ergm", "phi_ergm")
names(mean_nlls) <- c("alpha_nlls", "beta_X_nlls", "phi_nlls")
names(sd_ergm) <- c("alpha_ergm_sd", "beta_X_ergm_sd", "phi_ergm_sd")
names(sd_nlls) <- c("alpha_nlls_sd", "beta_X_nlls_sd", "phi_nlls_sd")

df <- bind_cols(params, mean_ergm, mean_nlls, sd_ergm, sd_nlls)

fig1 <- ggplot(data = df) +   
  geom_hline(aes(yintercept= phi)) +
  geom_point(aes(beta_X, phi_ergm, colour = "ergm")) +
  geom_point(aes(beta_X, phi_nlls, colour = "nlls"))  +
  geom_errorbar(aes(beta_X, ymin = phi_ergm - 1.96 * phi_ergm_sd, ymax = phi_ergm + 1.96 * phi_ergm_sd, colour = "ergm")) +
  geom_errorbar(aes(beta_X, ymin = phi_nlls - 1.96 * phi_nlls_sd, ymax = phi_nlls + 1.96 * phi_nlls_sd, colour = "nlls"))  +
  facet_wrap(c("phi", "psi"), ncol = 3) +
  scale_colour_manual(name="Type", values=c(ergm="red", nlls="blue"))+
  coord_cartesian( ylim = c(-0.3, .3)) +
  ylab("phi") + theme_minimal()
fig1

fig2 <- ggplot(data = df) +   
  geom_hline(aes(yintercept= beta_X)) +
  geom_point(aes(phi, beta_X_ergm, colour = "ergm")) +
  geom_point(aes(phi, beta_X_nlls, colour = "nlls")) +
  geom_errorbar(aes(phi, ymin = beta_X_ergm - 1.96 * beta_X_ergm_sd, ymax = beta_X_ergm + 1.96 * beta_X_ergm_sd, colour = "ergm")) +
  geom_errorbar(aes(phi, ymin = beta_X_nlls - 1.96 * beta_X_nlls_sd, ymax = beta_X_nlls + 1.96 * beta_X_nlls_sd, colour = "nlls"))  +
  facet_wrap(c("beta_X", "psi"), ncol = 3) +
  scale_colour_manual(name="Type", values=c(ergm="red", nlls="blue"))+
  coord_cartesian( ylim = c(-1, 12)) +
  ylab("beta") + theme_minimal()
fig2

fig3 <- ggplot(data = df) +   
  geom_hline(aes(yintercept= alpha)) +
  geom_point(aes(phi, alpha_ergm, colour = "ergm")) +
  geom_point(aes(phi, alpha_nlls, colour = "nlls")) +
  geom_errorbar(aes(phi, ymin = alpha_ergm - 1.96 * alpha_ergm_sd, ymax = alpha_ergm + 1.96 * alpha_ergm_sd, colour = "ergm")) +
  geom_errorbar(aes(phi, ymin = alpha_nlls - 1.96 * alpha_nlls_sd, ymax = alpha_nlls + 1.96 * alpha_nlls_sd, colour = "nlls"))  +
  facet_wrap(c("beta_X", "psi"), ncol = 3) +
  scale_colour_manual(name="Type", values=c(ergm="red", nlls="blue"))+
  coord_cartesian( ylim = c(-1, 1)) +
  ylab("alpha") + theme_minimal()
fig3
