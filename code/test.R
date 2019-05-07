library(tidyr)
library(dplyr)
library(ggplot2)
library(parallel)

library(miscTools)
library(econet)
library(Matrix)
library(statnet)

coef_ergm <- as.data.frame(result[[1]][[1]]$coef_ergm)[FALSE,]
coef_nlls_true <- coef_ergm
coef_nlls_alumni <- coef_ergm
coef_nlls_obs <- coef_ergm
coef_true <- coef_ergm

mean_true <- coef_ergm
mean_ergm <- coef_ergm
mean_nlls_true <- coef_ergm
mean_nlls_alumni <- coef_ergm
mean_nlls_obs <- coef_ergm
mean_true <- coef_ergm

sd_ergm <- coef_ergm
sd_nlls_true <- coef_ergm
sd_nlls_alumni <- coef_ergm
sd_nlls_obs <- coef_ergm
sd_true <- coef_ergm

for (i in 1:length(result)) {
  a <- result[[i]]
  
  
  for (j in 1:length(a)) {
    coef_ergm[j, ] <- a[[j]]$coef_ergm
    coef_nlls_true[j, ] <- a[[j]]$coef_nlls_true
    coef_nlls_alumni[j, ] <- a[[j]]$coef_nlls_alumni
    coef_nlls_obs[j, ] <- a[[j]]$coef_nlls_obs
    coef_true[j, ] <- a[[j]]$coef_true
  }
  
  mean_ergm[i, ] <- as.data.frame(t(sapply(coef_ergm, mean, na.rm = TRUE)))
  sd_ergm[i, ] <- as.data.frame(t(sapply(coef_ergm, sd, na.rm = TRUE)))
  
  mean_nlls_true[i, ] <- as.data.frame(t(sapply(coef_nlls_true, mean, na.rm = TRUE)))
  sd_nlls_true[i, ] <- as.data.frame(t(sapply(coef_nlls_true, sd, na.rm = TRUE)))
  
  mean_nlls_alumni[i, ] <- as.data.frame(t(sapply(coef_nlls_alumni, mean, na.rm = TRUE)))
  sd_nlls_alumni[i, ] <- as.data.frame(t(sapply(coef_nlls_alumni, sd, na.rm = TRUE)))
  
  mean_nlls_obs[i, ] <- as.data.frame(t(sapply(coef_nlls_obs, mean, na.rm = TRUE)))
  sd_nlls_obs[i, ] <- as.data.frame(t(sapply(coef_nlls_obs, sd, na.rm = TRUE)))
  
  mean_true[i, ] <- as.data.frame(t(sapply(coef_true, mean, na.rm = TRUE)))
  
}  

names(params)[2] <- "beta_X"

names(mean_ergm) <- c("alpha_ergm", "beta_X_ergm", "phi_ergm")
names(mean_nlls_true) <- c("alpha_true", "beta_X_true", "phi_true")
names(mean_nlls_alumni) <- c("alpha_alumni", "beta_X_alumni", "phi_alumni")
names(mean_nlls_obs) <- c("alpha_obs", "beta_X_obs", "phi_obs")
names(sd_ergm) <- c("alpha_ergm_sd", "beta_X_ergm_sd", "phi_ergm_sd")
names(sd_nlls_true) <- c("alpha_true_sd", "beta_X_true_sd", "phi_true_sd")
names(sd_nlls_alumni) <- c("alpha_alumni_sd", "beta_X_alumni_sd", "phi_alumni_sd")
names(sd_nlls_obs) <- c("alpha_obs_sd", "beta_X_obs_sd", "phi_obs_sd")


df <- bind_cols(params, mean_ergm, mean_nlls_true, mean_nlls_alumni, mean_nlls_obs, sd_ergm, sd_nlls_true, sd_nlls_alumni, sd_nlls_obs)

fig1 <- ggplot(data = df) +
  geom_point(aes(phi+0.005, phi_ergm-phi, colour = "ergm")) +
  geom_point(aes(phi-0.005, phi_true-phi, colour = "true")) +
  geom_point(aes(phi+0.01, phi_alumni-phi, colour = "alumni"))+
  geom_point(aes(phi-0.01, phi_obs-phi, colour = "obs")) +
  geom_errorbar(aes(phi+0.005, ymin = (phi_ergm - phi)- 1.96 * phi_ergm_sd, ymax = (phi_ergm - phi) + 1.96 * phi_ergm_sd, colour = "ergm")) +
  geom_errorbar(aes(phi-0.005, ymin = (phi_true-phi) - 1.96 * phi_true_sd, ymax = (phi_true-phi)+  1.96 * phi_true_sd, colour = "true"))  +
  geom_errorbar(aes(phi+0.01, ymin = (phi_alumni-phi) - 1.96 * phi_alumni_sd, ymax = (phi_alumni-phi)+ 1.96 * phi_alumni_sd, colour = "alumni"))  +
  geom_errorbar(aes(phi-0.01, ymin =  (phi_obs-phi)- 1.96 * phi_obs_sd, ymax = (phi_obs-phi)+  1.96 * phi_obs_sd, colour = "obs")) +
  ylab("phi-phi0") + theme_minimal()
fig1
