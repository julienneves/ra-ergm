library(tidyr)
library(dplyr)
library(ggplot2)
library(parallel)

library(miscTools)
library(econet)
library(statnet)

source('~/Dropbox/Work/Research/ergm/code/estimation-functions.R')
source('~/Dropbox/Work/Research/ergm/code/generative-functions.R')
source('~/Dropbox/Work/Research/ergm/code/ploting-functions.R')

## Parameters
params_net <- list(n = 50,
                   target_stats = c(10,100), 
                   terms = c("triangle", "edges"))

# Generate network
dgp_net <- GenerateNetwork(params_net)

# Plot true and observed network
png("output/network_03_21.png")
par(mfrow=c(1,2)) 
plot(dgp_net$G_true, coord = expand.grid(1:ceiling(sqrt(params_net$n)),1:ceiling(sqrt(params_net$n)))[1:params_net$n,])
title("True Network")
plot(dgp_net$G_obs,  coord = expand.grid(1:ceiling(sqrt(params_net$n)),1:ceiling(sqrt(params_net$n)))[1:params_net$n,])
title("Observed Network")
dev.off()


# Example 1
params_model_1 <- list(alpha = .5, 
                       beta_X = .5, 
                       phi = 0,
                       psi = 1,
                       sigma_e = 1)
result_1 <- mclapply(1:100, function(x, ...) EstimateModel(...), dgp_net, params_model_1)
fig_1 <- SimulationPlot(result_1)
ggsave("output/example_1_03_21.png")
fig_1

# Example 2
params_model_2 <- list(alpha = .5, 
                       beta_X = .5, 
                       phi = 0.8,
                       psi = 1,
                       sigma_e = 1)
result_2 <- mclapply(1:100, function(x, ...) EstimateModel(...), dgp_net, params_model_2)
fig_2 <- SimulationPlot(result_2)
ggsave("output/example_2_03_21.png")
fig_2