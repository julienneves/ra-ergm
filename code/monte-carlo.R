library(tidyr)
library(dplyr)
library(ggplot2)
library(parallel)

library(miscTools)
library(econet)
library(Matrix)
library(statnet)

source('~/GitHub/ra-ergm/code/estimation-functions.R')
source('~/GitHub/ra-ergm/code/generative-functions.R')
source('~/GitHub/ra-ergm/code/ploting-functions.R')

# Set seed
set.seed(123)

## Create network
# Set parameters
params_net <- list(n = 50,
                   target_xi = 85,
                   target_stats = c(20,150), 
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

## Run Monte Carlo simulation
# Set parameters grid
alpha <- 0
beta <- c(0,1,10)
phi <- c(-0.2, -0.1, 0, 0.1, 0.2)
psi <- c(1, 10, 100)
sigma_e <- 1

# Set number of replications
repl <- 100

# Create a data frame to hold parameters and p-values
params <- expand.grid(alpha, beta, phi, psi, sigma_e)
colnames(params) <-c("alpha", "beta", "phi", "psi", "sigma_e")

result <- vector("list", nrow(params)) 


numCores <- detectCores()
cl <- makeCluster(numCores)
clusterEvalQ(cl, {
  library(tidyr)
  library(dplyr)
  library(ggplot2)
  library(parallel)
  
  library(miscTools)
  library(econet)
  library(Matrix)
  library(statnet)
})

clusterExport(cl, {"EstimateModel"
  "GenerateData"
  "EstimateNLLS"
  "EstimateERGM"})


for (i in 1:nrow(params)){
  params_model <- list(alpha = params[i,"alpha"], 
                       beta_X = params[i,"beta"], 
                       phi = params[i,"phi"],
                       psi = params[i,"psi"],
                       sigma_e =  params[i,"sigma_e"])
  result[[i]] <- parLapply(cl, 1:repl, function(x, ...) EstimateModel(...), dgp_net, params_model)
  cat("Trial",i,"Out of",nrow(params))
}

## Print results

