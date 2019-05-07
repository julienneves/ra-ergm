library(tidyverse)

library(parallel)

library(miscTools)
library(econet)
library(Matrix)
library(statnet)

# source('~/Dropbox/Work/Research/ergm/code/estimation-functions.R')
# source('~/Dropbox/Work/Research/ergm/code/generative-functions.R')
# source('~/Dropbox/Work/Research/ergm/code/ploting-functions.R')

source('~/GitHub/ra-ergm/code/estimation-functions.R')
source('~/GitHub/ra-ergm/code/generative-functions.R')
source('~/GitHub/ra-ergm/code/ploting-functions.R')

# Set seed
set.seed(123)

## Create network
# Set parameters
params_net <- list(n = 50,
                   target_stats_alumni = c(0.1),
                   terms_alumni = c("density"),
                   target_stats_true = c(40, 30, 160), 
                   terms_true = c("hamming(G_alumni)", "triangle", "edges"))

# Generate network
dgp_net <- GenerateNetwork(params_net)

# Plot true and observed network
png("output/network_04_25.png")
par(mfrow=c(1,3)) 
plot(dgp_net$G_alumni,  coord = expand.grid(1:ceiling(sqrt(params_net$n)),1:ceiling(sqrt(params_net$n)))[1:params_net$n,])
title("Alumni Network")
plot(dgp_net$G_true, coord = expand.grid(1:ceiling(sqrt(params_net$n)),1:ceiling(sqrt(params_net$n)))[1:params_net$n,])
title("True Network")
plot(dgp_net$G_obs,  coord = expand.grid(1:ceiling(sqrt(params_net$n)),1:ceiling(sqrt(params_net$n)))[1:params_net$n,])
title("Observed Network")
dev.off()


params_model <- list(alpha = 0, 
                     beta_X = 0.2, 
                     phi = .5,
                     sigma_e =  1,
                     B = 99)
fit <- EstimateModel(dgp_net, params_model)


## Run Monte Carlo simulation
# Set parameters grid
alpha <- 0
beta <- c(0.5)
phi <- c(-0.2, -0.1, 0, 0.1, 0.2)
sigma_e <- c(1)
B <- 99

# Set number of replications
repl <- 250

# Create a data frame to hold parameters and p-values
params <- expand.grid(alpha, beta, phi, sigma_e, B)
colnames(params) <-c("alpha", "beta_X", "phi", "sigma_e", "B")
params <- tibble::rownames_to_column(params,"trial")

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

clusterExport(cl, c("EstimateModel", 
                    "GenerateData", 
                    "EstimateNLLS", 
                    "EstimateERGM"))


for (i in 1:nrow(params)){
  params_model <- list(alpha = params[i,"alpha"], 
                       beta_X = params[i,"beta_X"], 
                       phi = params[i,"phi"],
                       sigma_e =  params[i,"sigma_e"],
                       B =  params[i,"B"])
  result[[i]] <- parLapply(cl, 1:repl, function(x, ...) EstimateModel(...), dgp_net, params_model)
  cat("Trial",i,"Out of",nrow(params))
}
