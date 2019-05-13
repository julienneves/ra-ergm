#####
library(tidyverse)

library(econet)
library(statnet)

library(parallel)
library(miscTools)
library(Matrix)

# source('~/Dropbox/Work/Research/ergm/code/estimation-functions.R')
# source('~/Dropbox/Work/Research/ergm/code/generative-functions.R')
# source('~/Dropbox/Work/Research/ergm/code/ploting-functions.R')

source('~/GitHub/ra-ergm/code/estimation-functions.R')
source('~/GitHub/ra-ergm/code/generative-functions.R')
source('~/GitHub/ra-ergm/code/ploting-functions.R')

folder <- paste("output/result", Sys.Date(), sep = "-")
dir.create(folder)

#####
## Create network
# Set parameters
params_net <- list(n = 50,
                   xi = rnorm(50),
                   target_stats_alumni = 0.2,
                   terms_alumni = "density",
                   target_stats_true = c(.4, -.15, -1.7),
                   terms = c("hamming(G_alumni)", "triangle", "edges"),
                   target_stats_xi = 1)

# Generate network
dgp_net <- GenerateNetwork(params_net)

dgp_net$net_formation_est

plot(simulate(dgp_net$net_formation_est))
plot(dgp_net$xi,dgp_net$xi_hat/params_net$target_stats_xi, ylim = c(-3,3), asp = 1)


# Plot true and observed network
png("output/network_05_8.png")
par(mfrow=c(1,3)) 
plot(dgp_net$G_alumni,  coord = expand.grid(1:ceiling(sqrt(params_net$n)),1:ceiling(sqrt(params_net$n)))[1:params_net$n,])
title("Alumni Network")
plot(dgp_net$G_true, coord = expand.grid(1:ceiling(sqrt(params_net$n)),1:ceiling(sqrt(params_net$n)))[1:params_net$n,])
title("True Network")
plot(dgp_net$G_obs,  coord = expand.grid(1:ceiling(sqrt(params_net$n)),1:ceiling(sqrt(params_net$n)))[1:params_net$n,])
title("Observed Network")
dev.off()

params_model <- list(alpha = 0,
                     beta_X = 0,
                     phi = 0.2,
                     psi = 0,
                     sigma_e =  0.2,
                     B = 100)
fit <- EstimateModel(params_model, params_net, dgp_net)

#####
## Run Monte Carlo simulation
# Set parameters grid
alpha <- 0
beta <- c(0)
phi <- c(-0.2, -0.1, 0, 0.1, 0.2)
psi <- c(0)
sigma_e <- 0.1
B <- 100

# Set number of replications
repl <- 100

# Set parameters
params_net <- list(n = 50,
                   xi = rnorm(50),
                   target_stats_alumni = 0.2,
                   terms_alumni = "density",
                   target_stats_true = c(.4, -.15, -1.7),
                   terms = c("hamming(G_alumni)", "triangle", "edges"),
                   target_stats_xi = 1)

# Create a data frame to hold parameters and p-values
params <- expand.grid(alpha, beta, phi, psi, sigma_e, B)
colnames(params) <-c("alpha", "beta_X", "phi", "psi", "sigma_e", "B")
params <- tibble::rownames_to_column(params,"trial")

result <- vector("list", nrow(params)) 


numCores <- detectCores()


for (i in 1:nrow(params)){
cl <- makeCluster(numCores)
clusterEvalQ(cl, {
  library(tidyverse)
  library(econet)
  library(statnet)
  
  library(miscTools)
  library(Matrix)
})

clusterExport(cl, c("EstimateModel", 
                    "GenerateData", 
                    "EstimateNLLS", 
                    "EstimateERGM"))

  params_model <- list(alpha = params[i,"alpha"], 
                       beta_X = params[i,"beta_X"], 
                       phi = params[i,"phi"],
                       psi = params[i, "psi"],
                       sigma_e =  params[i,"sigma_e"],
                       B =  params[i,"B"])
  result[[i]] <- parLapply(cl, 1:repl, function(x, ...) EstimateModel(...), params_model, params_net, dgp_net)
  cat("Trial",i,"Out of",nrow(params), "\n")
  save.image("~/GitHub/ra-ergm/output/result_5_8.RData")
  
  stopCluster(cl)
}


# Plot results
PlotStat(result, folder)