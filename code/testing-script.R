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
source('~/GitHub/ra-ergm/code/plotting-functions.R')

folder <- paste("output/result", Sys.Date(), sep = "-")
dir.create(folder)

#####
## Create network
# Set parameters
params_net <- list(xi = rnorm(100),
                   X = matrix(rnorm(100), ncol = 1),
                   terms_alumni = "density",
                   target_stats_alumni = 0.2,
                   terms = c("hamming(G_alumni)", "triangle", "edges"),
                   target_stats_true = c(.4, -.15, -1.7),
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
                     B = 50)
fit <- EstimateModel(params_model, params_net, dgp_net)


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
  
  result <- parLapply(cl, 1:repl, function(x, ...) EstimateModel(...), params_model, params_net, dgp_net)
  cat("Trial",i,"Out of",nrow(params), "\n")
  save.image("~/GitHub/ra-ergm/output/result_5_8.RData")
  
  stopCluster(cl)
}


# Plot results
PlotStudent(result, folder)