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
params_net <- list(xi = rnorm(50),
                   X = matrix(rbinom(50, 1, .5), ncol = 1),
                   terms_alumni = "density",
                   target_stats_alumni = 0.1,
                   terms = c("hamming(G_alumni)", "triangle", "edges", "nodematch('X')"),
                   target_stats_true = c(.4, .15, -2, -.5),
                   target_stats_xi = 1)

dgp_net <- GenerateNetwork(params_net)

dgp_net$net_formation_est

plot(simulate(dgp_net$net_formation_est))
plot(dgp_net$xi,dgp_net$xi_hat/params_net$target_stats_xi, ylim = c(-3,3), asp = 1)
plot(dgp_net$G_true)

## Run Monte Carlo simulation
# Set parameters grid
alpha <- c(0,0.2)
beta <- c(0,0.2)
phi <- c(-0.2, -0.1, 0, 0.1, 0.2)
psi <- c(0,0.1, 1)
sigma_e <- c(0.001, 0.1, 1)
B <- 50

# Set number of replications
repl <- 100

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
                    "GenerateNetwork",
                    "EstimateNLLS", 
                    "EstimateERGM"))

  params_model <- list(alpha = params[i,"alpha"], 
                       beta_X = params[i,"beta_X"], 
                       phi = params[i,"phi"],
                       psi = params[i, "psi"],
                       sigma_e =  params[i,"sigma_e"],
                       B =  params[i,"B"])
  result[[i]] <- parLapply(cl, 1:repl, function(x, ...) tryCatch(EstimateModel(...), error = function(e) return(NULL)), params_model, params_net)
  cat("Trial",i,"Out of",nrow(params), "\n")
  save.image(paste(folder, "/result.RData", sep = ""))
  
  stopCluster(cl)
}


# Plot results
PlotStudent(result, folder)