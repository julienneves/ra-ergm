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

####Test####
## Create network
# Set parameters
params_net <- list(xi = rnorm(50),
                   X = matrix(rbinom(50, 1, .5), ncol = 1),
                   terms_alumni = "density",
                   target_stats_alumni = 0.1,
                   terms = c("hamming(G_alumni)", "triangle", "edges", "nodematch('X')"),
                   target_stats_true = c(.4, .15, -2, -.5),
                   target_stats_xi = 1)

# Generate network
dgp_net <- GenerateNetwork(params_net)

dgp_net$net_formation_est

plot(simulate(dgp_net$net_formation_est))
plot(dgp_net$xi,dgp_net$xi_hat/params_net$target_stats_xi, ylim = c(-3,3), asp = 1)


params_model <- list(alpha = 0.1,
                     beta_X = 0.1,
                     phi = 0,
                     psi = 0.1,
                     sigma_e =  0,
                     B = 20)
fit <- EstimateModel(params_model, params_net, dgp_net)

repl <- 100

numCores <- detectCores()
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

result <- parLapply(cl, 1:repl, function(x, ...) tryCatch(EstimateModel(...), error = function(e) return(NULL)), params_model, params_net)

stopCluster(cl)

tval_phi = as.data.frame(result[[1]]$tval["phi",])[FALSE,]
tval_alpha = as.data.frame(result[[1]]$tval["alpha",])[FALSE,]
tval_beta = as.data.frame(result[[1]]$tval["beta_X",])[FALSE,]

a <- result
  
est = a[[1]]$est
std = a[[1]]$std
tval = abs(a[[1]]$tval) <= 1.645
  
  for (j in 2:length(a)) {
    
    if(is.null(a[[j]])){
      
      a[[j]] <- a[[j-1]]
    }
    est = est + a[[j]]$est
    std = std + a[[j]]$std
    tval = tval + (abs(a[[j]]$tval) <= 1.645)
  }
  
  est = est / length(a)
  std = std / length(a)
  tval = tval / length(a)
  
  tval_alpha = tval["alpha",]
  tval_beta = tval["beta_X",]
  tval_phi = tval["phi",]

tval_phi
tval_alpha
tval_beta