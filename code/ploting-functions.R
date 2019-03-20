SimulationPlot <- function(result){
  coef_ergm <- t(sapply(result, function(x) colMeans(x$coef_ergm)))
  coef_nlls <- t(sapply(result, function(x) colMeans(x$coef_nlls)))
  coef_true <- t(sapply(result, function(x) colMeans(x$coef_true)))
  
  fig <- ggplot(NULL, aes(value)) + 
  geom_density(aes(value, fill = "ergm"), data = gather(as.data.frame(coef_ergm)), kernel = "gaussian", alpha = 0.4) +
  geom_density(aes(value, fill = "nlls"), data = gather(as.data.frame(coef_nlls)), kernel = "gaussian", alpha = 0.2) +
  geom_vline(aes(xintercept = value), data = gather(as.data.frame(coef_true))) +
  facet_wrap(~ key, scales = "free_x") + 
  theme(legend.position = "right")
  fig
  return(fig)
}
