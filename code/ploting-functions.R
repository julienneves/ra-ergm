SimulationPlot <- function(result){
  fig <- ggplot(NULL, aes(value)) + 
    geom_density(aes(value, y=..scaled.., fill = "ergm"), data = gather(coef_ergm), alpha = 0.4) +
    geom_density(aes(value, y=..scaled.., fill = "nlls"), data = gather(coef_nlls), alpha = 0.2) +
    geom_vline(aes(xintercept = value), data = gather(coef_true)) +
    facet_wrap(~ key + value, scales = "free_x") +
    theme(legend.position = "right")
  return(fig)
}
