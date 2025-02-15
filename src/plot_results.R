library(ggplot2)


plot_gamma_estimates <- function(results_df, save = NULL) {
    # Ensure 'method' is treated as factor in the order encountered
    results_df$method <- factor(results_df$method, levels = unique(results_df$method))
  
    # Create and return the plot
    p <- ggplot(results_df, aes(x = gamma_hat, y = method)) +
        geom_errorbarh(aes(xmin = LCL, xmax = UCL), height = 0.2) +
        geom_point(color = "red") +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
        theme_bw() +
        xlab("Point Estimate") +
        ylab("Method")
    
    if (!is.null(save)) {
        ggsave(paste("results/", save, ".png", sep = ""), plot = p)
    }
    
    return(p)
}

