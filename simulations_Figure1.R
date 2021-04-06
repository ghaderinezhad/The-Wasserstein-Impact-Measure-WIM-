#' Poisson model: simulations for Figure 1
#' comparison of a Gamma(2.5, 2.5) prior with a Gamma(alpha, beta) prior

library(foreign)
rm(list = ls())
library(ggplot2)
lapply(file.path('~/Desktop/Re-do Poisson/Poisson/functions_Poisson', dir('functions_Poisson')), source)
#setwd("~/Desktop/Re-do Poisson/Poisson")
set.seed( 0 )

n     <- 10              # sample size (n)
mu    <- c(1, 5, 20, 50) # different means/variances (mu)
alpha <- seq(0, 5, 0.25) # P2: alpha
beta  <- seq(0, 5, 0.25) # P2: beta

WS <- data.frame(mu = NA, alpha = NA, beta = NA, stein = NA, vallender = NA)
a  <- 1
for (i in 1:length(mu)) {
  for (j in 1:length(alpha)) {
    for (k in 1:length(beta)) {
      
      # replicate process S times and average over S replicates
      S <- 1000
      d <- matrix(NA, nrow = S, ncol = 2)
      for (s in 1:S) {
        x <- rpois(n, mu[i])
        X <- sum(x)
        d[s,1] <- poisson.GammaGamma.stein(n, X, 2.5, 2.5, alpha[j], beta[k])
        d[s,2] <- poisson.GammaGamma.vallender(n, X, 2.5, 2.5, alpha[j], beta[k])$DW
      }
      WS[a,] <- c(mu[i], alpha[j], beta[k], apply(d, 2, mean))
      a      <- a + 1
    }
  }
}

# save data
save(WS, file = 'data_figure1.rda')

## Figure 1
g1 <- ggplot(WS) +
  geom_raster(aes(x = alpha, y = beta, fill = vallender - stein)) +
  guides(fill = guide_legend(title = 'Difference')) +
  facet_grid(. ~ mu) +
  labs(x = expression(alpha[2]), y = expression(beta[2])) + 
  geom_vline(xintercept = 2.5, lty = 2) +
  geom_hline(yintercept = 2.5, lty = 2)
print(g1)
ggsave('figure1.png', g1, dpi = 300, device = 'png')

