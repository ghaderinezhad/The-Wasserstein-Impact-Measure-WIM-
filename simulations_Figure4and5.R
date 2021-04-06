#' Poisson model: simulations for Figure 4 and 5
#' Wasserstein distance vs. Neutrality and MOPESS
#' 
# (1) Flat (alpha = 1, beta = 0) ==> reference prior
# (2) exponential prior (alpha = 1, beta)
# (3) gamma prior (beta, beta)


rm(list = ls())
library(ggplot2)
lapply(file.path('functions_Poisson', dir('functions_Poisson')), source)

set.seed( 0 )

n  <- c(10, 50, 100, 200)  # different sample sizes (n)
mu <- c(1, 10, 20)         # different means (mu)
a1 <- 1                    # P1: alpha
b1 <- 0                    # P1: beta
a2 <- 1                    # P2: alpha
b2 <- seq(0, 5, 0.5)       # P2: beta

WS <- data.frame(N = NA, mu = NA, beta = NA, prior = NA, 
                 stein = NA, neutrality = NA, MOPESS = NA)
a  <- 1
for (i in 1:length(n)) {
  for (j in 1:length(mu)) {
    for (k in 1:length(b2)) {
      
      # replicate process S times and average over S replicates
      S  <- 10
      dE <- dG <- matrix(NA, nrow = S, ncol = 3)
      for (s in 1:S) {
        x       <- rpois(n[i], mu[j])
        X       <- sum(x)
        dE[s,1] <- poisson.GammaGamma.stein(n[i], X, a1, b1, a2, b2[k])
        dG[s,1] <- poisson.GammaGamma.stein(n[i], X, a1, b1, b2[k], b2[k])
        dE[s,2] <- poisson.GammaGamma.neutrality(n[i], X, a1, b1, a2, b2[k])$N2
        dG[s,2] <- poisson.GammaGamma.neutrality(n[i], X, a1, b1, b2[k], b2[k])$N2
        dE[s,3] <- poisson.GammaGamma.mopess(x, 20, a1, b1, a2, b2[k])$MOPESS
        dG[s,3] <- poisson.GammaGamma.mopess(x, 20, a1, b1, b2[k], b2[k])$MOPESS
      }
      WS[a,] <- c(n[i], mu[j], b2[k], 'Exponential(b)', apply(dE, 2, mean))
      a      <- a + 1
      WS[a,] <- c(n[i], mu[j], b2[k], 'Gamma(b,b)', apply(dG, 2, mean))
      a      <- a + 1
    }
  }
}

WS$N          <- as.numeric(WS$N)
WS$mu         <- as.numeric(WS$mu)
WS$beta       <- as.numeric(WS$beta)
WS$stein      <- as.numeric(WS$stein)
WS$neutrality <- as.numeric(WS$neutrality)
WS$MOPESS     <- as.numeric(WS$MOPESS)

# save data
save(WS, file = 'data_figure4and5.rda')


## Figure 4
g1 <- ggplot(WS) +
  geom_line(aes(x = stein, y = neutrality, color = beta), lty = 1, lwd = 1,
            data = WS[WS$prior == 'Exponential(b)',]) +
  geom_point(aes(x = stein, y = neutrality, color = beta),
             data = WS[WS$prior == 'Gamma(b,b)',]) +
  facet_grid(mu ~ N, scales = 'free') + 
  guides(color = guide_legend(title = expression(beta))) +
  labs(x = 'WIM', y = 'Neutrality')
print(g1)
ggsave('figure4.png', g1, dpi = 300, device = 'png')


## Figure 5
g2 <- ggplot(WS) +
  geom_line(aes(x = stein, y = MOPESS, color = beta), lty = 1, lwd = 1,
            data = WS[WS$prior == 'Exponential(b)',]) +
  geom_point(aes(x = stein, y = MOPESS, color = beta),
             data = WS[WS$prior == 'Gamma(b,b)',]) +
  facet_grid(mu ~ N, scales = 'free') + 
  guides(color = guide_legend(title = expression(beta))) +
  labs(x = 'WIM', y = 'MOPESS')
print(g2)
ggsave('figure5.png', g2, dpi = 300, device = 'png')


