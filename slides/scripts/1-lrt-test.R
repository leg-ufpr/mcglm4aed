#-----------------------------------------------------------------------
# Packages.

rm(list = ls())

library(mvtnorm)
library(latticeExtra)

#=======================================================================
# Wilks' lambda - likelihood ratio test.

#-----------------------------------------------------------------------
# Creating a function to run the simulation study.

# Simulation settings.
m <- 3
rept <- 4
k <- 4

simul <- function(m, k, rept, N = 1000) {
    trt <- gl(4, rept)
    mu <- rep(0, m)
    n <- length(trt)
    Sigma <- 0.5 * (diag(m) + matrix(1, m, m))
    # Sample of the statistic.
    smp <- replicate(N, {
        Y <- rmvnorm(n = length(trt),
                     mean = mu,
                     sigma = Sigma)
        # Sigma of the null model.
        m0 <- lm(Y ~ 1)
        S0 <- var(residuals(m0)) * (n - 1)/n
        # Sigma of the full model.
        m1 <- lm(Y ~ trt)
        S1 <- var(residuals(m1)) * (n - 1)/n
        # Wilks' statistic.
        stat <- -n * log(det(S1)/det(S0))
        return(stat)
    })
    return(smp)
}

# A sample of the Wilks' lambda statistic.
smp <- simul(m = m, rept = rept, k = k)

# Simulated distribution vs theoretical.
plot(ecdf(smp))
curve(pchisq(x, df = m * (k - 1)), add = TRUE, col = 2)

#-----------------------------------------------------------------------
# Design settings for a small simulation study.

# k:    number of treatment levels (fixed);
# N:    number of simulations (fixed);
# m:    number of responses (varying);
# rept: repetition per treatment level (varying).
k <- 4
N <- 1000
ds <- expand.grid(m = c(2, 5, 8),
                  rept = c(3, 5, 10, 20))

# Simulation for each design setting case.
res <- by(data = ds,
          INDICES = seq_len(nrow(ds)),
          FUN = function(case) {
              smp <- with(case,
                          simul(m = m,
                                rept = rept,
                                k = k,
                                N = N))
              data.frame(m = case$m,
                         rept = case$rept,
                         smp = smp)
          })

# Stack to one data frame.
res <- do.call(rbind, res)
res$n <- res$rept * k

# combineLimits(
    useOuterStrips(
        ecdfplot(~smp | factor(m) + factor(n),
                 as.table = TRUE,
                 # scales = list(x = "free"),
                 xlab = "Wilks' lambda",
                 ylab = "Empirical CDF",
                 data = res),
        strip = strip.custom(
            var.name = "Responses",
            strip.names = TRUE,
            sep = " = "),
        strip.left = strip.custom(
            var.name = "Subjects",
            strip.names = TRUE,
            sep = " = ")) +
    layer({
        panel.curve(pchisq(x, df = ds$m[which.packet()[1]] * k),
                    lty = 2)
    })
# )

# Show the F approximation to the Wilks' lambda test.
#-----------------------------------------------------------------------
