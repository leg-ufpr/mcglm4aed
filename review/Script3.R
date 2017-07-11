# Multivariate linear hypotheses tests -------------------------------
# Author: Wagner Hugo Bonat LEG/UFPR ---------------------------------
# Date: 13/06/2017 ---------------------------------------------------
rm(list=ls())
# Loading extra packages ---------------------------------------------
library(lattice)
library(latticeExtra)
library(ellipse)
library(mcglm)
require(mvtnorm)

# Dataset ------------------------------------------------------------
y1 <- c(463, 438, 494, 496, 448, 603, 596, 616, 633, 608, 471, 481, 449,
        443, 456)/100

y2 <- c(950, 890, 1010, 1230, 940, 1080, 1050, 1080, 1190, 1080, 960,
        930, 870, 820, 910)/1000

trt <- gl(n = 3, k = 5, labels = c("Test", "TurFer", "TurNat"))
da <- data.frame(trt, y1, y2)

# Likelihood ratio test ----------------------------------------------

# Approach 1 - Computing the likelihood ratio test -------------------
# mcglm
mod_null <- mcglm(c(y1 ~ 1, y2 ~ 1), list(mc_id(da), mc_id(da)), 
                  data = da, control_algorithm = list(correct = FALSE))
mod_alt <- mcglm(c(y1 ~ trt, y2 ~ trt), list(mc_id(da), mc_id(da)), 
                 data = da, control_algorithm = list(correct = FALSE))
plogLik(mod_null)
plogLik(mod_alt)

# lm
mod_lm_null <- lm(cbind(y1, y2) ~ 1,  data = da)
mod_lm_alt <- lm(cbind(y1, y2) ~ trt,  data = da)
logLik(mod_lm_null)
logLik(mod_lm_alt)

Lk <- -2*(logLik(mod_lm_null)[1] - logLik(mod_lm_alt)[1])
pchisq(Lk, df = 4, lower.tail = FALSE)


# Approach 2 - MANOVA + Wilk's test ----------------------------------
summary(manova(cbind(y1,y2) ~ trt, data = da), test = "Wilk")

# Approach 3 - Computing S under the null and the ratio of determinant
TT = 15 # Sample size
M = 2 # Number of responses
K = 3 # Number of beta's for each response
N = 2 # Number of beta under test for each response
FF <-  rbind(c(0, 1, 0), c(0, 0, 1))
G <- diag(1, ncol = M, nrow = N)
E <- matrix(0, nrow = 2, ncol = 2)

B <- matrix(coef(mod_alt)[1:6,1], ncol = M, nrow = K)
X <- model.matrix(~ trt, data = da)

#
S <- t(cbind(da$y1, da$y2) - fitted(mod_alt))%*%(cbind(da$y1, da$y2) - fitted(mod_alt))
W0 <- S + S%*%G%*%solve(t(G)%*%S%*%G)%*%t(FF%*%B%*%G-E)%*%solve(FF%*%solve(t(X)%*%X)%*%t(FF))%*%(FF%*%B%*%G-E)%*%solve(t(G)%*%S%*%G)%*%t(G)%*%S
W00 <- t(cbind(da$y1, da$y2) - fitted(mod_null))%*%(cbind(da$y1, da$y2) - fitted(mod_null))
delta1 <- det(S)/det(W00)

# Approach 4 - Alternative for the above approach --------------------
A <- t(FF%*%B%*%G - E)%*%solve(FF%*%solve(t(X)%*%X)%*%t(FF))%*%(FF%*%B%*%G - E)
B <- t(G)%*%S%*%G
delta2 <- det(B)/det(A+B)
Lk2 <- -TT*log(delta2)

vv <- eigen(A%*%solve(B))$values
delta3 <- prod(1/(1+vv))
c(delta1,delta2,delta3) # Perfect!

pchisq(Lk, df = 4, lower.tail = FALSE) 
pchisq(Lk2, df = 4, lower.tail = FALSE)

# Verifing the asymptotic approximation for the test statistics ------
lambda <- c()
W0 <- matrix(W0, 2,2)
for(i in 1:1000) {
  Y_simul_H0 <- rmvnorm(n = 15, mean = coef(mod_lm_null), sigma = W0)
  Y <- Y_simul_H0
  mod_nul <- lm(Y_simul_H0 ~ 1)
  mod_alt <- lm(Y_simul_H0 ~ rep(da$trt,1))
  B <- coef(mod_alt)
  X <- model.matrix(mod_alt)
  S <- t(Y - X%*%B)%*%(Y-X%*%B)
  A <- t(FF%*%B%*%G - E)%*%solve(FF%*%solve(t(X)%*%X)%*%t(FF))%*%(FF%*%B%*%G - E)
  B <- t(G)%*%S%*%G
  lambda[i] <- det(B)/det(A+B)
}

# Empirical distribution
plot(ecdf(-15*log(lambda)))
# Theoretical distribution
curve(pchisq(x, df = 4), 0, 50, add = TRUE, col = "red", lty = 2, lwd = 2)

# --------------------------------------------------------------------

# Test Wald ----------------------------------------------------------
mod_alt <- mcglm(c(y1 ~ trt, y2 ~ trt), list(mc_id(da), mc_id(da)), 
                 data = da, control_algorithm = list(correct = FALSE))

# Estimated variance-covariance matrix
S <- t(cbind(da$y1, da$y2) - fitted(mod_alt))%*%(cbind(da$y1, da$y2) - fitted(mod_alt))/15

C <- kronecker(G, FF)
beta <- coef(mod_alt)[1:6,1] # Stacked vector of beta's
cc <- rep(0,4)

t(C%*%beta - cc)%*%solve(C%*%(kronecker(S,solve(t(X)%*%X)))%*%t(C))%*%(C%*%beta - cc)
sum(diag(A%*%solve(B)))

kronecker(S,solve(t(X)%*%X))
vcov(mod_alt)[1:6,1:6]
vcov(mod_lm_alt)

t_wald <- as.numeric(t(C%*%beta - cc)%*%solve(C%*%(vcov(mod_alt)[1:6,1:6])%*%t(C))%*%(C%*%beta - cc))/15
pchisq(t_wald, df = 4, lower.tail = FALSE)

summary(manova(cbind(y1,y2) ~ trt, data = da), test = "Hotelling-Lawley")
pf(87.817, 4, 20, lower.tail = FALSE)

