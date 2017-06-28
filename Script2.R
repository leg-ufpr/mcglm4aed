# Course: Multivariate covariance generalized linear models for ------
# analysis of experimental data --------------------------------------
# Author: Wagner Hugo Bonat and Walmes Marques Zeviani LEG/UFPR ------
# Date: 12/06/2017 ---------------------------------------------------
rm(list = ls())

# Loading extra packages ---------------------------------------------
library(lattice)
library(latticeExtra)
library(ellipse)
library(mcglm)

# Dataset ------------------------------------------------------------

y1 <- c(463, 438, 494, 496, 448, 603, 596, 616, 633, 608, 471, 481, 449,
        443, 456)/100

y2 <- c(950, 890, 1010, 1230, 940, 1080, 1050, 1080, 1190, 1080, 960,
        930, 870, 820, 910)/1000

trt <- gl(n = 3, k = 5, labels = c("Test", "TurFer", "TurNat"))
da <- data.frame(trt, y1, y2)
da


# Exploratory analysis -----------------------------------------------
xyplot(y2 ~ y1, groups = trt, data = da, xlab = "y1", ylab = "y2",
       auto.key = list(corner = c(0.05, 0.95))) + 
  glayer({ panel.rug(x = x, y = y, col = col.line)
    panel.abline(v = mean(x), h = mean(y), col = col.line, lty = 3)
    xy <- cbind(x, y)
    v <- var(xy)
    m <- colMeans(xy)
    ell <- ellipse::ellipse(v, centre = m, level = 0.8) 
    panel.lines(ell[, 1], ell[, 2], col = col.line, lty = 2)
  })

# Multivariate (gaussian) linear model - manova approach -------------
m0 <- lm(cbind(y1, y2) ~ trt, data = da)

# mcglm approach -----------------------------------------------------
mc0 <- mcglm(linear_pred = c(y1 ~ trt, y2 ~ trt),
             matrix_pred = list(mc_id(da), mc_id(da)), data = da,
             control_algorithm = list(correct = FALSE))

# Estimates
coef(m0)
matrix(coef(mc0, type = "beta")[,1], 3,2)

# Variance-covariance matrix
round(vcov(m0),4)
round(vcov(mc0)[1:6,1:6],4)

# A Wald F test for treatment effect in each response ----------------
V_lm <- vcov(m0)
B_lm <- matrix(coef(m0), dimnames = list(rownames(V_lm)))
L <- rbind(c(0, 1, 0), c(0, 0, 1))

# For y1.
U1 <- kronecker(rbind(1:0), L)
UB1_lm <- U1 %*% B_lm
UVU1_lm <- U1 %*% V_lm %*% t(U1)
(t(UB1_lm) %*% solve(UVU1_lm) %*% UB1_lm)/nrow(U1)
anova(aov(y1 ~ trt,data = da))

# For y2.
U2 <- kronecker(rbind(0:1), L)
UB2_lm <- U2 %*% B_lm
UVU2_lm <- U2 %*% V_lm %*% t(U2)
(t(UB2_lm) %*% solve(UVU2_lm) %*% UB2_lm)/nrow(UB2_lm)
anova(aov(y2 ~ trt,data = da))

# Joint test using lm results ----------------------------------------
U_joint <- rbind(U1, U2)
U_j_lm <- U_joint %*% B_lm
UVUj_lm <- U_joint %*% V_lm %*% t(U_joint)
(t(U_j_lm) %*% solve(UVUj_lm) %*% U_j_lm)/nrow(U_j_lm)

summary(manova(cbind(y1, y2) ~ trt, data = da), test = "Pillai")
summary(manova(cbind(y1, y2) ~ trt, data = da), test = "Wilks")
summary(manova(cbind(y1, y2) ~ trt, data = da), test = "Roy")
summary(manova(cbind(y1, y2) ~ trt, data = da), test = "Hotelling-Lawley")

# NOTE: these F statistics are the same from univariate anova, as it
# should be.

# Doing the same for mcglm -------------------------------------------
V_mc <- vcov(mc0)[1:6,1:6]
B_mc <- matrix(coef(mc0)[1:6,1], dimnames = list(rownames(V_mc)))
L <- rbind(c(0, 1, 0), c(0, 0, 1))

# For y1.
U1 <- kronecker(rbind(1:0), L)
UB1_mc <- U1 %*% B_mc
UVU1_mc <- U1 %*% V_mc %*% t(U1)
(t(UB1_mc) %*% solve(UVU1_mc) %*% UB1_mc)/nrow(U1)
anova(aov(y1 ~ trt,data = da))

# For y2.
U2 <- kronecker(rbind(0:1), L)
UB2_mc <- U2 %*% B_mc
UVU2_mc <- U2 %*% V_mc %*% t(U2)
(t(UB2_mc) %*% solve(UVU2_mc) %*% UB2_mc)/nrow(UB2_mc)
anova(aov(y2 ~ trt,data = da))

# Joint test using mcglm results -------------------------------------
U_joint <- rbind(U1, U2)
U_j_mc <- U_joint %*% B_mc
UVUj_mc <- U_joint %*% V_mc %*% t(U_joint)
(t(U_j_mc) %*% solve(UVUj_mc) %*% U_j_mc)/nrow(U_j_mc)


# Checking results from article --------------------------------------
X <- model.matrix(~ trt,data = da)
Y <- with(da, cbind(y1, y2))
B <- coef(m0)
S <- t(Y - X%*%B)%*%(Y-X%*%B)
FF <- L
G <- diag(1, ncol = 2, nrow = 2)
E <- matrix(0, nrow = 2, ncol = 2)

H0 <- FF%*%B%*%G
E
# Likelihood ratio test
W0 <- S + S%*%G%*%solve(t(G)%*%S%*%G)%*%t(FF%*%B%*%G-E)%*%solve(FF%*%solve(t(X)%*%X)%*%t(FF))%*%(FF%*%B%*%G-E)%*%solve(t(G)%*%S%*%G)%*%t(G)%*%S
W0
# Wilks test
delta <- det(S)/det(W0)
delta
summary(manova(cbind(y1, y2) ~ trt, data = da), test = "Wilks") # OK

# Or equivalent to
m00 <- lm(cbind(y1,y2) ~ 1)
X0 <- model.matrix(~ 1, data = da)
W00 <- t(Y - X0%*%coef(m00))%*%(Y - X0%*%coef(m00))
W00

# Another alternative
A <- t(FF%*%B%*%G - E)%*%solve(FF%*%solve(t(X)%*%X)%*%t(FF))%*%(FF%*%B%*%G - E)
B <- t(G)%*%S%*%G
det(B)/det(A+B)

# Distribution of test statistics under H0 ---------------------------
mod_nul <- lm(cbind(y1,y2) ~ 1, data = da)
mod_alt <- lm(cbind(y1,y2) ~ trt, data = da)
W0 <- t(Y - X0%*%coef(mod_nul))%*%(Y - X0%*%coef(mod_nul))

lambda <- c()
for(i in 1:10000) {
Y_simul_H0 <- rmvnorm(n = 15, mean = coef(mod_nul), sigma = W0)
mod_nul <- lm(Y_simul_H0 ~ 1)
mod_alt <- lm(Y_simul_H0 ~ da$trt)
B <- coef(mod_alt)
S <- t(Y - X%*%B)%*%(Y-X%*%B)
A <- t(FF%*%B%*%G - E)%*%solve(FF%*%solve(t(X)%*%X)%*%t(FF))%*%(FF%*%B%*%G - E)
B <- t(G)%*%S%*%G
lambda[i] <- det(B)/det(A+B)
}

hist(-2*log(lambda), prob = TRUE)
curve(dchisq(x, df = 4), 0, 15, add = TRUE)
-2*log(delta)
pchisq(-2*log(delta), df = 4, lower.tail = FALSE)

