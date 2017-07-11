# Example 1 - MANOVA type test for McGLMs ----------------------------
# Author: Wagner Hugo Bonat LEG/UFPR ---------------------------------
# Date: 14/06/2017 ---------------------------------------------------
rm(list=ls())

# Loading packages ---------------------------------------------------
require(Matrix)
require(mcglm)
require(mvtnorm)
require(bindata)

# Loading extra functions --------------------------------------------
source("functions.R")

# Example 1 ----------------------------------------------------------

# Loading data set ---------------------------------------------------
y1 <- c(463, 438, 494, 496, 448, 603, 596, 616, 633, 608, 471, 481, 449,
        443, 456)/100

y2 <- c(950, 890, 1010, 1230, 940, 1080, 1050, 1080, 1190, 1080, 960,
        930, 870, 820, 910)/1000

trt <- gl(n = 3, k = 5, labels = c("Test", "TurFer", "TurNat"))
da <- data.frame(trt, y1, y2)

# Fitting McGLMs -----------------------------------------------------
ex1_lm <- manova(cbind(y1, y2) ~ trt, data = da)
ex1_mc <- mcglm(c(y1~ trt, y2 ~ trt), list(mc_id(da), mc_id(da)), 
              data = da, control_algorithm = list(correct = FALSE))

summary(ex1_lm, test = "Hotelling-Lawley")
manova.mcglm(ex1_mc)


# Example 2: Iris data set -------------------------------------------
# Standard MANOVA fit
ex2_lm <- manova(cbind(Sepal.Length, Sepal.Width, Petal.Length, Petal.Width) ~ Species, 
               data = iris)

# McGLM fit
Z0 <- mc_id(iris)
ex2_mc <- mcglm(c(Sepal.Length ~ Species, Sepal.Width ~ Species,
                  Petal.Length ~ Species, Petal.Width ~ Species),
                list(Z0,Z0,Z0,Z0), data = iris, 
                control_algorithm = list(correct = FALSE))
summary(ex2_lm, test = "Hotelling-Lawley")
manova.mcglm(ex2_mc)

# Simulation study ---------------------------------------------------
# Gaussian case
Sigma0 <- matrix(c(1,0.8,0.8,1),2,2)
beta0 <- 10
beta1 <- 20
n = 1000
n_trat <- 50
n_simul = 100

result <- matrix(NA, ncol = 4, nrow = n_simul)
for(i in 1:n_simul) {
  data_simul <- rmvnorm(n = n, mean = c(beta0,beta1), sigma = Sigma0)
  trt <- rep(1:n_trat, each = n/n_trat)
  data_fit <- data.frame(data_simul, as.factor(trt))
  names(data_fit) <- c("y1","y2","trt")
  Z0 <- mc_id(data_fit)
  fit_temp <- mcglm(c(y1 ~ trt, y2 ~ trt), list(Z0,Z0), data = data_fit, 
                    control_algorithm = list(correct = FALSE))
  result[i,] <- as.numeric(manova.mcglm(fit_temp)[2,2:5])
  print(i)
}
result[,1]
plot(ecdf(result[,3]))
curve(pchisq(x, df = 98), 0, 400, add = TRUE, col = "red", lty = 2, lwd = 2)


# Binomial case
Sigma0 <- matrix(c(1,0.5,0.5,1),2,2)
beta0 <- 0.5
beta1 <- 0.25
n = 2000
n_trat <- 10
n_simul = 100

result <- matrix(NA, ncol = 4, nrow = n_simul)
for(i in 1:n_simul) {
  data_simul <- rmvbin(n, margprob = c(beta0,beta1), 
                       bincorr = Sigma0)
  trt <- rep(1:n_trat, each = n/n_trat)
  data_fit <- data.frame(data_simul, as.factor(trt))
  names(data_fit) <- c("y1","y2","trt")
  Z0 <- mc_id(data_fit)
  fit_temp <- mcglm(c(y1 ~ trt, y2 ~ trt), list(Z0,Z0), data = data_fit,
                    link = c("logit","logit"), 
                    variance = c("binomialP","binomialP"),
                    control_algorithm = list(correct = FALSE))
  result[i,] <- as.numeric(manova.mcglm(fit_temp)[2,2:5])
  print(i)
}
result[,1]
plot(ecdf(result[,3]))
curve(pchisq(x, df = 18), 0, 40, add = TRUE, col = "red", lty = 2, lwd = 2)

