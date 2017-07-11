# Course: Multivariate covariance generalized linear models for ------
# analysis of experimental data --------------------------------------
# Author: Wagner Hugo Bonat and Walmes Marques Zeviani LEG/UFPR ------
# Date: 12/06/2017 ---------------------------------------------------
rm(list=ls())

# Loading extra packages ---------------------------------------------
require(devtools)
require(Matrix)
require(mcglm)
require(lattice)
require(car)
require(mvtnorm)
# Installing package with the data sets from github ------------------ 
install_github("walmes/EACS")
require(EACS)

# Loading extra functions
source("functions.R")

# Loading data set - teca_qui ----------------------------------------
data(teca_qui)

# Four response variables --------------------------------------------
scatterplotMatrix(~ p + k + ca + mg | cam, diagonal = "density",
                  data = teca_qui, smooth = FALSE, reg.line = FALSE, 
                  ellipse = TRUE, by.groups = TRUE, 
                  legend.pos = "topright")

# Standard analysis using manova -------------------------------------
mod1 <- manova(cbind(p,k,ca,mg) ~ cam, data = teca_qui)
coef(mod1)
summary(mod1, test = "Hotelling-Lawley")

# manova by mcglm package --------------------------------------------
Z0 <- mc_id(teca_qui)
# MLE
mod2 <- mcglm(c(p ~ cam, k ~ cam, ca ~ cam, mg ~ cam),
              matrix_pred = list(Z0,Z0,Z0,Z0), 
              control_algorithm = list(correct = FALSE), 
              data = teca_qui)
manova.mcglm(mod2)

# RMLE
mod3 <- mcglm(c(p ~ cam, k ~ cam, ca ~ cam, mg ~ cam),
              matrix_pred = list(Z0,Z0,Z0,Z0), 
              control_algorithm = list(correct = TRUE), 
              data = teca_qui)
manova.mcglm(mod3)

# Comparing log-likelihood values ------------------------------------
logLik.mlm(mod1)
plogLik(mod2)
plogLik(mod3)

# How to take into account the repeated measures structure -----------





