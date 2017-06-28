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

# Installing package with the data sets from github ------------------ 
install_github("walmes/EACS")
require(EACS)

# Loading data set - gen_teca ----------------------------------------
data(gen_teca)
da <- gen_teca[,c(1,2,3,4,5)]
da <- na.exclude(da)
da$dias <- as.numeric(da$data - min(da$data) + 2)

# Exploratory analysis - response HEIGHT -----------------------------
xyplot(alt ~ dias | gen,
       data = da,
       groups = dose,
       type = "l",
       xlab = "Time after first evaluation (days)",
       ylab = "Height (mm)",
       as.table = TRUE,
       scales = list(y = list(log = FALSE)))

xyplot(dac ~ dias | gen,
       data = da,
       groups = dose,
       type = "l",
       xlab = "Time after first evaluation (days)",
       ylab = "Lap's diameter (mm)",
       as.table = TRUE,
       scales = list(y = list(log = FALSE)))

# Covariates ---------------------------------------------------------
x <- unique(sort(da$dose))
da$dos <- factor(da$dose,
                 levels = x,
                 labels = seq_along(x) - 1)
da$ue <- with(da, interaction(gen, dos, drop = TRUE))
da$dosep <- da$dose^0.3
# Linear predictors --------------------------------------------------
form_alt <- alt ~ dias + I(dias^2) + gen * (dosep + I(dosep^2))
form_dac <- dac ~ dias + I(dias^2) + gen * (dosep + I(dosep^2))

# Univariate models
# Height
fit_alt <- mcglm(c(form_alt), list(mc_id(da)), data = da)
fit_lm <- lm(form_alt, data = da)

# Comparing estimates
par(mfrow = c(2,2))
plot(coef(fit_alt, type = "beta")[,1] ~ coef(fit_lm), 
     ylab = "Estimates-mcglm", xlab = "Estimates-lm")
plot(round(coef(fit_alt, type = "beta")[,1]/coef(fit_lm),4) ~ c(1:89),
     xlab = "beta index", ylab = "Ratio between estimates")
# Comparing standard errors
plot(coef(fit_alt, type = "beta", std.error = TRUE)[,2] ~ sqrt(diag(vcov(fit_lm))),
     ylab = "Std-mcglm", xlab = "Std-lm")
plot(round(coef(fit_alt, type = "beta", std.error = TRUE)[,2]/sqrt(diag(vcov(fit_lm))),4),
     ylab = "Ratio", xlab = "Beta index")

# dac
fit_dac <- mcglm(c(form_dac), list(mc_id(da)), data = da)
fit_lm_dac <- lm(form_dac, data = da)

# Comparing estimates
par(mfrow = c(2,2))
plot(coef(fit_dac, type = "beta")[,1] ~ coef(fit_lm_dac), 
     ylab = "Estimates-mcglm", xlab = "Estimates-lm")
plot(round(coef(fit_dac, type = "beta")[,1]/coef(fit_lm_dac),4) ~ c(1:89),
     xlab = "beta index", ylab = "Ratio between estimates")
# Comparing standard errors
plot(coef(fit_dac, type = "beta", std.error = TRUE)[,2] ~ sqrt(diag(vcov(fit_lm_dac))),
     ylab = "Std-mcglm", xlab = "Std-lm")
plot(round(coef(fit_dac, type = "beta", std.error = TRUE)[,2]/sqrt(diag(vcov(fit_lm_dac))),4),
     ylab = "Ratio", xlab = "Beta index")

# Bivariate models ---------------------------------------------------
fit_joint <- mcglm(c(form_alt, form_dac), list(mc_id(da), mc_id(da)), 
                   data = da)
fit_lm_joint <- lm(cbind(alt, dac) ~ dias + I(dias^2) + gen * (dosep + I(dosep^2)), data = da)


plot(coef(fit_joint, type = "beta", response = 1)[,1]~ coef(fit_lm_joint)[,1])
plot(coef(fit_joint, type = "beta", response = 2)[,1]~ coef(fit_lm_joint)[,2])

plot(coef(fit_joint, type = "beta", response = 1, std.error = TRUE)[,2] ~ summary(fit_lm_joint)[[1]][[4]][,2])
plot(coef(fit_joint, type = "beta", response = 2, std.error = TRUE)[,2] ~ summary(fit_lm_joint)[[2]][[4]][,2])

manova(fit_lm_joint)

