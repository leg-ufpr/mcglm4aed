#-----------------------------------------------------------------------
# Packages.

library(latticeExtra)
library(car)
library(reshape2)
source("../config/_setup.R")

#-----------------------------------------------------------------------
# Getting the dataset.

# Online documentation of the EACS::teca_qui dataset.
u <- "http://leg.ufpr.br/~walmes/pacotes/EACS/reference/teca_qui.html"
browseURL(u)

csv <- "https://raw.githubusercontent.com/walmes/EACS/master/data-raw/teca_qui.csv"
teca <- read.csv2(file = csv, dec = ".")
str(teca)

teca$cam <- factor(teca$cam, labels = 1:3)

#-----------------------------------------------------------------------

# The cations.
v <- c("k", "ca", "mg")
summary(teca[v])

# Only the cations.
densityplot(~k + ca + I(mg + 0.1),
            outer = TRUE,
            groups = cam,
            scales = list(relation = "free",
                          x = list(log = 10)),
            as.table = TRUE,
            data = teca)

scatterplotMatrix(~log(k) + log(ca) + log(mg + 0.1) | cam,
                  data = teca,
                  gap = 0,
                  smooth = FALSE,
                  reg.line = FALSE,
                  ellipse = TRUE,
                  by.groups = TRUE,
                  diagonal = "qqplot")

teca <- transform(teca,
                  lk = log(k),
                  lca = log(ca),
                  lmg = log(mg + 0.1))

# Transformed variables.
v <- tail(names(teca), n = 3)

#-----------------------------------------------------------------------
# Repeated measures design.

# Outer factor: none.
# Inner factor: soil layer (`cam`)
# Responses: 3 cations x 3 layers = 9 conditions.

# Long format.
tecal <- melt(data = teca[c("loc", "cam", v)],
              measure.vars = v,
              variable.name = "res")
str(tecal)

bwplot(value ~ cam | res,
       pch = "|",
       data = tecal) +
    layer(panel.xyplot(x = x,
                       y = y,
                       jitter.x = TRUE,
                       type = c("p", "a")))

# Combine 3 responses x 3 layers = 9 response variables.
tecal$res.cam <- with(tecal, paste(res, cam, sep = "."))

# Wide format.
tecaw <- dcast(data = tecal,
               formula = loc ~ res.cam,
               value = "value")
str(tecaw)

#-----------------------------------------------------------------------
# Repeated measures analysis.

dput(names(tecaw)[-1])

m0 <- lm(as.matrix(tecaw[, 2:10]) ~ 1)
m0

summary(m0)
summary.aov(m0)
anova(m0)

r <- residuals(m0)

# Checking the models assumptions on the residuals.
scatterplotMatrix(r,
                  gap = 0,
                  smooth = FALSE,
                  reg.line = FALSE,
                  ellipse = TRUE,
                  diagonal = "qqplot")

# Inner factors data design.
idata <- expand.grid(res = levels(tecal$res),
                     cam = levels(tecal$cam))

Anova(m0,
      idata = idata,
      idesign = ~res * cam)

# TODO How interpret these results? The difference in means among soil
# layers depends on the cation?

# How the calculations is done? How construct the L matrix to test the
# hypothesis?

an0 <- lm(value ~ res * cam, data = tecal)
anova(an0)

#-----------------------------------------------------------------------
install_github("wbonat/mcglm", ref = "devel")
require(mcglm)

# Subset
data <- teca[,c(1,2,16,17,18)]
data = data[order(data$cam),]
data = data[order(data$loc),]
data$cam <- as.factor(data$cam)

# Linear predictors
form.lk <- lk ~ cam
form.lca <- lca ~ cam
form.lmg <- lmg ~ cam

# Matrix linear predictors
Z0 <- mc_id(data)

# Unstructured model for covariance between cam
Z_ns <- mc_ns(data, id = "loc")

# Moving average first order
Z_ma1 <- mc_ma(id = "loc", time = "cam", data = data, order = 1)

# Distance based
Z_dist <- mc_dist(id = "loc", time = "cam", data = data, 
                  method = "euclidean")

# Random walk
Z_rw <- mc_rw(id = "loc", time = "cam", data = data, order = 1, 
              proper = TRUE)

# Fitting
# Standard MANOVA
fit1 <- mcglm(linear_pred = c(form.lk, form.lca, form.lmg), 
              matrix_pred = list(Z0,Z0,Z0), data = data)

# MANOVA for all response
Z00 <- mc_id(tecaw)
linear_pred <- c(lca.1 ~ 1, lca.2 ~ 1, lca.3 ~ 1, 
                 lk.1 ~ 1, lk.2 ~ 1, lk.3 ~ 1,
                 lmg.1 ~ 1, lmg.2 ~ 1, lmg.3 ~ 1)

fit1.1 <- mcglm(linear_pred = linear_pred, 
                matrix_pred = list(Z00,Z00,Z00,Z00,Z00,Z00,Z00,Z00,Z00),
                data = tecaw)

# MANOVA + repeated measures using unstructured matrix
fit2 <- mcglm(linear_pred = c(form.lk, form.lca, form.lmg),
              matrix_pred = list(c(Z0,Z_ns), c(Z0,Z_ns), c(Z0,Z_ns)),
              control_algorithm = list(tunning = 0.9, max_iter = 100),
              data = data)

# MANOVA + repeated measures using moving average first order
fit3 <- mcglm(linear_pred = c(form.lk, form.lca, form.lmg),
              matrix_pred = list(c(Z0,Z_ma1), c(Z0,Z_ma1), c(Z0,Z_ma1)),
              control_algorithm = list(tunning = 0.8),
              data = data)

# MANOVA + repeated measures using distance based
fit4 <- mcglm(linear_pred = c(form.lk, form.lca, form.lmg),
              matrix_pred = list(c(Z0,Z_dist), c(Z0,Z_dist), c(Z0,Z_dist)),
              control_algorithm = list(tunning = 0.8),
              data = data)

# MANOVA + repeated measures using distance based + expm covariance link function
fit5 <- mcglm(linear_pred = c(form.lk, form.lca, form.lmg),
              matrix_pred = list(c(Z0,Z_dist), c(Z0,Z_dist), c(Z0,Z_dist)),
              covariance = c("expm","expm","expm"),
              control_algorithm = list(tunning = 0.8),
              data = data)

# MANOVA + repeated measures using distance based + expm covariance link function
fit6 <- mcglm(linear_pred = c(form.lk, form.lca, form.lmg),
              matrix_pred = list(c(Z_rw), c(Z_rw), c(Z_rw)),
              covariance = c("inverse","inverse","inverse"),
              control_algorithm = list(tunning = 0.5, max_iter = 100),
              data = data)

# Comparing models fit
rbind(gof(fit1), gof(fit1.1), gof(fit2), gof(fit3), 
      gof(fit4), gof(fit5), gof(fit6))

# Multivariate hypotheses test
source("../../review/functions.R")
manova.mcglm(fit1)
manova.mcglm(fit2)

# Correlation between outcomes
summary(fit2, print = "Correlation")

# Correlation between cam
# Response lk
COR_lk <- matrix(NA, 3, 3)
COR_lk[lower.tri(COR_lk)] <- fit2$Covariance[5:7]
COR_lk[upper.tri(COR_lk)] <- fit2$Covariance[5:7]
diag(COR_lk) <- fit2$Covariance[4]
cov2cor(COR_lk)

# Response lca
COR_lca <- matrix(NA, 3, 3)
COR_lca[lower.tri(COR_lca)] <- fit2$Covariance[9:11]
COR_lca[upper.tri(COR_lca)] <- fit2$Covariance[9:11]
diag(COR_lca) <- fit2$Covariance[8]
cov2cor(COR_lca)

# Response lmg
COR_lmg <- matrix(NA, 3, 3)
COR_lmg[lower.tri(COR_lmg)] <- fit2$Covariance[13:15]
COR_lmg[upper.tri(COR_lmg)] <- fit2$Covariance[13:15]
diag(COR_lmg) <- fit2$Covariance[12]
cov2cor(COR_lmg)


