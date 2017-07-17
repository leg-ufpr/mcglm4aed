#-----------------------------------------------------------------------
# Packages.

rm(list = ls())
library(lattice)
library(latticeExtra)
library(car)
library(candisc)

#=======================================================================
# Warming up with the `iris` dataset.

scatterplotMatrix(~Sepal.Length + Sepal.Width +
                      Petal.Length + Petal.Width | Species,
                  data = iris,
                  smooth = FALSE,
                  reg.line = FALSE,
                  ellipse = TRUE,
                  by.groups = TRUE,
                  gap = 0,
                  diagonal = "density")

#-----------------------------------------------------------------------
# Petal.Length & Petal.Width.

scatterplotMatrix(~Petal.Length + Petal.Width | Species,
                  data = iris,
                  smooth = FALSE,
                  reg.line = FALSE,
                  ellipse = TRUE,
                  by.groups = TRUE,
                  gap = 0,
                  diagonal = "density")

# Iris full model.
m1 <- lm(cbind(Petal.Length,
               Petal.Width) ~ Species,
         data = iris)

# Iris null model.
m0 <- update(m1, . ~ 1)

# MANOVA with Pillai test.
anova(m1)

# Error SSP of the full and null models.
S_full <- crossprod(residuals(m1))
S_null <- crossprod(residuals(m0))

# Extra SSP due the hypothesis.
S_extr <- S_null - S_full

# Eigen decomposition -> canonical analysis.
ei <- eigen(solve(S_full, S_extr))
ei

# Cumulated proportion.
cumsum(ei$values)/sum(ei$values)

# Scores (canonical variables).
scores <- as.matrix(m1$model[1]) %*% ei$vectors[, 1:2]

# Plot of the scores.
xyplot(scores[, 2] ~ scores[, 1],
       groups = iris$Species,
       grid = TRUE,
       auto.key = TRUE,
       aspect = "iso") +
    layer(panel.ellipse(...))

#--------------------------------------------
# Using the `candisc` package.

c1 <- candisc(m1, term = "Species")
str(c1)

summary(c1)

ei$vectors[, 1:2]
c1$coeffs.raw
c1$coeffs.raw/ei$vectors[, 1:2]

scatterplotMatrix(~Can1 + Can2| Species,
                  data = c1$scores,
                  smooth = FALSE,
                  reg.line = FALSE,
                  ellipse = TRUE,
                  by.groups = TRUE,
                  gap = 0,
                  diagonal = "density")


# Scores: canonical variates.
xy_can <-
    xyplot(Can2 ~ Can1,
           groups = Species,
           data = c1$scores,
           xlab = NULL,
           ylab = NULL,
           aspect = "iso",
           auto.key = list(columns = 3),
           grid = TRUE) +
    layer(panel.ellipse(...)) +
    layer(panel.arrows(0, 0,
                       c1$coeffs.raw[, 1],
                       c1$coeffs.raw[, 2],
                       length = 0.1))

xy_ori <-
    xyplot(Petal.Length ~ Petal.Width,
           groups = Species,
           aspect = "iso",
           data = iris,
           grid = TRUE) +
    layer(panel.ellipse(...))

c(xy_can, xy_ori, layout = c(1, 2))

#--------------------------------------------
# Univariate analysis of each canical variate.

an0 <- lm(scores ~ Species, data = iris)
summary(an0)

#-----------------------------------------------------------------------
# Sepal.Length & Sepal.Width.

scatterplotMatrix(~Sepal.Length + Sepal.Width | Species,
                  data = iris,
                  smooth = FALSE,
                  reg.line = FALSE,
                  ellipse = TRUE,
                  by.groups = TRUE,
                  gap = 0,
                  diagonal = "density")

# Iris full model.
m1 <- lm(cbind(Sepal.Length,
               Sepal.Width) ~ Species,
         data = iris)

# Iris null model.
m0 <- update(m1, . ~ 1)

# MANOVA with Pillai test.
anova(m1)

# Error SSP of the full and null models.
S_full <- crossprod(residuals(m1))
S_null <- crossprod(residuals(m0))

# Extra SSP due the hypothesis.
S_extr <- S_null - S_full

# Eigen decomposition -> canonical analysis.
ei <- eigen(solve(S_full, S_extr))
ei

# Cumulated proportion.
cumsum(ei$values)/sum(ei$values)

# Scores (canonical variables).
scores <- as.matrix(m1$model[1]) %*% ei$vectors[, 1:2]

# Plot of the scores.
xyplot(scores[, 2] ~ scores[, 1],
       groups = iris$Species,
       grid = TRUE,
       auto.key = TRUE,
       aspect = "iso") +
    layer(panel.ellipse(...))

#--------------------------------------------
# Using the `candisc` package.

c1 <- candisc(m1, term = "Species")
str(c1)

summary(c1)

ei$vectors[, 1:2]
c1$coeffs.raw
c1$coeffs.raw/ei$vectors[, 1:2]

scatterplotMatrix(~Can1 + Can2| Species,
                  data = c1$scores,
                  smooth = FALSE,
                  reg.line = FALSE,
                  ellipse = TRUE,
                  by.groups = TRUE,
                  gap = 0,
                  diagonal = "density")


# Scores: canonical variates.
xy_can <-
    xyplot(Can2 ~ Can1,
           groups = Species,
           data = c1$scores,
           xlab = NULL,
           ylab = NULL,
           aspect = "iso",
           auto.key = list(columns = 3),
           grid = TRUE) +
    layer(panel.ellipse(...)) +
    layer(panel.arrows(0, 0,
                       c1$coeffs.raw[, 1],
                       c1$coeffs.raw[, 2],
                       length = 0.1))

xy_ori <-
    xyplot(Sepal.Length ~ Sepal.Width,
           groups = Species,
           aspect = "iso",
           data = iris,
           grid = TRUE) +
    layer(panel.ellipse(...))

c(xy_can, xy_ori, layout = c(1, 2))

#--------------------------------------------
# Univariate analysis of each canical variate.

an0 <- lm(scores ~ Species, data = iris)
summary(an0)

#=======================================================================
# Multivariate analysis of the soybean dataset.

data(soybeanwp, package = "wzRfun")
str(soybeanwp)
head(soybeanwp)

# A copy with shorter names.
swp <- soybeanwp
# dput(names(swp))
names(swp) <- c("pts", "wtr", "blk",
                "yield", "w100", "Kconc", "tg", "nip", "nvp")

# Convert the controled numeric factors to categorical factors.
swp <- transform(swp,
                 Pts = factor(pts),
                 Wtr = factor(wtr))

# Creates the proportion of viable pods.
swp$pvp <- with(swp, nvp/(nvp + nip))
str(swp)

#-----------------------------------------------------------------------
# Data visulization.

# The 74 observation is an outlier to yield and tg.

combineLimits(
    useOuterStrips(
        xyplot(yield + tg + w100 + Kconc + pvp ~ pts | Wtr,
               outer = TRUE,
               as.table = TRUE,
               # data = swp[-74, ],
               data = swp,
               pch = 19,
               xlab = "Potassium",
               ylab = "Values",
               scales = list(y = "free"),
               groups = blk))) +
    layer(panel.xyplot(x = x, y = y, type = "a", col = 1))

#-----------------------------------------------------------------------
# Fitting the multivariate Gaussian linear model.

m0 <- lm(cbind(yield, tg, w100, Kconc, pvp) ~ blk + Wtr * Pts,
         data = swp[-74, ])

# MANOVA.
anova(m0)

# ANOVAs.
summary.aov(m0)

# Univariate model summaries.
summary(m0)

# Estimated parameter.
coef(m0)

# Raw residuals.
r <- residuals(m0)
var(r)
cor(r)

# Change the 3rd color of the palette used for the ellipses.
oldpal <- palette()
palette(c("black", "blue", "red"))
scatterplotMatrix(r,
                  gap = 0,
                  smooth = FALSE,
                  reg.line = FALSE,
                  ellipse = TRUE,
                  diagonal = "qqplot")
palette(oldpal)

#-----------------------------------------------------------------------
# Canonical analysis.

cd_WP <- candisc(m0, term = "Wtr:Pts")
cd_WP
cd_WP$coeffs.raw
head(cd_WP$scores)

cd_P <- candisc(m0, term = "Pts")
cd_P
cd_P$coeffs.raw
head(cd_P$scores)

cd_W <- candisc(m0, term = "Wtr")
cd_W
cd_W$coeffs.raw
head(cd_W$scores)

par(mfrow = c(1, 3))
plot(cd_W)
plot(cd_P)
plot(cd_WP)
layout(1)

#=======================================================================

# full <- cbind(yield, tg, w100, Kconc, pvp) ~ blk + Wtr * Pts
# null <- . ~ blk + Wtr + Pts
# data <- swp[-74, ]

gugudada <- function(full, null, data) {
    # Nested models.
    m_full <- lm(full, data)
    m_null <- update(m_full, null)
    # SSP.
    S_full <- crossprod(residuals(m_full))
    S_null <- crossprod(residuals(m_null))
    S_extr <- S_null - S_full
    # Eigen decomposition.
    ei <- eigen(solve(S_full, S_extr))
    print(ei)
    scores <- as.matrix(m_full$model[1]) %*% ei$vectors[, 1:2]
    plot(scores, asp = 1)
    print(anova(update(m_full, scores[, 1] ~ .)))
    print(anova(update(m_full, scores[, 2] ~ .)))
}

# Effect of interaction.
par(mfrow = c(1, 2))
plot(cd_WP)
grid()
gugudada(full = cbind(yield, tg, w100, Kconc, pvp) ~ blk + Wtr * Pts,
         null = . ~ blk + Wtr + Pts,
         data = swp[-74, ])
grid()
layout(1)

# Effect of Wtr.
par(mfrow = c(1, 2))
plot(cd_W)
grid()
gugudada(full = cbind(yield, tg, w100, Kconc, pvp) ~ blk + Pts + Wtr,
         null = . ~ blk + Pts,
         data = swp[-74, ])
grid()
layout(1)

# ATTENTION: This requires more study to fully understand.

#-----------------------------------------------------------------------
