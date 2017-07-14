#-----------------------------------------------------------------------
# Packages.

rm(list = ls())
library(lattice)
library(latticeExtra)
library(car)
library(corrplot)
library(candisc)
library(doBy)
library(multcomp)

#-----------------------------------------------------------------------
# Loading the `EACS::rp_eucal` dataset.

# URL of the dataset documentation (in portuguese).
u <- "http://leg.ufpr.br/~walmes/pacotes/EACS/reference/rp_eucal.html"
browseURL(u)

u <- "https://raw.githubusercontent.com/walmes/EACS/master/data-raw/rp_eucal.csv"
euc <- read.csv2(file = u, dec = ".")
str(euc)

# Layout of the dataset.
xtabs(~manejo + cam, data = euc)

#-----------------------------------------------------------------------
# Data visualization.

useOuterStrips(
    xyplot(umid0 + umid6 + rp + dens ~ cam | manejo,
           outer = TRUE,
           scales = list(y = list(relation = "free")),
           data = euc,
           jitter.x = TRUE,
           type = c("p", "smooth")))

# ATTENTION: The independence among soil layers in the same sample unit
# is a strong assumption. Thus, for the sake of simplicity, only one
# soil layer will be used in the analysis.

combineLimits(
    useOuterStrips(
        xyplot(umid0 + umid6 + rp + dens ~ manejo | factor(cam),
               outer = TRUE,
               scales = list(y = list(relation = "free")),
               data = euc,
               jitter.x = TRUE,
               as.table = TRUE,
               type = c("p", "smooth"))))

xyplot(umid0 + umid6 + rp + dens ~ manejo,
       outer = TRUE,
       layout = c(4, 1),
       scales = list(y = list(relation = "free")),
       data = euc,
       subset = cam == 5,
       jitter.x = TRUE,
       type = c("p", "a"))

splom(~euc[4:7] | factor(cam),
      as.matrix = TRUE,
      groups = manejo,
      layout = c(3, 3),
      as.table = TRUE,
      data = euc)

splom(~euc[4:7],
      as.matrix = TRUE,
      groups = manejo,
      subset = cam == 5,
      data = euc)

#-----------------------------------------------------------------------
# Multivariate analysis.

eucs <- subset(euc, cam == 5)
str(eucs)

# The full model.
m_full <- lm(cbind(umid0, umid6, rp, dens) ~ bloc + manejo,
             data = eucs)

# A little of checking model assumptions.
r <- residuals(m_full)
scatterplotMatrix(r,
                  smooth = FALSE,
                  reg.line = FALSE,
                  ellipse = TRUE,
                  gap = 0,
                  diagonal = "qqplot")

# corrplot(cor(r),
#          method = "ellipse",
#          mar = c(5, 5, 5, 5),
#          type = "lower")

anova(m_full, test = "Pillai")
anova(m_full, test = "Roy")

#-----------------------------------------------------------------------
# Exploring linear hypothesis tests.

L <- LSmatrix(m_full, effect = "manejo")
grid <- attr(L, "grid")

# Estimated means for each response.
L %*% m_full[[1]]

# A set of structured contrasts.
rownames(L) <- grid$manejo
rownames(L)

# List with the planned contrasts.
K <- list("SvsB"       = L[4, ] - colMeans(L[-4, ]),
          "BSLvsBS+BL" = L[3, ] - colMeans(L[1:2, ]),
          "BLvsBS"     = L[1, ] - L[2, ])

# Test the first contrast.
lh <- linearHypothesis(model = m_full,
                       hypothesis.matrix = K[[1]],
                       test = "Pillai")
print(lh, SSP = TRUE)

# Eval the planned contrasts serially.
k <- lapply(K,
            FUN = linearHypothesis,
            model = m_full,
            test = "Pillai")
for (x in k) {
    print(x, SSP = FALSE)
}

# All pairwise contrasts.
A <- wzRfun::apc(L)
A <- by(A, rownames(A), as.matrix)

# Eval the contrasts serially.
a <- lapply(A,
            FUN = linearHypothesis,
            model = m_full,
            test = "Pillai")
sapply(a,
       FUN = function(x) {
           print(x, SSP = FALSE)
           invisible()
       })

#-----------------------------------------------------------------------

# Contrast BSLvsBS+BL in 4 responses.
linearHypothesis(model = m_full,
                 hypothesis.matrix = K[[2]],
                 test = "Pillai",
                 P = diag(4))

# Contrast BSLvsBS+BL in the third response (rp).
linearHypothesis(model = m_full,
                 hypothesis.matrix = K[[2]],
                 test = "Pillai",
                 P = cbind(c(0, 0, 1, 0)))

# Contrast BSLvsBS+BL in the first two responses (umid0 and umid6).
linearHypothesis(model = m_full,
                 hypothesis.matrix = K[[2]],
                 test = "Pillai",
                 P = cbind(c(1, 0, 0, 0),
                           c(0, 1, 0, 0)))

#-----------------------------------------------------------------------
# Which is the variable that most contribute to differentiate the
# groups?

cd <- candisc(m_full, term = "manejo")
cd
summary(cd)

cd$coeffs.raw
cd$coeffs.std
cd$scores

plot(cd)

#-----------------------------------------------------------------------
# Univariate model for the first canonical variate.

splom(~cbind(eucs[4:7], Can1 = cd$scores$Can1),
      pch = 19,
      as.matrix = TRUE,
      auto.key = TRUE,
      groups = manejo,
      data = eucs)

c(xyplot(umid0 + umid6 + rp + dens ~ manejo,
         outer = TRUE,
         scales = list(y = list(relation = "free")),
         data = eucs,
         as.table = TRUE,
         type = c("p", "a")),
  Can1 = xyplot(Can1 ~ manejo,
                data = cd$scores,
                type = c("p", "a")),
  layout = c(NA, 1))

an <- lm(Can1 ~ bloc + manejo,
         data = cd$scores)
anova(an)

g <- glht(an, linfct = mcp(manejo = "Tukey"))
summary(g)
cld(g)

plot(g)

#-----------------------------------------------------------------------
