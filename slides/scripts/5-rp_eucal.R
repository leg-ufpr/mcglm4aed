#-----------------------------------------------------------------------
# Packages.

rm(list = ls())
library(lattice)
library(latticeExtra)
library(car)
library(corrplot)
library(candisc)
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
#

# Apply linearHypothesis to all pairwise contrasts. Correct the pvalue
# with bonferroni or another method of corrections. Focus on the
# construction of the L matrix (means) and the K matrix (contrasts).

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
