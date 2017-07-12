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
