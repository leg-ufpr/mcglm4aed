#-----------------------------------------------------------------------
# Packages.

library(lattice)

#-----------------------------------------------------------------------
# Pimentel Gomes, page 254.

da <- data.frame(trt = gl(n = 3,
                          k = 5,
                          labels = c("Test", "TurFer", "TurNat")),
                 y1 = c(463, 438, 494, 496, 448, 603, 596, 616, 633,
                        608, 471, 481, 449, 443, 456)/100,
                 y2 = c(950, 890, 1010, 1230, 940, 1080, 1050, 1080,
                        1190, 1080, 960, 930, 870, 820, 910)/1000)
str(da)

# ATTENTION: This is a toy example, small dataset easy to handle.
xyplot(y1 ~ y2, groups = trt, data = da)

# Fitting the manova model.
m_full <- lm(cbind(y1, y2) ~ trt, data = da)

# Separated univariate models.
summary(m_full)

# Separated ANOVA tables.
summary.aov(m_full)

# MANOVA table.
anova(m_full)

# Estimated parameters.
coef(m_full)

# Estimated covariance matrix (least squares).
var(residuals(m_full))

# Null model.
m_null <- update(m_full, . ~ 1)

#-----------------------------------------------------------------------
# Doing the math.

# Error SSP of the full model.
E_full <- crossprod(residuals(m_full))
E_full

# Error SSP of the null model.
E_null <- crossprod(residuals(m_null))
E_null

# Extra error SSP = SSPnull - SSPfull.
E_extr <- E_null - E_full
E_extr

# Wilks' lambda.
det(E_full)/det(E_full + E_extr)
anova(m_full, test = "Wilks")

# Roy's largest root.
eigen(solve(E_full, E_extr))$values[1]
anova(m_full, test = "Roy")

# Hotelling-Lawley (Wald)
sum(diag(solve(E_full, E_extr)))
anova(m_full, test = "Hotelling-Lawley")

# Pillai.
sum(diag(E_extr %*% solve(E_extr + E_full)))
anova(m_full, test = "Pillai")

#-----------------------------------------------------------------------
# Calculating using matricial algebra.

# Matrix of responses.
Y <- as.matrix(da[, 2:3])
colSums(Y)

# Total SSP = Y'Y.
t(Y) %*% Y

# Design matrix (full model).
X <- model.matrix(~trt, data = da)
dim(X)

# Estimation of the model parameters.
XlX <- crossprod(X)
XlY <- crossprod(X, Y)
B <- solve(XlX, XlY)
B

# Just to confirm.
coef(m_full)

# To create projection matrices.
proj <- function(X) {
    tX <- t(X)
    X %*% solve(tX %*% X) %*% tX
}

# Projection (hat) matrix for the full model.
Hfull <- proj(X)

# Projection (hat) matrix for the null model.
Hnull <- proj(X[, 1])

# Identity.
I <- diag(nrow(X))

# Error SSP of the full model.
E_full <- t(Y) %*% (I - Hfull) %*% Y
E_full

# Error SSP of the null model.
E_null <- t(Y) %*% (I - Hnull) %*% Y
E_null

# Extra error SSP.
E_extr <- t(Y) %*% (Hfull - Hnull) %*% Y
E_extr

# Wilks' lambda.
det(E_full)/det(E_full + E_extr)
anova(m_full, test = "Wilks")

# Roy's largest root.
eigen(solve(E_full, E_extr))$values[1]
anova(m_full, test = "Roy")

# Hotelling-Lawley (Wald)
sum(diag(solve(E_full, E_extr)))
anova(m_full, test = "Hotelling-Lawley")

# Pillai.
sum(diag(E_extr %*% solve(E_extr + E_full)))
anova(m_full, test = "Pillai")

#-----------------------------------------------------------------------
