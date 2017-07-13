#-----------------------------------------------------------------------
# Packages.

library(lattice)

#-----------------------------------------------------------------------
# Pimentel Gomes, page 254.

y1 <- c(463, 438, 494, 496, 448, 603, 596, 616, 633, 608, 471, 481, 449,
        443, 456)/100
y2 <- c(950, 890, 1010, 1230, 940, 1080, 1050, 1080, 1190, 1080, 960,
        930, 870, 820, 910)/1000
trt <- gl(n = 3, k = 5, labels = c("Test", "TurFer", "TurNat"))
da <- data.frame(trt, y1, y2)
str(da)

# ATTENTION: This is a toy example, small dataset easy to handle.
xyplot(y1 ~ y2, groups = trt, data = da)

# Fitting the manova model.
m0 <- lm(cbind(y1, y2) ~ trt, data = da)

# Separated univariate models.
summary(m0)

# Separated ANOVA tables.
summary.aov(m0)

# MANOVA table.
anova(m0)

# Estimated parameters.
coef(m0)

# Estimated covariânce matrix (least squares).
var(residuals(m0))

#-----------------------------------------------------------------------
# Doing the math.

# Matrix of responses.
Y <- as.matrix(da[, grep("^y\\d+", names(da))])
colSums(Y)

# Total SSP = Y'Y.
t(Y) %*% Y

# Design matrix.
X <- model.matrix(~trt, data = da)
dim(X)

# Estimation of the model parameters.
XlX <- crossprod(X)
XlY <- crossprod(X, Y)
B <- solve(XlX, XlY)
B

proj <- function(X) {
    tX <- t(X)
    X %*% solve(tX %*% X) %*% tX
}

# Projection matrices.
H <- proj(X)
Hmu <- proj(X[, 1])
I <- diag(nrow(H))

# SSP of the full model.
SSP_full <- t(Y) %*% (I - H) %*% Y

# SSP of the restricted model.
SSP_rest <- t(Y) %*% (I - Hmu) %*% Y

# Extra SSP = SSPrest - SSPfull.
SSP_extr <- t(Y) %*% (H - Hmu) %*% Y

# Just to confirm.
SSP_extr - (SSP_rest - SSP_full)

# Wilks' lambda.
det(SSP_full)/det(SSP_full + SSP_extr)
anova(m0, test = "Wilks")

# Roy's largest root.
eigen(solve(SSP_full, SSP_extr))$values[1]
anova(m0, test = "Roy")

# Hotelling-Lawley (Wald)
sum(diag(solve(SSP_full, SSP_extr)))
anova(m0, test = "Hotelling-Lawley")

# Pillai.
sum(diag(SSP_extr %*% solve(SSP_extr + SSP_full)))
anova(m0, test = "Pillai")

#-----------------------------------------------------------------------
