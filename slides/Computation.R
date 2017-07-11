data(iris)
head(iris)

# Response variable matrix
Y <- as.matrix(iris[,1:4],ncol = 4, nrow = 150)

# Design matrix
X = model.matrix(~Species, data = iris)

# Regression coefficients
B <- solve(tcrossprod(t(X)))%*%t(X)%*%Y
B

# Covariance matrix
Sigma <- t(Y - X%*%B)%*%(Y - X%*%B)/dim(iris)[1]
Sigma

# Correlation matrix
cov2cor(Sigma)

# Stacked vector of response variable
YY <- c(Y)

# Design matrix
XX <- bdiag(X,X,X,X)

# MLE for beta
MLE_beta <- solve(tcrossprod(t(XX)))%*%t(XX)%*%YY
MLE_beta

# MLE for Omega
Omega <- (YY - XX%*%MLE_beta)%*%t(YY - XX%*%MLE_beta)
t(XX)%*%solve(Omega)%*%XX

II <- Diagonal(150,1)
R <- (YY - XX%*%MLE_beta)
Omega <- kronecker(Sigma, II)
solve(Omega)%*%solve(Omega)

F_Sigma <- solve(Sigma)%*%solve(Sigma)/2
solve(F_Sigma)/150

require(mcglm)
Z0 <- mc_id(iris)
tt = mcglm(c(Sepal.Length ~ Species, Sepal.Width ~ Species, 
        Petal.Length ~ Species, Petal.Width ~ Species),
        matrix_pred = list(Z0,Z0,Z0,Z0), data = iris)
summary(tt)
vcov(tt)[19:22,19:22]
