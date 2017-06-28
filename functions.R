# Auxiliary functions ------------------------------------------------
# Author: Wagner Hugo Bonat LEG/UFPR ---------------------------------
# Date: 14/06/2017 ---------------------------------------------------

# logLik function for object of mlm class ----------------------------

logLik.mlm <- function(object,...)
{
  resids <- residuals(object)
  n <- nrow(resids)
  Sigma_ML <- crossprod(resids) /n
  ans <- sum(dmvnorm(resids, sigma=Sigma_ML, log=T))
  
  df <- length(coef(object)) + nrow(Sigma_ML) * (nrow(Sigma_ML) + 1) / 2
  attr(ans, "nobs") <- n
  attr(ans, "df") <- df
  class(ans) <- "logLik"
  ans
}

# Build F matrix for Wald multivariate test --------------------------
build_F <- function(vector) {
  FF <- diag(length(vector))
  FF_list <- by(FF, vector, as.matrix)
  return(FF_list)
}

# Auxiliary Kronecker product ----------------------------------------
aux_kronecker <- function(FF, G) {
  kronecker(G,FF)
}


# manova for objects of mcglm class ----------------------------------
manova.mcglm <- function(obj, ...) {
  beta <- coef(obj, type = "beta")[,1]
  n_beta <- length(beta)
  VCOV <- vcov(obj)[1:n_beta, 1:n_beta]
  ## Conditional variance-covariance model
  FF <- build_F(vector = attr(obj$list_X[[1]], "assign"))
  G <- Diagonal(length(obj$mu_list), 1)
  CC <- lapply(FF, function(x, G){kronecker(G,x)}, G = G)
  N <- obj$n_obs
  test_W <- c()
  df <- c()
  p_value <- c()
  for(i in 1:length(CC)) {
    test_W[i] <- as.numeric(t(CC[[i]]%*%beta)%*%solve(CC[[i]]%*%VCOV%*%t(CC[[i]]))%*%(CC[[i]]%*%beta))
    df[i] <- dim(CC[[i]])[1]
    p_value[i] <- pchisq(test_W[i], df = df[i], lower.tail = FALSE)
  }
  names <- c("Intercept", attr(terms(obj$linear_pred[[1]]), "term.labels"))
  out <- data.frame("Effects" = names, "Df" = df, 
                    "Hotelling-Lawley" = round(test_W/N,3),
                    "Qui-square" = round(test_W,3),
                    "p_value" = p_value)
  return(out)
}


