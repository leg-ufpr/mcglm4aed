# Example 3: Comparing aesthetic eyelid surgery techniques -----------
data <- read.table("plastica.csv", header = TRUE, sep = ",")
names(data) <- c("ID", "Technique", "Time", "CEN", "LAT", "MED")
data$time <- as.factor(rep(c(0,3,12,120), 60))
data$Time <- as.factor(data$Time)
data$ID <- as.factor(data$ID)
levels(data$Technique) <- c("Blepharoplasty", "Endoscopic", 
                            "Endoscopic + Blepharoplasty")
levels(data$Time) <- c("Preoperative", "30 days",
                       "90 days", "10 Years")

# Dropping missing data 
cod_na <- as.numeric(rownames(data[is.na(data$MED),]))
data_na <- data[-cod_na,]
summary(data_na)

# Matrix linear predictor
Z0 <- mc_id(data_na)

# Technique
ex3_lm1 <- manova(cbind(LAT, MED, CEN) ~ Technique, 
                  data = data_na)
ex3_mc1 <- mcglm(c(LAT ~ Technique, MED ~ Technique, CEN ~ Technique), 
                 matrix_pred = list(Z0,Z0,Z0), data = data_na, 
                 control_algorithm = list(correct = FALSE))

summary(ex3_lm1, test = "Hotelling-Lawley")
manova.mcglm(ex3_mc1)

# Technique + Time
ex3_lm2 <- manova(cbind(LAT, MED, CEN) ~ Technique + Time, 
                  data = data_na)
ex3_mc2 <- mcglm(c(LAT ~ Technique + Time, MED ~ Technique + Time, 
                   CEN ~ Technique + Time), 
                 matrix_pred = list(Z0,Z0,Z0), data = data_na, 
                 control_algorithm = list(correct = FALSE))

summary(ex3_lm2, test = "Hotelling-Lawley")
manova.mcglm(ex3_mc2)

## Technique + Time + Technique*Time
ex3_lm3 <- manova(cbind(LAT, MED, CEN) ~ Technique*Time, data = data_na)
ex3_mc3 <- mcglm(c(LAT ~ Technique*Time, MED ~ Technique*Time, 
                   CEN ~ Technique*Time), matrix_pred = list(Z0,Z0,Z0), 
                 data = data_na, control_algorithm = list(correct = FALSE))

summary(ex3_lm3, test = "Hotelling-Lawley")
manova.mcglm(ex3_mc3)

c(coef(ex3_lm1),rep(NA,27))/c(coef(ex3_mc1, type = "beta")[,1],rep(NA,27))
c(coef(ex3_lm2),rep(NA,18))/c(coef(ex3_mc2, type = "beta")[,1],rep(NA,18))
c(coef(ex3_lm3))/coef(ex3_mc3, type = "beta")[,1]

c(coef(ex3_lm1),rep(NA,27))/c(coef(ex3_lm3))
c(coef(ex3_lm2),rep(NA,18))/c(coef(ex3_lm3))

c(coef(ex3_mc1,type = "beta")[,1],rep(NA,27))/c(coef(ex3_mc3,type = "beta")[,1])
c(coef(ex3_mc2,type = "beta")[,1],rep(NA,18))/c(coef(ex3_mc3,type = "beta")[,1])

c(coef(ex3_mc1,type = "beta")[,1],rep(NA,27))/c(coef(ex3_lm3))
c(coef(ex3_mc2,type = "beta")[,1],rep(NA,18))/c(coef(ex3_lm3))
c(coef(ex3_mc3,type = "beta")[,1])/c(coef(ex3_lm3))

# By hand
X <- model.matrix(~ Technique*Time, data = data_na)
S <- t(cbind(data_na$LAT, data_na$MED, data_na$CEN) - fitted(ex3_mc3))%*%(cbind(data_na$LAT, data_na$MED, data_na$CEN) - fitted(ex3_mc3))/213
VVV <- kronecker(S, solve(t(X)%*%X))
round(VVV - vcov(ex3_mc3)[1:36,1:36], 2)
round(VVV - vcov(ex3_lm3), 2)

# Technique effect
FF <- build_F(vector = c(0,1,1,2,2,2,3,3,3,3,3,3))
CC <- kronecker(diag(3), FF[[2]])
beta <- coef(ex3_mc3, type = "beta")[,1]
t(CC%*%beta)%*%solve(CC%*%(kronecker(S,solve(t(X)%*%X)))%*%t(CC))%*%(CC%*%beta)/dim(data_na)[1]

# Using another approach
B <- coef(ex3_lm3)
FF <- FF[[2]]
G <- diag(3)
A <- t(FF%*%B%*%G)%*%solve(FF%*%solve(t(X)%*%X)%*%t(FF))%*%(FF%*%B%*%G)
B <- t(G)%*%S%*%G
sum(diag(A%*%solve(B)))/213 # Exactly what I got !!

summary(ex3_lm3, test = "Hotelling-Lawley")
manova.mcglm(ex3_mc3)

# Conditional test
Sigma <- vcov(ex3_mc3)[1:36,1:36]
fixed <- c(1,13,25) # Intercepts
test <- c(2:3,14:15,26:27)
beta_test <- coef(ex3_mc3)[test,1]
Sigma11 <- Sigma[test, test]
Sigma12 <- Sigma[test, fixed]
Sigma21 <- Sigma[fixed, test]
Sigma22 <- Sigma[fixed, fixed]
Sigma.cond <- Sigma11 - Sigma12 %*% solve(Sigma22) %*% Sigma21

# Time effect
FF <- build_F(vector = c(0,1,1,2,2,2))[[3]]
CC <- kronecker(diag(3), FF)
beta <- coef(ex3_mc2, type = "beta")[,1]
t(CC%*%beta)%*%solve(CC%*%(kronecker(S,solve(t(X)%*%X)))%*%t(CC))%*%(CC%*%beta)/dim(data_na)[1]

summary(ex3_lm2, test = "Hotelling-Lawley")
manova.mcglm(ex3_mc2)