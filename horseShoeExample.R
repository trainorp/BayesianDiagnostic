library(lars)
library(rstan)
setwd("~/gdrive/Dissertation/Aim3")

data('diabetes')

# Divide diabetes dataset into a training set and a disjoint test
# set of szie 100
set.seed(1234)
test.ixs <- sample(nrow(diabetes), 100)
test.data <- diabetes[test.ixs,]
train.data <- diabetes[-test.ixs,]

# OLS ###
ols.fit <- lm(y ~ x2, data = train.data)
coef(ols.fit)
mean((predict(ols.fit, newdata = test.data) - test.data$y)^2)

# LASSO ####
lasso.fit <- with(train.data, glmnet::cv.glmnet(x2, y))
coef(lasso.fit) # Same as s = "lambda.1se"
mean((predict(lasso.fit, test.data$x2) - test.data$y)^2)

# Horeshoe (Stan) ####
# This seems to produce poor results; for example, low n_eff (sometimes below 5)
# and (not surprisingly) significant variation in the mean beta from run to run.
# seed=123 => beta[3]: mean=0.96, n_eff=21
# seed=124 => beta[3]: mean=202.17, n_eff=2
stan.fit <- stan(file = 'horseshoe.stan', 
                 data = list(n = length(train.data$y), 
                             y = train.data$y,
                             p = ncol(train.data$x2), 
                             X = train.data$x2),
                 seed = 123)
