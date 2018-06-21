library(tidyverse)
library(rstan)
library(brms)

rstan_options(auto_write=TRUE)
options(mc.cores=parallel::detectCores()) # Run on multiple cores

set.seed(3875)

ir<-data.frame(scale(iris[,-5]),Species=iris[,5])

ptm<-proc.time()
b1<-brm(Species~Petal.Length+Petal.Width+Sepal.Length+Sepal.Width,data=ir,
        family="categorical",chains=4,iter=2000)
proc.time()-ptm

launch_shinystan(b1)
