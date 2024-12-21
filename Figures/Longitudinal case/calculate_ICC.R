library(data.table)
set.seed(1)
n = 100
data.dt = data.table(cluster = sample(c(1,2,3), n, replace = T))
data.dt[,outcome := rnorm(n, mean = cluster, sd = 1)]
data.dt[,cluster := factor(cluster)]

summary_aov = summary(aov(outcome ~ cluster,data=data.dt))

summary_aov[[1]][1,2]/sum(summary_aov[[1]][,2])

#------------------------------------------------

#install.packages("irr")
library(irr)

data("anxiety", package = "irr")
head(anxiety, 4)

icc(
  anxiety, model = "twoway", 
  type = "agreement", unit = "single"
)