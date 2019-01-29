library(ggplot2)
library(cowplot)

# sample a discrete distribution from start to finish
sampleDist <- function(n, start, finish){
  leng <- finish - start + 1
  sample(x = seq(start, finish, by = 1), n, replace = T, prob = rep(1/leng, times = leng))
}
