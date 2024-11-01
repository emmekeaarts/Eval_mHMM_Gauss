library(ggplot2)

a <- 2.5
b <- 2.25
sample <- 1/rgamma(10000, shape = a, rate = b)
# x11()
ggplot(data.frame(x = sample), aes(x=x)) + geom_density() +
  ggtitle(paste0("Shape ", a, ", rate ", b)) + xlim(c(0, 10)) + xlab("Variance")

median(sample)
sqrt(median(sample))

mean(sample)
sqrt(mean(sample))
