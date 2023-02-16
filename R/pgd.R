n <- 100
set.seed(1234)
x1 <- rnorm(n, 3)
x2 <- rbinom(n, 1, 0.5)
b <- c(2, 0.5, 1.5)

y <- b[1] + b[2]*x1 + b[3]*x2 + rnorm(n, 0.75)
hist(y)
reg <- lm(y ~ x1 + x2)
reg

