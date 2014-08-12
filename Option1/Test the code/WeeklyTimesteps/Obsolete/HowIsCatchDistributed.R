Obs <- log(rowSums(noisy.Catch))
Mod <- log(rowSums(Catch))

x <- Obs - Mod

lik.fct <- function(par) -sum(dnorm(x, mean = par[1], sd = par[2], log = TRUE))

solution <- optim(c(0,100), lik.fct)
print(solution)

hist(x, freq = FALSE)
curve(dnorm(x, mean = solution$par[1], sd = solution$par[2]), col = "blue", add = TRUE)
