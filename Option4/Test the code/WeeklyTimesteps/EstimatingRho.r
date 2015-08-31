# CREATED  11 February 2013
# MODIFIED 20 August 2013

# AUTHOR marco.kienzle@gmail.com

# PURPOSE estimate growth parameter rho for EKP to use in a weekly delay difference model

# REFERENCE Deriso, 1980 Harvesting strategies and parameter estimation for an age-structured model
#           Schnute, 1985 A general theory for analysis of catch and effort data
#           Hilborn and Walters, 1992 Quantitative fisheries stock assessment
#           Quinn and Deriso, 1999 Quantitative fish dynamics

## # Schnute's growth model (Eq. 1.14)
Eq114 <- function(par, x, k) {par[1] + (par[2] - par[1]) * (1-par[3]^(1+x-k))/(1-par[3])}

# Von Bertalanffy growth parameters and length-weight relationship from Gribble and Dredge (1994)
# 1990 data
Linf.m <- 37.8; Linf.f <- 46.3 # Carapace length, mm
K.m <- 0.0528; K.f <- 0.0436; # per week

# Length-weight relationship (from carapace length to grams)
a.m <-  0.0022; a.f <-  0.0023
b.m <-  2.7557; b.f <- 2.7292

# Von Bertalanffy growth function for length
vbgf <- function(t, Linf, k, t0=0) Linf * (1-exp(-k * (t -t0)))

# Create a vector of length at age
length.at.age <- data.frame(Female = rep(NA,200), Male = rep(NA,200)); dimnames(length.at.age)[[1]] <- seq(1,200);
length.at.age$Female <- vbgf(seq(1,200), Linf.f, K.f, t0 = 0)
length.at.age$Male <- vbgf(seq(1,200), Linf.m, K.m, t0 = 0)

# Weight at age
weight.at.age <- data.frame(Female = rep(NA,200), Male = rep(NA,200)); dimnames(weight.at.age)[[1]] <- seq(1,200);
weight.at.age$Female <- a.f * length.at.age$Female ^ b.f
weight.at.age$Male <- a.m * length.at.age$Male ^ b.m

average.weight.at.age <- rowMeans(weight.at.age)
#print(average.weight.at.age)

# Plot weight at age
plot(seq(1,200), weight.at.age$Female, pch = "F", xlab = "Time (months)", ylab = "Weight (g)", las = 1)
points( seq(1,200), weight.at.age$Male, pch = "M")
points(seq(1,200), average.weight.at.age, pch = 19)

# Estimate Ford-Walford parameters using the whole range of ages ?
# NOTE that Ford-Walford graphs are used for length data ( see Ricker, 1958, Eq. 9.8)
(lm1 <- lm(y~x, data = data.frame(x = average.weight.at.age[seq(1, 29)], y = average.weight.at.age[seq(2, 30)])));

plot(average.weight.at.age[seq(1, 199)], y = average.weight.at.age[seq(2, 200)],
xlab = "Average weight at month t", ylab = "Average weight at month t+1", xlim = c(0,80), ylim = c(0,80), las = 1)
abline(0,1)

abline(lm1, col = "blue");

# What are the parameter of this linear relationship applied for age >= recruitment?
rec.age <- 22 # here the age at recruitment is chosen so that size at recruitment ~= 7.7 grams

(lm2 <- lm(y~x, data = data.frame(x = average.weight.at.age[seq(rec.age, 199)], y = average.weight.at.age[seq(rec.age+1, 200)])));
abline(lm2, col = "green");

# Now, let's apply Schnute method to the TigerPrawn dataset

# Schnute's growth model (Eq. 1.14)
ssq.Eq114.TigerPrawn <- function(par){ sum((average.weight.at.age[seq(rec.age, 200)] - Eq114(par, seq(rec.age, 200), k = rec.age))^2)}
(Schnute.results.for.TigerPrawn <- optim(c(1, 1, 0.5), ssq.Eq114.TigerPrawn))
# NOTE: you end-up with an estimate of the weight at recruitment which is quite differenct that average.weight.at.age[rec.age] NOT GOOD !

# Let's fit Schnute fixing weight at recruitment
Eq114.bis <- function(par, x, k) {par[1] + (average.weight.at.age[rec.age] - par[1]) * (1-par[2]^(1+x-k))/(1-par[2])}
ssq.Eq114.TigerPrawn.2par <- function(par){ sum((average.weight.at.age[seq(rec.age, 200)] - Eq114.bis(par, seq(rec.age, 200), k = rec.age))^2)}
(Schnute.results.for.TigerPrawn.2par <- optim(c(1, 0.5), ssq.Eq114.TigerPrawn.2par))

# Let's fit Schnute fixing weight at recruitment and fixing rho
Eq114.tris <- function(par, x, k, rho) {par[1] + (average.weight.at.age[rec.age] - par[1]) * (1-rho^(1+x-k))/(1-rho)}
ssq.Eq114.TigerPrawn.1par <- function(par){ sum((average.weight.at.age[seq(rec.age, 200)] - Eq114.tris(par, seq(rec.age, 200), k = rec.age, rho = exp(-0.1794067)))^2)}
(Schnute.results.for.TigerPrawn.1par <- optimize(ssq.Eq114.TigerPrawn.1par,lower = -10, upper = 1e3))

## How does it compares to Von Bertalanffy estimates
# The curve are sigmoidal so we have to consider only age above which Von Bertalanffy growth function for length might apply (age >= 3)

# Von Bertalanffy growth function for weight
vbgf <- function(t, Linf, k, t0=0) 0.5*(a.m + a.f) * (Linf * (1-exp(-k * (t -t0))))^(0.5 * (b.m+b.f))

# Estimate parameters
(nls.model <- nls( y ~ vbgf( x, a, b,c), data = data.frame(x = seq(rec.age,200), y = average.weight.at.age[seq(rec.age, 200)]),
start = list(a = 50, b = 0.03, c= 1)))

# Plot weight at age
plot(seq(1,200), weight.at.age$Female, pch = 1, xlab = "Time (weeks)", ylab = "Weight (g)", las = 1)
points( seq(1,200), weight.at.age$Male, pch = 2)
points(seq(1,200), average.weight.at.age, pch = 19)

curve( vbgf(x, coef(nls.model)[1], coef(nls.model)[2], coef(nls.model)[3]), from = rec.age, to = 200, lwd = 2, col = "blue", add = TRUE)
#curve( Eq114(Schnute.results.for.TigerPrawn$par, x, k=rec.age), from = rec.age, to = 200, lwd = 2, col = "grey", add = TRUE)
curve( Eq114.bis(Schnute.results.for.TigerPrawn.2par$par, x, k=rec.age), from = rec.age, to = 200, lwd = 2, col = "red", add = TRUE)
#curve( Eq114.tris(Schnute.results.for.TigerPrawn.1par$minimum, x, k=rec.age, rho = exp(-0.1794067)), from = rec.age, to = 200, lwd = 2, col = "green", add = TRUE)

legend(55, 25, pch = c(2,1, 19, NA, NA, NA), lty = c(NA, NA, NA, 1, 1,1), col = c("black", "black", "black", "blue", "red","yellow"), lwd = c(NA,NA,NA,2,2,2), 
legend = c("VBGF for male", "VBGF for female", "VBGF both sex combined", "VBGF  model for combined sex", "Schnute", "Deriso"))

# Use Deriso formula (Equation 1)

Deriso.W <- function(rho, t, weight.at.recruitment, age.at.recruitment){

weight.at.recruitment * ( 1 - rho ^ (1 + t - age.at.recruitment)) / (1 - rho)
}

ssq.Deriso.Eq <- function(par){
 sum((average.weight.at.age[seq(rec.age+1, 200)] - Deriso.W(par, seq(rec.age+1, 200), weight.at.recruitment = average.weight.at.age[rec.age],
 age.at.recruitment = rec.age))^2)
}
Deriso.opt <- optimize(ssq.Deriso.Eq, lower = 0.1, upper = 0.99)

curve(Deriso.W(rho =  Deriso.opt$minimum, x,  average.weight.at.age[rec.age], rec.age), from = rec.age, to =200, col = "yellow", lwd = 2, add = T)


### Parameters for the Delay difference model
print("Parameters for the Delay difference model")
print(paste("Rho",  round(Schnute.results.for.TigerPrawn.2par$par[2],5)))
print(paste("Weight at age", rec.age, "months is", round(average.weight.at.age[rec.age],3)))

print(paste("And parameter weight 1 month before recruitment", round(Schnute.results.for.TigerPrawn.2par$par[1],3)))

