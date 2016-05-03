# CREATED  10 October 2014
# MODIFIED 10 October 2014


# Write the fixed parameters file for input into the delay difference
write.inputfile <- function(par){

    cat("# Brody growth coefficients", file = "FixedWeeklyParameters.txt", sep="\n")
    cat(par[1], file = "FixedWeeklyParameters.txt", append = TRUE, sep="\n")
    cat("# Estimated weight at recruitment at 22 weeks (in kg)", file = "FixedWeeklyParameters.txt", append = TRUE, sep="\n")
    cat(par[2], file = "FixedWeeklyParameters.txt", append = TRUE, sep="\n")    
    cat("# Parameter defining weight one timestep before recruitment (at week 21, in kg)", file = "FixedWeeklyParameters.txt", append = TRUE, sep="\n")
    cat(par[3], file = "FixedWeeklyParameters.txt", append = TRUE, sep="\n")    
    cat("# Catchability scaling factor", file = "FixedWeeklyParameters.txt", append = TRUE, sep="\n")
    cat("1e-4", file = "FixedWeeklyParameters.txt", append = TRUE, sep="\n")
    cat("# Biomass scaling factor", file = "FixedWeeklyParameters.txt", append = TRUE, sep="\n")
    cat("1e5", file = "FixedWeeklyParameters.txt", append = TRUE, sep="\n")
    cat("# Recruitment scaling factor", file = "FixedWeeklyParameters.txt", append = TRUE, sep="\n")
    cat("1e7", file = "FixedWeeklyParameters.txt", append = TRUE, sep="\n")
    cat("# Number of intervals in a year", file = "FixedWeeklyParameters.txt", append = TRUE, sep="\n")
    cat("52", file = "FixedWeeklyParameters.txt", append = TRUE, sep="\n")

}

#  This function create somatic growth parameters based on some input
GenerateSomaticGrowthParameters <- function(vb.Linf.mult.fac = 1, vb.k.mult.fac = 1, rec.age = 22){

# Von Bertalanffy growth parameters and length-weight relationship from Gribble and Dredge (1994)
# 1990 data
Linf.m <- 37.8 * vb.Linf.mult.fac; Linf.f <- 46.3 * vb.Linf.mult.fac# Carapace length, mm
K.m <- 0.0528 * vb.k.mult.fac; K.f <- 0.0436 * vb.k.mult.fac; # per week

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

# Let's fit Schnute fixing weight at recruitment
Eq114.bis <- function(par, x, k) {par[1] + (average.weight.at.age[rec.age] - par[1]) * (1-par[2]^(1+x-k))/(1-par[2])}
ssq.Eq114.TigerPrawn.2par <- function(par){ sum((average.weight.at.age[seq(rec.age, 200)] - Eq114.bis(par, seq(rec.age, 200), k = rec.age))^2)}
(Schnute.results.for.TigerPrawn.2par <- optim(c(1, 0.5), ssq.Eq114.TigerPrawn.2par))

### Parameters for the Delay difference model
print("Parameters for the Delay difference model")
print(paste("Rho",  round(Schnute.results.for.TigerPrawn.2par$par[2],5)))
print(paste("Weight at age", rec.age, "months is", round(average.weight.at.age[rec.age],3)))

print(paste("And parameter weight 1 timestep before recruitment", round(Schnute.results.for.TigerPrawn.2par$par[1],3)))

return(c(rho=round(Schnute.results.for.TigerPrawn.2par$par[2],5), wk=1e-3 * as.numeric(round(average.weight.at.age[rec.age],3)), wk.1=1e-3 * round(Schnute.results.for.TigerPrawn.2par$par[1],3)))
}
