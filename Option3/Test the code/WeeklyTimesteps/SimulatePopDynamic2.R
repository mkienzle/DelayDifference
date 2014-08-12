# CREATED  11 May 2012
# MODIFIED 29 May 2014

# PURPOSE simulate fishery data to test the weekly delay difference model

# Useful libraries
library(chron)
library("RcppArmadillo")
library ("inline")
source("RcppArmadilloFastLoop.R")

##### Parameters
years <- seq(1956,2010)
intra.year.timesteps <- paste("Week", 1:52, sep = "")
nb.age.groups <- seq(22,52 * 3) # Age-groups from 22 weeks old to 3 years old

random.catchability <- runif(1, min = 1, max = 10)
catchability.q <- random.catchability * 1e-4
#catchability.q <- 5 * 1e-5
write.table(paste(" ", round(catchability.q * 1e4,2), " "), file = "Results/SimulatedParameters.txt", append = T, eol="\n", col.names = FALSE, row.names = FALSE, quote = FALSE)

write.table(round(catchability.q * 1e4,2), file = "Results/SimPar.txt", eol="\n", col.names = FALSE, row.names = FALSE, quote = FALSE)
#nat.mortality <- 0
#nat.mortality <- 0.195 * 12 / 52 # per week
nat.mortality <- runif(1, min = 0.025, max = 0.065)
write.table(round(nat.mortality,4), file = "Results/SimPar.txt", append = TRUE, eol="\n", col.names = FALSE, row.names = FALSE, quote = FALSE)

# Proportion of spawners in each month (fixed to 60%)
prop.spawners <- 0.6

recruitment.by.week <- rep(NA, length(intra.year.timesteps));
# Recruitment pattern generated using the von Mises dist with rbar = 13 and k = 5 (using the WeeklyDelayDifference.m script)
#rec.pat <- c(0.0013, 0.0023, 0.0042, 0.0072, 0.0121, 0.0194, 0.0298, 0.0432, 0.0591, 0.0757, 0.0906, 0.1010, 0.1048, 0.1010, 0.0906, 0.0757, 0.0591, 0.0432, 0.0298, 0.0194, 0.0121, 0.0072, 0.0042, 0.0023, 0.0013, 0.0007, 0.0004, 0.0002, 0.0001, 0.0001, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0001, 0.0001, 0.0002, 0.0004, 0.0007)
#rec.pat <- c(1,rep(0,51))
#rec.pat <- c(rep(0,12), 1, rep(0,39))
# Recruitment pattern generated using the von Mises dist with rbar = 50 and k = 3 (using the WeeklyDelayDifference.m script)
#rec.pat <- c(0.0651, 0.0561, 0.0465, 0.0372, 0.0288, 0.0217, 0.0159, 0.0114, 0.0081, 0.0057, 0.0039, 0.0027, 0.0019, 0.0014, 0.0010, 0.0007, 0.0005, 0.0004, 0.0003, 0.0003, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0002, 0.0003, 0.0003, 0.0004, 0.0005, 0.0007, 0.0010, 0.0014, 0.0019, 0.0027, 0.0039, 0.0057, 0.0081, 0.0114, 0.0159, 0.0217, 0.0288, 0.0372, 0.0465, 0.0561, 0.0651, 0.0725, 0.0774, 0.0791, 0.0774, 0.0725)
#recruitment.by.week[1:52] <- rec.pat / sum(rec.pat) # Must sum to 1

# Parameters
#A <- 0
#B <- 10
A <- runif(1, min = -pi, max = pi)
B <- runif(1, min = 1, max = 60)

write.table(round(A,2), append = TRUE, file = "Results/SimPar.txt", eol="\n", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(round(B,2), append = TRUE, file = "Results/SimPar.txt", eol="\n", col.names = FALSE, row.names = FALSE, quote = FALSE)

library(circular)
# Interval of time you are interested in
boundaries <- seq(-pi, pi, length =53)

precision <- 1e3
res <- rep(NA, 52); 
for(i in 1:52) res[i] <- sum(dvonmises(seq(boundaries[i], boundaries[i+1] - 1 / precision, 1 / precision), A, B))/precision
recruitment.by.week <- res # / sum(res)

#total.recruitment.in.year <- 1e1 * rep(1, length(years))
#total.recruitment.in.year <- 1e9 * c(rep(1, length(years)-10), seq(3,0.2, length = 5), seq(0.2, 1, length = 5))
rand.rec <- runif(10, min = 0.1, max = 10)
total.recruitment.in.year <- 1e7 * c(rep(1, length(years)-10), rand.rec)
write.table(round(rand.rec,2), append = TRUE, file = "Results/SimPar.txt", eol="\n", col.names = FALSE, row.names = FALSE, quote = FALSE)
#total.recruitment.in.year <- 1e9 * c(rep(1, length(years)-10), seq(3,0.2, length = 5), seq(0.2, 1, length = 5))
#total.recruitment.in.year <- 1e9 * c(rep(1, length(years)-10), c(3,0.25,2,0.2,2,0.1,4,3,2,1))
#total.recruitment.in.year <- 1e9 * rep(1, length(years))

# Weight at age according to Schnute's model
# round(Eq114.bis(Schnute.results.for.TigerPrawn.2par$par, seq(13,3*52), k=rec.age),2)
source("EstimatingRho.r")
weight.at.age <- round(Eq114.bis(Schnute.results.for.TigerPrawn.2par$par, seq(22,3*52), k=rec.age),2)
#weight.at.age <- c(7.66,  9.54, 11.36, 13.12, 14.82, 16.47, 18.07, 19.61, 21.10, 22.55, 23.95, 25.30,
#26.61, 27.88, 29.10, 30.29, 31.44, 32.55, 33.62, 34.66, 35.67, 36.64, 37.58, 38.49,
#39.38, 40.23, 41.05, 41.85, 42.63, 43.37, 44.10, 44.80, 45.48, 46.13, 46.77, 47.38,
#47.98, 48.55, 49.11, 49.65, 50.17, 50.67, 51.16, 51.63, 52.09, 52.53, 52.96, 53.37,
#53.77, 54.16, 54.53, 54.89, 55.25, 55.58, 55.91, 56.23, 56.54, 56.84, 57.12, 57.40,
#57.67, 57.93, 58.19, 58.43, 58.67, 58.90, 59.12, 59.33, 59.54, 59.74, 59.93, 60.12,
#60.30, 60.48, 60.65, 60.81, 60.97, 61.13, 61.28, 61.42, 61.56, 61.69, 61.83, 61.95,
#62.07, 62.19, 62.31, 62.42, 62.53, 62.63, 62.73, 62.83, 62.92, 63.01, 63.10, 63.19,
#63.27, 63.35, 63.42, 63.50, 63.57, 63.64, 63.71, 63.77, 63.84, 63.90, 63.96, 64.02,
#64.07, 64.13, 64.18, 64.23, 64.28, 64.32, 64.37, 64.41, 64.46, 64.50, 64.54, 64.58,
#64.61, 64.65, 64.68, 64.72, 64.75, 64.78, 64.81, 64.84, 64.87, 64.90, 64.93, 64.95,
#64.98, 65.00, 65.03, 65.05, 65.07, 65.09, 65.11, 65.13, 65.15, 65.17, 65.19, 65.21)

##### Data: use Moreton Tiger Prawn data as the basis of simulated datasets

# Moreton Bay brown tiger prawn data
dataset <- read.csv("Data/MBTigerCatchAndEffortByBiologicalYearAndWeek.csv") 

## Using constant effort from year to year
effort <- as.matrix(rep(subset(dataset, Biological.Year == 2004)$Total.Tiger.Effort, length(years)), ncol = 1)
#effort <- matrix(rep(0, 2860), ncol = 1)

##### Define the matrices for the population dynamic

Nb.at.age <- data.frame(matrix(1e-12, nrow = length(years) * length(intra.year.timesteps), ncol = length(nb.age.groups)))
dimnames(Nb.at.age) <- list(apply(expand.grid(intra.year.timesteps, years),1, paste, collapse = "."), nb.age.groups)

Weight.at.age <- outer(rep(1, length(years) * length(intra.year.timesteps)), weight.at.age)
dimnames(Weight.at.age) <- list(apply(expand.grid(intra.year.timesteps, years),1, paste, collapse = "."), nb.age.groups)

Nat.mortality <- matrix(nat.mortality, nrow = length(years) * length(intra.year.timesteps), ncol = length(nb.age.groups))
dimnames(Nat.mortality) <- list(apply(expand.grid(intra.year.timesteps, years),1, paste, collapse = "."), nb.age.groups)

Effort <- outer(effort[,1], rep(1, length(nb.age.groups)))
dimnames(Effort) <- list(apply(expand.grid(intra.year.timesteps, years),1, paste, collapse = "."), nb.age.groups)

Fish.mortality <- catchability.q * Effort

Recruitment <- matrix(0, nrow = length(years) * length(intra.year.timesteps), ncol = length(nb.age.groups))
dimnames(Recruitment) <- list(apply(expand.grid(intra.year.timesteps, years),1, paste, collapse = "."), nb.age.groups)
#Recruitment[,1] <- rep(total.recruitment.in.year, each = length(intra.year.timesteps)) * rep(recruitment.by.week, length(years))
#Recruitment[,1] <- rep(1e9 * recruitment.by.week, length(years))

# Varying recruitment
Recruitment[,1] <- c(outer(recruitment.by.week, total.recruitment.in.year))

## ##### Population dynamic
## for(i in seq(1,length(years) * length(intra.year.timesteps) - 1))
## {

##   # Apply mortality
##   Nb.at.age[i+1, seq(2, length(nb.age.groups))] <- Nb.at.age[i, seq(1, length(nb.age.groups) - 1)] * exp( - (Nat.mortality + Fish.mortality)[i, seq(1, length(nb.age.groups) - 1)])
##   # Add recruitment
##   Nb.at.age[i+1,] <- Nb.at.age[i+1,] + Recruitment[i+1,]
## #print(paste("##### ",i))
## }

##### Population dynamic
tmp <- FastLoop(as.matrix(Nb.at.age), as.matrix(Nat.mortality), as.matrix(Fish.mortality), as.matrix(Recruitment))
tmp <- as.data.frame(tmp)
dimnames(tmp) <- dimnames(Nb.at.age)
Nb.at.age <- tmp

# Calculate biomass in kg
Biomass <- Nb.at.age * 1e-3 * Weight.at.age

# Spawning stock biomass in kg
SSB <- tapply(rowSums(Nb.at.age * 1e-3 * Weight.at.age * prop.spawners * 0.5), rep(years, each = 52), sum)

# Calculate catch in kg
Catch <- Fish.mortality / (Fish.mortality + Nat.mortality) * Nb.at.age * 1e-3 * Weight.at.age * (1 - exp(-(Fish.mortality + Nat.mortality)))

# Monthly prawn number
par(mfrow = c(3,1))

#plot(seq(1,dim(Nb.at.age)[1]), rowSums(Nb.at.age), type = "b", xlab = "Time (months)", ylab = "Nb of prawn", axes = F)
#axis(1, at = seq(1,dim(Nb.at.age)[1])[seq(1, dim(Nb.at.age)[1],6)], label = dimnames(Nb.at.age)[[1]][seq(1, dim(Nb.at.age)[1],6)])
#axis(2, las = 1)
#box()
#abline(v = seq(1,200,12))

# Plot biomass
plot(rowSums(1e-3 * Biomass), type = "b", xlab = "Time (weeks)",
  ylab = "Biomass (tons)", axes = F, pch = 19)
axis(1, at = seq(1,dim(Nb.at.age)[1])[seq(1, dim(Nb.at.age)[1],6)], label = dimnames(Nb.at.age)[[1]][seq(1, dim(Nb.at.age)[1],6)])
axis(2, las = 1)
box()
abline(v = seq(1,3e3,length(intra.year.timesteps)))
#lines(with(subset(dataset, year %in% 1988:2000), CATCH/effort), col = "red", pch = "b")


# Plot catches
plot(rowSums(1e-3 * Catch), axes = F, xlab = "Time (weeks)", ylab = "Catch (tons)", type = "b", pch = 19)
axis(1, at = seq(1,dim(Nb.at.age)[1])[seq(1, dim(Nb.at.age)[1],6)], label = dimnames(Nb.at.age)[[1]][seq(1, dim(Nb.at.age)[1],6)])
axis(2, las = 1)
box()

# Plot CPUE
plot(rowSums(catchability.q * Biomass), type = "b", xlab = "Time (weeks)",
  ylab = "Catch per unit effort (kg per boat-days)", axes = F, pch = 19)
axis(1, at = seq(1,dim(Nb.at.age)[1])[seq(1, dim(Nb.at.age)[1],6)], label = dimnames(Nb.at.age)[[1]][seq(1, dim(Nb.at.age)[1],6)])
axis(2, las = 1)
box()
abline(v = seq(1,3e3,length(intra.year.timesteps)))
#lines(with(subset(dataset, year %in% 1988:2000), CATCH/effort), col = "red", pch = "b")

# Other graphical output
#library(lattice)
#xyplot(rowSums(Catch) ~ rep(1:12, 55) | substr(dimnames(Catch)[[1]], 5, 8))


# Output the data from 1958 to 2010

# Add a bit of random variability
noisy.Catch <- Catch * matrix( runif(length(years) * length(intra.year.timesteps) * length(nb.age.groups), min = 0.9, max = 1.1), nrow = length(years) * length(intra.year.timesteps), ncol = length(nb.age.groups))

noisy.Effort <- Effort * matrix( runif(length(years) * length(intra.year.timesteps) * length(nb.age.groups), min = 0.9, max = 1.1), nrow = length(years) * length(intra.year.timesteps), ncol = length(nb.age.groups))

#write.csv( file = "Data/SimulatedCPUE.csv", rowSums(catchability.q * Biomass)[-seq(1,2 * length(intra.year.timesteps))])
write.csv( file = "Data/SimulatedCPUE.csv", rowSums(noisy.Catch)[-seq(1, 2 * length(intra.year.timesteps))] / noisy.Effort[-seq(1,2 * length(intra.year.timesteps)), 1])
#write.csv( file = "Data/SimulatedEffort.csv", noisy.Effort[-seq(1, 2 * length(intra.year.timesteps)), 1])
write.csv( file = "Data/SimulatedEffort.csv", Effort[-seq(1, 2 * length(intra.year.timesteps)), 1])
#write.csv( file = "Data/SimulatedCatch.csv", rowSums(noisy.Catch)[-seq(1, 2 * length(intra.year.timesteps))])
write.csv( file = "Data/SimulatedCatch.csv", rowSums(Catch)[-seq(1, 2 * length(intra.year.timesteps))])
write.csv( file = "Data/SimulatedBiomass.csv", rowSums(Biomass))

# Output a recruitment survey with a bit of noise
#write.csv( file = "Data/SimulatedFebruaryRecruitmentSurvey.csv", Recruitment[seq(2, length(years) * 12,12),1])
write.csv( file = "Data/SimulatedFebruaryRecruitmentSurvey.csv", Recruitment[seq(2, length(years) * 12,12),1] * runif( length(seq(2, length(years) * 12,12)), min = 0.8, max = 1.2))


### NOTE its a 4 columns format, the last one representing non-targeted effort was set to 0

write.table(file = "Data/SimData3.txt", cbind(seq(1, 9*52), rowSums(noisy.Catch)[seq(53, 10*52)], Effort[seq(53, 10*52),1], 0), row.names = FALSE, col.names = FALSE)

# last 4 years
write.table(file = "Data/SimData4.txt", data.frame(seq(1, 4*52), "SimYear", rep(years, each = length(intra.year.timesteps))[seq(length(rowSums(noisy.Catch)) - 4*52 + 1, length(rowSums(noisy.Catch)))], rep(1:52, 4), rowSums(noisy.Catch)[seq(length(rowSums(noisy.Catch)) - 4*52 + 1, length(rowSums(noisy.Catch)))], Effort[seq(length(rowSums(noisy.Catch)) - 4*52 + 1, length(rowSums(noisy.Catch))),1], 0), row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(file = "Data/SimData4bis.txt", cbind(seq(1, 4*52), rowSums(Catch)[seq(length(rowSums(Catch)) - 4*52 + 1, length(rowSums(Catch)))], Effort[seq(length(rowSums(Catch)) - 4*52 + 1, length(rowSums(Catch))),1], 0), row.names = FALSE, col.names = FALSE)

# last 23 years
write.table(file = "Data/SimData5.txt", cbind(seq(1, 23*52), rowSums(noisy.Catch)[seq(length(rowSums(noisy.Catch)) - 23*52 + 1, length(rowSums(noisy.Catch)))], Effort[seq(length(rowSums(noisy.Catch)) - 23*52 + 1, length(rowSums(noisy.Catch))),1], 0), row.names = FALSE, col.names = FALSE)

write.table(file = "Data/SimData6.txt", cbind(seq(1, 23*52), rowSums(Catch)[seq(length(rowSums(Catch)) - 23*52 + 1, length(rowSums(Catch)))], Effort[seq(length(rowSums(Catch)) - 23*52 + 1, length(rowSums(Catch))),1],0), row.names = FALSE, col.names = FALSE)

