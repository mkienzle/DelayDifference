# CREATED  11 May 2012
# MODIFIED 21 May 2014

# PURPOSE simulate fishery data to test the weekly delay difference model

# Useful libraries
library(chron)
library("RcppArmadillo")
library ("inline")
source("RcppArmadilloFastLoop.R")

##### Parameters
years <- seq(1956,2010)
intra.year.timesteps <- paste("Month", 1:12, sep = "")
nb.age.groups <- seq(5,12 * 3) # Age-groups from 5 months old to 3 years old

random.catchability <- runif(1, min = 2, max = 8)
catchability.q <- random.catchability * 1e-4
#catchability.q <- 5 * 1e-5
write.table(paste(" ", round(catchability.q * 1e4,2), " "), file = "Results/SimulatedParameters.txt", append = T, eol="\n", col.names = FALSE, row.names = FALSE, quote = FALSE)

write.table(paste("Targeted catchability,",round(catchability.q * 1e4,2)), file = "Results/SimPar.txt", eol="\n", col.names = FALSE, row.names = FALSE, quote = FALSE)
#nat.mortality <- 0
nat.mortality <- 0.195 # per month

# Proportion of spawners in each month (fixed to 60%)
prop.spawners <- 0.6

###
## Recruitment distribution
###
recruitment.by.week <- rep(NA, length(intra.year.timesteps));

# Parameters
#A <- 0
#B <- 10
A <- runif(1, min = -3.1, max = 3.1)
B <- runif(1, min = 1, max = 70)

write.table(paste("von Mises mean,", round(A,2)), append = TRUE, file = "Results/SimPar.txt", eol="\n", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(paste("von Mises sigma,", round(B,2)), append = TRUE, file = "Results/SimPar.txt", eol="\n", col.names = FALSE, row.names = FALSE, quote = FALSE)

library(circular)
# Interval of time you are interested in
boundaries <- seq(-pi, pi, length =13)

precision <- 1e3
res <- rep(NA, 12); 
for(i in 1:12) res[i] <- sum(dvonmises(seq(boundaries[i], boundaries[i+1] - 1 / precision, 1 / precision), A, B))/precision
recruitment.by.week <- res # / sum(res)

###
# Total recruitment each year
###

rand.rec <- runif(10, min = 0.1, max = 10)
total.recruitment.in.year <- 1e7 * c(rep(1, length(years)-10), rand.rec)
write.table(paste("Recruitment,", round(rand.rec,2)), append = TRUE, file = "Results/SimPar.txt", eol="\n", col.names = FALSE, row.names = FALSE, quote = FALSE)

###
# Weight at age according to Schnute's model
###
source("EstimatingRho.r")
weight.at.age <- round(Eq114.bis(Schnute.results.for.TigerPrawn.2par$par, seq(rec.age,3*12), k=rec.age),2)

##### Data: use Moreton Tiger Prawn data as the basis of simulated datasets

# Moreton Bay brown tiger prawn data
#dataset <- read.csv("Data/MBTigerCatchAndEffortByBiologicalYearAndWeek.csv") 
dataset <- read.csv("Data/MBTigerCatchAndEffortByBiologicalYearAndMonth.csv")

## Using constant effort from year to year
effort <- as.matrix(rep(subset(dataset, Biological.Year == 2004)$Total.Tiger.Effort, length(years)), ncol = 1)

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

# Varying recruitment
Recruitment[,1] <- c(outer(recruitment.by.week, total.recruitment.in.year))

##### Population dynamic
tmp <- FastLoop(as.matrix(Nb.at.age), as.matrix(Nat.mortality), as.matrix(Fish.mortality), as.matrix(Recruitment))
tmp <- as.data.frame(tmp)
dimnames(tmp) <- dimnames(Nb.at.age)
Nb.at.age <- tmp

# Calculate biomass in kg
Biomass <- Nb.at.age * 1e-3 * Weight.at.age

# Spawning stock biomass in kg
SSB <- tapply(rowSums(Nb.at.age * 1e-3 * Weight.at.age * prop.spawners * 0.5), rep(years, each = 12), sum)

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
write.table(file = "Data/SimData4.txt", data.frame(seq(1, 4*12), "SimYear", rep(years, each = length(intra.year.timesteps))[seq(length(rowSums(noisy.Catch)) - 4*12 + 1, length(rowSums(noisy.Catch)))], rep(1:12, 4), rowSums(noisy.Catch)[seq(length(rowSums(noisy.Catch)) - 4*12 + 1, length(rowSums(noisy.Catch)))], Effort[seq(length(rowSums(noisy.Catch)) - 4*12 + 1, length(rowSums(noisy.Catch))),1], 0), row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(file = "Data/SimData4bis.txt", cbind(seq(1, 4*12), rowSums(Catch)[seq(length(rowSums(Catch)) - 4*12 + 1, length(rowSums(Catch)))], Effort[seq(length(rowSums(Catch)) - 4*12 + 1, length(rowSums(Catch))),1], 0), row.names = FALSE, col.names = FALSE)

# last 23 years
write.table(file = "Data/SimData5.txt", cbind(seq(1, 23*12), rowSums(noisy.Catch)[seq(length(rowSums(noisy.Catch)) - 23*12 + 1, length(rowSums(noisy.Catch)))], Effort[seq(length(rowSums(noisy.Catch)) - 23*12 + 1, length(rowSums(noisy.Catch))),1], 0), row.names = FALSE, col.names = FALSE)

write.table(file = "Data/SimData6.txt", cbind(seq(1, 23*12), rowSums(Catch)[seq(length(rowSums(Catch)) - 23*12 + 1, length(rowSums(Catch)))], Effort[seq(length(rowSums(Catch)) - 23*12 + 1, length(rowSums(Catch))),1],0), row.names = FALSE, col.names = FALSE)
