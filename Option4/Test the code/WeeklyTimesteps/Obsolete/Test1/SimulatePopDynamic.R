# CREATED  11 May 2012
# MODIFIED 20 Aug 2013

# PURPOSE simulate fishery data to test the weekly delay difference model

# BACKGROUND here we are using fixed parameter to check that all calculations are correct

# Useful libraries
library(chron)

##### Parameters
years <- seq(1956,2010)
intra.year.timesteps <- paste("Week", 1:52, sep = "")
nb.age.groups <- seq(13,52 * 3) # Age-groups from 13 weeks old to 3 years old

#random.catchability <- runif(1, min = 2, max = 8)
#catchability.q <- random.catchability * 1e-5
catchability.q <- 5 * 1e-5
write.table(paste(" ", round(catchability.q * 1e5,2), " "), file = "Results/SimulatedParameters.txt", append = T, eol="\n", col.names = FALSE, row.names = FALSE, quote = FALSE)
#nat.mortality <- 0
nat.mortality <- 0.195 * 12 / 52 # per week

# Proportion of spawners in each month (fixed to 60%)
prop.spawners <- 0.6

recruitment.by.week <- rep(NA, length(intra.year.timesteps));
rec.pat <- c(rep(0,12), 1, rep(0,39))

recruitment.by.week[1:52] <- rec.pat / sum(rec.pat) # Must sum to 1

total.recruitment.in.year <- 1e9 * c(rep(1, length(years)-10), seq(3,0.2, length = 5), seq(0.2, 1, length = 5))

# Weight at age according to Schnute's model
# round(Eq114.bis(Schnute.results.for.TigerPrawn.2par$par, seq(13,3*52), k=rec.age),2)
weight.at.age <- c(7.66,  9.54, 11.36, 13.12, 14.82, 16.47, 18.07, 19.61, 21.10, 22.55, 23.95, 25.30,
26.61, 27.88, 29.10, 30.29, 31.44, 32.55, 33.62, 34.66, 35.67, 36.64, 37.58, 38.49,
39.38, 40.23, 41.05, 41.85, 42.63, 43.37, 44.10, 44.80, 45.48, 46.13, 46.77, 47.38,
47.98, 48.55, 49.11, 49.65, 50.17, 50.67, 51.16, 51.63, 52.09, 52.53, 52.96, 53.37,
53.77, 54.16, 54.53, 54.89, 55.25, 55.58, 55.91, 56.23, 56.54, 56.84, 57.12, 57.40,
57.67, 57.93, 58.19, 58.43, 58.67, 58.90, 59.12, 59.33, 59.54, 59.74, 59.93, 60.12,
60.30, 60.48, 60.65, 60.81, 60.97, 61.13, 61.28, 61.42, 61.56, 61.69, 61.83, 61.95,
62.07, 62.19, 62.31, 62.42, 62.53, 62.63, 62.73, 62.83, 62.92, 63.01, 63.10, 63.19,
63.27, 63.35, 63.42, 63.50, 63.57, 63.64, 63.71, 63.77, 63.84, 63.90, 63.96, 64.02,
64.07, 64.13, 64.18, 64.23, 64.28, 64.32, 64.37, 64.41, 64.46, 64.50, 64.54, 64.58,
64.61, 64.65, 64.68, 64.72, 64.75, 64.78, 64.81, 64.84, 64.87, 64.90, 64.93, 64.95,
64.98, 65.00, 65.03, 65.05, 65.07, 65.09, 65.11, 65.13, 65.15, 65.17, 65.19, 65.21)

##### Data: use Moreton Tiger Prawn data as the basis of simulated datasets

# Moreton Bay brown tiger prawn data
dataset <- read.csv("Data/MBTigerCatchAndEffortByBiologicalYearAndWeek.csv") 

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
for(i in seq(1,length(years) * length(intra.year.timesteps) - 1))
{

  # Apply mortality
  Nb.at.age[i+1, seq(2, length(nb.age.groups))] <- Nb.at.age[i, seq(1, length(nb.age.groups) - 1)] * exp( - (Nat.mortality + Fish.mortality)[i, seq(1, length(nb.age.groups) - 1)])
  # Add recruitment
  Nb.at.age[i+1,] <- Nb.at.age[i+1,] + Recruitment[i+1,]
}

# Calculate biomass in kg
Biomass <- Nb.at.age * 1e-3 * Weight.at.age

# Spawning stock biomass in kg
SSB <- tapply(rowSums(Nb.at.age * 1e-3 * Weight.at.age * prop.spawners * 0.5), rep(years, each = 52), sum)

# Calculate catch in kg
Catch <- Fish.mortality / (Fish.mortality + Nat.mortality) * Nb.at.age * 1e-3 * Weight.at.age * (1 - exp(-(Fish.mortality + Nat.mortality)))

# Monthly prawn number
par(mfrow = c(3,1))

# Plot biomass
plot(rowSums(1e-3 * Biomass), type = "b", xlab = "Time (weeks)",
  ylab = "Biomass (tons)", axes = F, pch = 19)
axis(1, at = seq(1,dim(Nb.at.age)[1])[seq(1, dim(Nb.at.age)[1],6)], label = dimnames(Nb.at.age)[[1]][seq(1, dim(Nb.at.age)[1],6)])
axis(2, las = 1)
box()
abline(v = seq(1,3e3,length(intra.year.timesteps)))

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

# Output the data from 1958 to 2010

# Add a bit of random variability
noisy.Catch <- Catch * matrix( runif(length(years) * length(intra.year.timesteps) * length(nb.age.groups), min = 0.8, max = 1.2), nrow = length(years) * length(intra.year.timesteps), ncol = length(nb.age.groups))

noisy.Effort <- Effort * matrix( runif(length(years) * length(intra.year.timesteps) * length(nb.age.groups), min = 0.8, max = 1.2), nrow = length(years) * length(intra.year.timesteps), ncol = length(nb.age.groups))

#write.csv( file = "Data/SimulatedCPUE.csv", rowSums(catchability.q * Biomass)[-seq(1,2 * length(intra.year.timesteps))])
write.csv( file = "Data/SimulatedCPUE.csv", rowSums(noisy.Catch)[-seq(1, 2 * length(intra.year.timesteps))] / noisy.Effort[-seq(1,2 * length(intra.year.timesteps)), 1])
#write.csv( file = "Data/SimulatedEffort.csv", noisy.Effort[-seq(1, 2 * length(intra.year.timesteps)), 1])
write.csv( file = "Data/SimulatedEffort.csv", Effort[-seq(1, 2 * length(intra.year.timesteps)), 1])
#write.csv( file = "Data/SimulatedCatch.csv", rowSums(noisy.Catch)[-seq(1, 2 * length(intra.year.timesteps))])
write.csv( file = "Data/SimulatedCatch.csv", rowSums(Catch)[-seq(1, 2 * length(intra.year.timesteps))])
write.csv( file = "Data/SimulatedBiomass.csv", rowSums(Biomass))

# Format output for C++
tmp <- cbind(rowSums(Catch)[-seq(1, 2 * length(intra.year.timesteps))], Effort[-seq(1, 2 * length(intra.year.timesteps)), 1])
write.table( file = "Data/SimulatedCatchAndEffort.txt", cbind(seq(1,dim(tmp)[1]), round(tmp,1)), quote = FALSE, col.names = FALSE, row.names = FALSE)

tmp2 <- cbind(rowSums(noisy.Catch)[-seq(1, 2 * length(intra.year.timesteps))], noisy.Effort[-seq(1, 2 * length(intra.year.timesteps)), 1])
write.table( file = "Data/SimulatedCatchAndEffort2.txt", cbind(seq(1,dim(tmp2)[1]), round(tmp2,1)), quote = FALSE, col.names = FALSE, row.names = FALSE)

# Output a recruitment survey with a bit of noise
#write.csv( file = "Data/SimulatedFebruaryRecruitmentSurvey.csv", Recruitment[seq(2, length(years) * 12,12),1])
write.csv( file = "Data/SimulatedFebruaryRecruitmentSurvey.csv", Recruitment[seq(2, length(years) * 12,12),1] * runif( length(seq(2, length(years) * 12,12)), min = 0.8, max = 1.2))

