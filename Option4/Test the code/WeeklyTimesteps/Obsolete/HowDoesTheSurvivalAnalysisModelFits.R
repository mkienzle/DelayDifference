# CREATED   21 June 2012
# MODIFIED  21 June 2012

# PURPOSE parameter estimates using survival analysis

### Load the data
# Plot the data
source("PlotSimulatedData.R")

## Load the data
sim.catch <- read.csv("Data/SimulatedCatch.csv")
sim.effort <- read.csv("Data/SimulatedEffort.csv")

# Reformat into matrix year x month
sim.catch.mat <- with(sim.catch, matrix(x, ncol = 12, byrow = T)); dimnames(sim.catch.mat) <- with(sim.catch, list( unique(substr(X,5,8)), unique(substr(X,1,3))))
sim.effort.mat <- with(sim.effort, matrix(x, ncol = 12, byrow = T)); dimnames(sim.effort.mat) <- with(sim.effort, list( unique(substr(X,5,8)), unique(substr(X,1,3))))

mat.dim <- dim(sim.catch.mat)

# Take only the last 10 years of data
#catch <- sim.catch.mat[ seq(mat.dim[1] - 10 + 1, mat.dim[1]),]
#effort <- sim.effort.mat[ seq(mat.dim[1] - 10  + 1, mat.dim[1]),]
#catch <- sim.catch.mat[ mat.dim[1],]
#effort <- sim.effort.mat[ mat.dim[1],] / 30
catch <- sim.catch[seq(636-15+1,636-3),2]
effort <- sim.effort[seq(636-15+1,636-3),2] / 30

# weight.at.age <- c(12.19267, 19.34308, 26.17824, 32.23908, 37.37359, 41.59756, 45.00507, 47.71736, 49.85630, 51.53215, 52.83919,
# 53.85537, 54.64368, 55.25434, 55.72693, 56.09248, 56.37515, 56.59373, 56.6, 56.6, 56.6, 56.6)

# nls( y ~ vbgf(x, a, b, c), data = data.frame(x = seq(3, 24), y = weight.at.age), start = list( a = 70, b = 0.2, c=0))

### The model
# Load libraries and useful functions
library(chron)
source("~/mystuff/Work/DEEDI/Moreton\ Bay\ Prawn\ Trawl\ fishery/Analysis/Scripts/Mortality\ estimates/Ccode/SurvivalLib.R")

# Suppose each year is broken down into months
time.breaks <- c(seq(0,330,30), 365)

##### Likelihood function to estimate q using catch and effort data grouped by month
# log-likelihood of catching 1 kg of prawn in a given month
LogLik.estimator.of.q <- function(q, nat.mort){

# ARGUMENTS parameter 1 is q and parameter 2 is natural mortality
# ARGUMENTS VBGF in weight is assumed to be known

# A vector to hold the proportion of catch in each month
prop <- rep(NA, length = (length(time.breaks)-1))

# Calculate the proportion of catch in each month using a truncated distribution of the biomass caught
for(i in 1:(length(time.breaks)-1)){

prop[i] <- CW(seq(time.breaks[i],time.breaks[i+1]-1), time.breaks, nat.mort, 
1e-5 * q * effort, Winf = 57 , k = 0.27/30, t0 = -0.08)  / 
CW(seq(time.breaks[1], time.breaks[length(time.breaks)] - 1), time.breaks, nat.mort, 
1e-5 * q * effort, Winf = 57, k = 0.27/30, t0 = -0.08)

}

# Remove the zeroes to calculate the log-likelihood
zeroes <- which(prop == 0)

# return log-likelihood
ifelse(length(zeroes) > 0, -sum(catch[-zeroes] * log(prop[-zeroes])), 
  -sum(catch * log(prop)))
}

# Search for maximum likelihood estimator of catchability
solution <- optimize(LogLik.estimator.of.q , interval = c(1e-2,3e1), nat.mort = 0.15/30)
print(paste("Estimated catchability is:", round(solution$minimum,2)))

# Plot data overlaid with the model
catch <- as.table(catch); names(catch) <- months(1:12 * 28)
# A vector to hold the proportion of catch in each month
prop <- rep(NA, length = length(time.breaks)-1)

# Calculate the proportion of catch in each month using a truncated distribution of the biomass caught
for(i in 1:(length(time.breaks)-1))
prop[i] <- CW(seq(time.breaks[i],time.breaks[i+1]-1), time.breaks, nat.mort = 0.15/30, 
1e-5 * solution$minimum * effort, Winf = 57, k = 0.27/30, t0 = -0.08) / 
CW(seq(time.breaks[1], time.breaks[length(time.breaks)] - 1), time.breaks, nat.mort = 0.15/30, 
1e-5 * solution$minimum * effort, Winf = 57, k = 0.27/30, t0 = -0.08)


my.width <- 1
my.space <- 0.2
barplot(catch, las = 1, xlab = "Month", ylab = "Catch", ylim = c(0, 1.2 * max(catch)))
lines( 0.5 * (my.width + my.width * my.space) + my.width * (1 + my.space) * (seq(1,12) - 1), sum(catch) * prop, type = "b", pch = 19)

