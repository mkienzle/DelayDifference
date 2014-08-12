# This script was used to debug the program by comparing the result of FitMonthlyDelayDifference_v0.2 Data/SimData4.txt to simulated data

sim.biomass <- read.csv("Data/SimulatedBiomass.csv")
sim.biomass <- sim.biomass[seq(nrow(sim.biomass) - 4 * 12 + 1, nrow(sim.biomass)),]

est.quantities <- read.csv("Results/EstimatedFisheriesQuantities.txt")

plot(sim.biomass[,2], est.quantities$EstimatedBiomass, xlab = "Simulated", ylab = "Estimated", pch = 19)
abline(0,1)

print(cbind(sim.biomass[,2],est.quantities$EstimatedBiomass))
