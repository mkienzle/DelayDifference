# CREATED  14 November 2013
# MODIFIED 13 June 2015

# PURPOSE display the results of the simulation to test the delay difference model option 4 (written in C++) using data aggregated on weekly timesteps

# BACKGROUND since the number of simulation

# Number of simulations passed as an argument to the script
args <- commandArgs(trailingOnly = TRUE)
nb.sim <- as.numeric(args[2])

# cases where MINUIT converged to boundaries or stopped prematuraly
failed.minimization <- unique(read.table(paste("Results/", args[1],"/FailedFit.txt", sep=""))$V1)

sim.results <- data.frame(matrix(ncol = 29, nrow = nb.sim))

dimnames(sim.results)[[2]] <- list("Sim.targeted.q", "Est.targeted.q", "sd.Est.targeted.q", 
			      "Sim.vm.mean1", "Est.vm.mean1", "sd.Est.vm.mean1", 
			      "Sim.vm.mean2", "Est.vm.mean2", "sd.Est.vm.mean2", 
			      "Sim.vm.mean3", "Est.vm.mean3", "sd.Est.vm.mean3", 
			      "Sim.vm.mean4", "Est.vm.mean4", "sd.Est.vm.mean4", 
			      "Sim.vm.sigma", "Est.vm.sigma", "sd.Est.vm.sigma",
                              "Sim.B1", "Est.B1", "sd.Est.B1",
			      "Sim.B2", "Est.B2", "sd.Est.B2",
			      "Sim.NatMort", "Est.NatMort", "sd.Est.NatMort","Sigma","SigmaError")

for(i in setdiff(1:nb.sim, failed.minimization)){
sim.data <- read.table(paste("Results/", args[1], "/", i, "/SimData4.txt", sep = ""))

sim.par <- read.table(paste("Results/", args[1], "/", i, "/SimPar.txt", sep = ""))
sim.results$Sim.targeted.q[i] <- sim.par$V1[1]
sim.results$Sim.NatMort[i] <- sim.par$V1[2]
sim.results$Sim.vm.mean1[i] <- sim.par$V1[54]
sim.results$Sim.vm.mean2[i] <- sim.par$V1[55]
sim.results$Sim.vm.mean3[i] <- sim.par$V1[56]
sim.results$Sim.vm.mean4[i] <- sim.par$V1[57]

sim.results$Sim.vm.sigma[i] <- sim.par$V1[58]

par.est <- read.csv(paste("Results/", args[1], "/", i, "/ParameterEstimates.txt", sep = ""), header = FALSE)

sim.results$Est.NatMort[i] <- par.est$V2[1]
sim.results$sd.Est.NatMort[i] <- par.est$V3[1]

sim.results$Est.targeted.q[i] <- par.est$V2[2]
sim.results$sd.Est.targeted.q[i] <- par.est$V3[2]


#sim.results$Est.vm.mean1[i] <- par.est$V2[12]
sim.results$Est.vm.mean1[i] <- ifelse(par.est$V2[12]<0,-1,1) * abs(par.est$V2[12]) %% pi
sim.results$sd.Est.vm.mean1[i] <- par.est$V3[12]

sim.results$Est.vm.mean2[i] <- ifelse(par.est$V2[13]<0,-1,1) * abs(par.est$V2[13]) %% pi
sim.results$sd.Est.vm.mean2[i] <- par.est$V3[13]

sim.results$Est.vm.mean3[i] <- ifelse(par.est$V2[14]<0,-1,1) * abs(par.est$V2[14]) %% pi
sim.results$sd.Est.vm.mean3[i] <- par.est$V3[14]

sim.results$Est.vm.mean4[i] <- ifelse(par.est$V2[15]<0,-1,1) * abs(par.est$V2[15]) %% pi
sim.results$sd.Est.vm.mean4[i] <- par.est$V3[15]

sim.results$Est.vm.sigma[i] <- par.est$V2[7]
sim.results$sd.Est.vm.sigma[i] <- par.est$V3[7]

sim.results$Sigma[i] <- par.est$V2[4]
sim.results$SigmaError[i] <- par.est$V3[4]


##### Biomass
print(i)
sim.biomass <- read.csv(paste("Results/", args[1], "/", i, "/SimulatedBiomass.csv", sep = ""), header = TRUE)
sim.results$Sim.B1[i] <- sim.biomass[seq(nrow(sim.biomass) - 4*52 + 1, nrow(sim.biomass) - 4*52 + 2),"x"][1]
sim.results$Est.B1[i] <- par.est$V2[5] * 1e5
sim.results$sd.Est.B1[i] <- par.est$V3[5] * 1e5

sim.results$Sim.B2[i] <- sim.biomass[seq(nrow(sim.biomass) - 4*52 + 1, nrow(sim.biomass) - 4*52 + 2),"x"][2]
sim.results$Est.B2[i] <- par.est$V2[6] * 1e5
sim.results$sd.Est.B2[i] <- par.est$V3[6] * 1e5

}
print("... Finished to load the data for plotting")

#### Output summary
sim.summary <- cbind("WhiteNoiseLevel"=args[1], sim.results)
write.csv(sim.summary, file = paste("Results/Tables/", args[1], "SimulationsSummary.csv", sep = ""))

##### Plot
#par(mfrow = c(3,2))

# all results
postscript(file = paste("Results/Graphics/", args[1], "-EstimateVsSimulate-NaturalMortality.ps", sep =""))
#with(sim.results, plot(Sim.NatMort, Est.NatMort, xlab = "Simulated natural mortality", ylab = "Estimated natural mortality (+- 2 SD)", main = paste("white noise level ", args[1], " %"), pch = 19, type = "n", las = 1, cex = 0.2))
with(sim.results, plot(Sim.NatMort, Est.NatMort, xlab = "Simulated natural mortality", ylab = "Estimated natural mortality (+- 2 SD)", pch = 19, type = "n", las = 1, cex = 0.2))
for(i in 1:nb.sim) with(sim.results, text( Sim.NatMort, Est.NatMort, pch = i))
abline(0,1, lty = 3, col = "black", lwd = 1.5)
with(sim.results, segments(Sim.NatMort, Est.NatMort - 2 * sd.Est.NatMort, Sim.NatMort, Est.NatMort + 2 * sd.Est.NatMort))
dev.off()

# a sample of a 100
print(nrow(sim.results))
small.sample <- sample(seq(1, nrow(sim.results))[complete.cases(sim.results)], 100)

postscript(file = paste("Results/Graphics/", args[1], "-EstimateVsSimulate-NaturalMortality-Subset.ps", sep =""))
par(fg = "black", col.main = "black", col.axis = "black", col.lab = "black", cex = 1.5, cex.lab = 1.8, cex.axis = 1.5, cex.main = 2.5, mai = c(1.5, 1.8, 1.02, 0.2), mgp = c(3.8, 1, 0))
with(sim.results[small.sample,], plot(Sim.NatMort, Est.NatMort, xlab = "Simulated", ylab = "Estimated (+- 2 SD)", main = "Natural mortality", pch = 19, type = "n", las = 1, cex = 0.2, col = "black"))
with(sim.results[small.sample,], points( Sim.NatMort, Est.NatMort, pch = 19, col = "black"))
with(sim.results[small.sample,], segments(Sim.NatMort, Est.NatMort - 2 * sd.Est.NatMort, Sim.NatMort, Est.NatMort + 2 * sd.Est.NatMort, col = "black"))
abline(0,1, lty = 3, col = "black", lwd = 1.5)
dev.off()

# all results
postscript(file = paste("Results/Graphics/", args[1], "-EstimateVsSimulate-TargetedCatchability.ps", sep = ""))
with(sim.results, plot(Sim.targeted.q, Est.targeted.q, main = paste("white noise level ", args[1], " %"), xlab = "Simulated targeted catchability", ylab = "Estimated targeted catchability (+- 2 SD)", pch = 19, type = "n", las = 1, cex = 0.2))
for(i in 1:nb.sim) with(sim.results, text( Sim.targeted.q, Est.targeted.q, pch = i))
abline(0,1, lty = 3, col = "black", lwd = 1.5)
with(sim.results, segments(Sim.targeted.q, Est.targeted.q - 2 * sd.Est.targeted.q, Sim.targeted.q, Est.targeted.q + 2 * sd.Est.targeted.q))
dev.off()

# subset
postscript(file = paste("Results/Graphics/", args[1], "-EstimateVsSimulate-TargetedCatchability-Subset.ps", sep = ""))
par(fg = "black", col.main = "black", col.axis = "black", col.lab = "black", cex = 1.5, cex.lab = 1.8, cex.axis = 1.5, cex.main = 2.5, mai = c(1.5, 1.8, 1.02, 0.2), mgp = c(3.8, 1, 0))
with(sim.results[small.sample,], plot(Sim.targeted.q, Est.targeted.q, xlab = "Simulated", ylab = "Estimated (+- 2 SD)", main = "Catchability", pch = 19, type = "n", las = 1, cex = 0.2, col = "black"))
with(sim.results[small.sample,], points( Sim.targeted.q, Est.targeted.q, pch = 19, col = "black"))
with(sim.results[small.sample,], segments(Sim.targeted.q, Est.targeted.q - 2 * sd.Est.targeted.q, Sim.targeted.q, Est.targeted.q + 2 * sd.Est.targeted.q, col = "black"))
abline(0,1, lty = 3, col = "black", lwd = 1.5)
dev.off()

# all results
postscript(file = paste("Results/Graphics/", args[1], "-EstimateVsSimulate-VonMisesMean.ps", sep = ""))
with(sim.results, plot(Sim.vm.mean1, Est.vm.mean1, main = paste("white noise level ", args[1], " %"), xlab = "Simulated von Mises mean", ylab = "Estimated von Mises mean (+- 2 SD)", pch = 19, type = "p", las = 1, cex = 0.2))
with(sim.results, segments(Sim.vm.mean1, Est.vm.mean1 - 2 * sd.Est.vm.mean1, Sim.vm.mean1, Est.vm.mean1 + 2 * sd.Est.vm.mean1))

with(sim.results, points(Sim.vm.mean2, Est.vm.mean2, pch = 19))
with(sim.results, segments(Sim.vm.mean2, Est.vm.mean2 - 2 * sd.Est.vm.mean2, Sim.vm.mean2, Est.vm.mean2 + 2 * sd.Est.vm.mean2))

with(sim.results, points(Sim.vm.mean3, Est.vm.mean3, pch = 19))
with(sim.results, segments(Sim.vm.mean3, Est.vm.mean3 - 2 * sd.Est.vm.mean3, Sim.vm.mean3, Est.vm.mean3 + 2 * sd.Est.vm.mean3))

with(sim.results, points(Sim.vm.mean4, Est.vm.mean4, pch = 19))
with(sim.results, segments(Sim.vm.mean4, Est.vm.mean4 - 2 * sd.Est.vm.mean4, Sim.vm.mean4, Est.vm.mean4 + 2 * sd.Est.vm.mean4))


abline(0,1, lty = 3, col = "black", lwd = 1.5)
dev.off()


# a subset
postscript(file = paste("Results/Graphics/", args[1], "-EstimateVsSimulate-VonMisesMean-Subset.ps", sep = ""))
par(fg = "black", col.main = "black", col.axis = "black", col.lab = "black", cex = 1.5, cex.lab = 1.8, cex.axis = 1.5, cex.main = 2.5, mai = c(1.5, 1.8, 1.02, 0.2), mgp = c(3.8, 1, 0))
with(sim.results[small.sample,], plot(Sim.vm.mean1, Est.vm.mean1, xlab = "Simulated", ylab = "Estimated (+- 2 SD)", main = "von Mises mean", pch = 19, type = "p", las = 1, cex = 0.2, col = "black"))
with(sim.results[small.sample,], segments(Sim.vm.mean1, Est.vm.mean1 - 2 * sd.Est.vm.mean1, Sim.vm.mean1, Est.vm.mean1 + 2 * sd.Est.vm.mean1, col = "black"))

with(sim.results[small.sample,], points(Sim.vm.mean2, Est.vm.mean2, pch = 19, col = "black"))
with(sim.results[small.sample,], segments(Sim.vm.mean2, Est.vm.mean2 - 2 * sd.Est.vm.mean2, Sim.vm.mean2, Est.vm.mean2 + 2 * sd.Est.vm.mean2, col = "black"))

with(sim.results[small.sample,], points(Sim.vm.mean3, Est.vm.mean3, pch = 19, col = "black"))
with(sim.results[small.sample,], segments(Sim.vm.mean3, Est.vm.mean3 - 2 * sd.Est.vm.mean3, Sim.vm.mean3, Est.vm.mean3 + 2 * sd.Est.vm.mean3, col = "black"))

with(sim.results[small.sample,], points(Sim.vm.mean4, Est.vm.mean4, pch = 19, col = "black"))
with(sim.results[small.sample,], segments(Sim.vm.mean4, Est.vm.mean4 - 2 * sd.Est.vm.mean4, Sim.vm.mean4, Est.vm.mean4 + 2 * sd.Est.vm.mean4, col = "black"))

abline(0,1, lty = 3, col = "black", lwd = 1.5)
dev.off()

# all results
postscript(file = paste("Results/Graphics/", args[1], "-EstimateVsSimulate-VonMisesStdDeviation.ps", sep=""))
with(sim.results, plot(Sim.vm.sigma, Est.vm.sigma, main = paste("white noise level ", args[1], " %"), xlab = "Simulated von Mises SD ", ylab = "Estimated von Mises SD (+- 2 SD)", pch = 19, type = "p", las = 1, xlim = c(0,60), ylim = c(0,60), cex = 0.2))
abline(0,1, lty = 3, col = "black", lwd = 1.5)
with(sim.results, segments(Sim.vm.sigma, Est.vm.sigma - 2 * sd.Est.vm.sigma, Sim.vm.sigma, Est.vm.sigma + 2 * sd.Est.vm.sigma))
dev.off()

# a subset
postscript(file = paste("Results/Graphics/", args[1], "-EstimateVsSimulate-VonMisesStdDeviation-Subset.ps", sep=""))
par(fg = "black", col.main = "black", col.axis = "black", col.lab = "black", cex = 1.5, cex.lab = 1.8, cex.axis = 1.5, cex.main = 2.5, mai = c(1.5, 1.8, 1.02, 0.2), mgp = c(3.8, 1, 0))
with(sim.results[small.sample,], plot(Sim.vm.sigma, Est.vm.sigma, xlab = "Simulated", ylab = "Estimated (+- 2 SD)", main = "von Mises SD", pch = 19, type = "p", las = 1, xlim = c(0,60), ylim = c(0,60), col = "black"))
with(sim.results[small.sample,], segments(Sim.vm.sigma, Est.vm.sigma - 2 * sd.Est.vm.sigma, Sim.vm.sigma, Est.vm.sigma + 2 * sd.Est.vm.sigma, col = "black"))
abline(0,1, lty = 3, col = "black", lwd = 1.5)
dev.off()

# all results
postscript(file = paste("Results/Graphics/", args[1], "-EstimateVsSimulate-Biomasses.ps", sep = ""))
with(sim.results, plot(1e-6 * Sim.B1, 1e-6 * Est.B1, main = paste("white noise level ", args[1], " %"), xlab = "Simulated initial biomasses (x 1000 tonnes)", ylab = "Estimated initial biomasses (+- 2 SD) (x 1000 tonnes)", pch = 19, type = "p", las = 1, cex = 0.2))
abline(0,1, lty = 3, col = "black", lwd = 1.5)
with(sim.results, segments(1e-6 * Sim.B1, 1e-6 * (Est.B1 - 2 * sd.Est.B1), 1e-6 * Sim.B1, 1e-6 * (Est.B1 + 2 * sd.Est.B1)))

with(sim.results, points(1e-6 * Sim.B2, 1e-6 * Est.B2, pch = 19, cex = 0.2))
with(sim.results, segments(1e-6 * Sim.B2, 1e-6 * (Est.B2 - 2 * sd.Est.B2), 1e-6 * Sim.B2, 1e-6 * (Est.B2 + 2 * sd.Est.B2)))
dev.off()

# a subset
postscript(file = paste("Results/Graphics/", args[1], "-EstimateVsSimulate-Biomasses-Subset.ps", sep = ""))
par(fg = "black", col.main = "black", col.axis = "black", col.lab = "black", cex = 1.5, cex.lab = 1.8, cex.axis = 1.5, cex.main = 2.5, mai = c(1.5, 1.8, 1.02, 0.2), mgp = c(3.8, 1, 0))
with(sim.results[small.sample,], plot(1e-6 * Sim.B1, 1e-6 * Est.B1, xlab = "Simulated (x 1000 tonnes)", ylab = "Estimated (+- 2 SD) (x 1000 tonnes)", main = "Initial biomasses", pch = 19, type = "p", las = 1, col = "black"))

with(sim.results[small.sample,], segments(1e-6 * Sim.B1, 1e-6 * (Est.B1 - 2 * sd.Est.B1), 1e-6 * Sim.B1, 1e-6 * (Est.B1 + 2 * sd.Est.B1), col = "black"))

with(sim.results[small.sample,], points(1e-6 * Sim.B2, 1e-6 * Est.B2, pch = 19, col = "black"))
with(sim.results[small.sample,], segments(1e-6 * Sim.B2, 1e-6 * (Est.B2 - 2 * sd.Est.B2), 1e-6 * Sim.B2, 1e-6 * (Est.B2 + 2 * sd.Est.B2), col = "black"))
abline(0,1, lty = 3, col = "black", lwd = 1.5)
dev.off()
