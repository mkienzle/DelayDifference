# CREATED  14 November 2013
# MODIFIED  9 December 2013

# PURPOSE display the results of the simulation test of the weekly delay difference model written in C++


nb.sim <- 100

# cases where MINUIT converged to boundaries or stopped prematuraly
failed.minimization <- c(14, 24, 25, 26, 31, 37, 48, 57, 58, 61, 66, 73, 80, 88, 92, 93)

sim.results <- data.frame(matrix(ncol = 18, nrow = nb.sim))

dimnames(sim.results)[[2]] <- list("Sim.targeted.q", "Est.targeted.q", "sd.Est.targeted.q", 
			      "Sim.vm.mean", "Est.vm.mean", "sd.Est.vm.mean", 
			      "Sim.vm.sigma", "Est.vm.sigma", "sd.Est.vm.sigma",
                              "Sim.B1", "Est.B1", "sd.Est.B1",
			      "Sim.B2", "Est.B2", "sd.Est.B2",
			      "Sim.NatMort", "Est.NatMort", "sd.Est.NatMort")

for(i in setdiff(1:nb.sim, failed.minimization)){
sim.data <- read.table(paste("Results/C++Version2/", i, "/SimData4.txt", sep = ""))

sim.par <- read.table(paste("Results/C++Version2/", i, "/SimPar.txt", sep = ""))
sim.results$Sim.targeted.q[i] <- sim.par$V1[1]
sim.results$Sim.NatMort[i] <- sim.par$V1[2]
sim.results$Sim.vm.mean[i] <- sim.par$V1[3]
sim.results$Sim.vm.sigma[i] <- sim.par$V1[4]

par.est <- read.csv(paste("Results/C++Version2/", i, "/ParameterEstimates.txt", sep = ""), header = FALSE)

sim.results$Est.NatMort[i] <- par.est$V2[1]
sim.results$sd.Est.NatMort[i] <- par.est$V3[1]

sim.results$Est.targeted.q[i] <- par.est$V2[2]
sim.results$sd.Est.targeted.q[i] <- par.est$V3[2]


sim.results$Est.vm.mean[i] <- par.est$V2[7]
sim.results$sd.Est.vm.mean[i] <- par.est$V3[7]

sim.results$Est.vm.sigma[i] <- par.est$V2[8]
sim.results$sd.Est.vm.sigma[i] <- par.est$V3[8]

##### Biomass
print(i)
sim.biomass <- read.csv(paste("Results/C++Version2/", i, "/SimulatedBiomass.csv", sep = ""), header = TRUE)
sim.results$Sim.B1[i] <- sim.biomass[seq(nrow(sim.biomass) - 4*52 + 1, nrow(sim.biomass) - 4*52 + 2),"x"][1]
sim.results$Est.B1[i] <- par.est$V2[5] * 1e5
sim.results$sd.Est.B1[i] <- par.est$V3[5] * 1e5

sim.results$Sim.B2[i] <- sim.biomass[seq(nrow(sim.biomass) - 4*52 + 1, nrow(sim.biomass) - 4*52 + 2),"x"][2]
sim.results$Est.B2[i] <- par.est$V2[6] * 1e5
sim.results$sd.Est.B2[i] <- par.est$V3[6] * 1e5

}

# Plot
par(mfrow = c(3,2))
 
with(sim.results, plot(Sim.NatMort, Est.NatMort, main = "Natural mortality", xlab = "Simulated", ylab = "Estimated (+- 2 SD)", pch = 19, type = "n", las = 1, cex = 0.2))
for(i in 1:nb.sim) with(sim.results, text( Sim.NatMort, Est.NatMort, pch = i))
abline(0,1, lty = 3)
with(sim.results, segments(Sim.NatMort, Est.NatMort - 2 * sd.Est.NatMort, Sim.NatMort, Est.NatMort + 2 * sd.Est.NatMort))

with(sim.results, plot(Sim.targeted.q, Est.targeted.q, main = "Targeted catchability", xlab = "Simulated", ylab = "Estimated (+- 2 SD)", pch = 19, type = "n", las = 1, cex = 0.2))
for(i in 1:nb.sim) with(sim.results, text( Sim.targeted.q, Est.targeted.q, pch = i))
abline(0,1, lty = 3)
with(sim.results, segments(Sim.targeted.q, Est.targeted.q - 2 * sd.Est.targeted.q, Sim.targeted.q, Est.targeted.q + 2 * sd.Est.targeted.q))

with(sim.results, plot(Sim.vm.mean, Est.vm.mean, main = "Von mises distribution's mean", xlab = "Simulated", ylab = "Estimated (+- 2 SD)", pch = 19, type = "p", las = 1, cex = 0.2))
abline(0,1, lty = 3)
with(sim.results, segments(Sim.vm.mean, Est.vm.mean - 2 * sd.Est.vm.mean, Sim.vm.mean, Est.vm.mean + 2 * sd.Est.vm.mean))

with(sim.results, plot(Sim.vm.sigma, Est.vm.sigma, main = "Von mises distribution's sigma", xlab = "Simulated", ylab = "Estimated (+- 2 SD)", pch = 19, type = "p", las = 1, xlim = c(0,80), ylim = c(0,80), cex = 0.2))
abline(0,1, lty = 3)
with(sim.results, segments(Sim.vm.sigma, Est.vm.sigma - 2 * sd.Est.vm.sigma, Sim.vm.sigma, Est.vm.sigma + 2 * sd.Est.vm.sigma))

with(sim.results, plot(Sim.B1, Est.B1, main = "Biomass(1) and biomass(2)", xlab = "Simulated", ylab = "Estimated (+- 2 SD)", pch = 19, type = "p", las = 1, cex = 0.2))
abline(0,1, lty = 3)
with(sim.results, segments(Sim.B1, Est.B1 - 2 * sd.Est.B1, Sim.B1, Est.B1 + 2 * sd.Est.B1))

with(sim.results, points(Sim.B2, Est.B2, pch = 19, cex = 0.2))
with(sim.results, segments(Sim.B2, Est.B2 - 2 * sd.Est.B2, Sim.B2, Est.B2 + 2 * sd.Est.B2))

dev.print(device=postscript)
