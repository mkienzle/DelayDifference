# CREATED  14 November 2013
# MODIFIED 18 May 2014

# PURPOSE display the results of the simulation test of the monthly delay difference model written in C++


nb.sim <- 100

# cases where MINUIT converged to boundaries or stopped prematuraly
failed.minimization <- c(5,8,10,29,30,32,35,39,48,54,66,67,78,80,81,84,89,92,97,98)

sim.results <- data.frame(matrix(ncol = 15, nrow = nb.sim))

dimnames(sim.results)[[2]] <- list("Sim.targeted.q", "Est.targeted.q", "sd.Est.targeted.q", 
			      "Sim.vm.mean", "Est.vm.mean", "sd.Est.vm.mean", 
			      "Sim.vm.sigma", "Est.vm.sigma", "sd.Est.vm.sigma",
                              "Sim.B1", "Est.B1", "sd.Est.B1",
			      "Sim.B2", "Est.B2", "sd.Est.B2")

for(i in setdiff(1:nb.sim, failed.minimization)){
sim.data <- read.table(paste("Results/", i, "/SimData4.txt", sep = ""))

#sim.par <- read.table(paste("Results/", i, "/SimPar.txt", sep = ""))
sim.par <- read.csv(paste("Results/", i, "/SimPar.txt", sep = ""), head = FALSE)
sim.results$Sim.targeted.q[i] <- sim.par$V2[1]
sim.results$Sim.vm.mean[i] <- sim.par$V2[2]
sim.results$Sim.vm.sigma[i] <- sim.par$V2[3]

par.est <- read.csv(paste("Results/", i, "/ParameterEstimates.txt", sep = ""), header = FALSE)
sim.results$Est.targeted.q[i] <- par.est$V2[1]
sim.results$sd.Est.targeted.q[i] <- par.est$V3[1]

sim.results$Est.vm.mean[i] <- par.est$V2[6]
sim.results$sd.Est.vm.mean[i] <- par.est$V3[6]

sim.results$Est.vm.sigma[i] <- par.est$V2[7]
sim.results$sd.Est.vm.sigma[i] <- par.est$V3[7]

##### Biomass
print(i)
sim.biomass <- read.csv(paste("Results/", i, "/SimulatedBiomass.csv", sep = ""), header = TRUE)
sim.results$Sim.B1[i] <- sim.biomass[seq(nrow(sim.biomass) - 4*12 + 1, nrow(sim.biomass) - 4*12 + 2),"x"][1]
sim.results$Est.B1[i] <- par.est$V2[4] * 1e5
sim.results$sd.Est.B1[i] <- par.est$V3[4] * 1e5

sim.results$Sim.B2[i] <- sim.biomass[seq(nrow(sim.biomass) - 4*12 + 1, nrow(sim.biomass) - 4*12 + 2),"x"][2]
sim.results$Est.B2[i] <- par.est$V2[5] * 1e5
sim.results$sd.Est.B2[i] <- par.est$V3[5] * 1e5

}

# Plot
par(mfrow = c(2,2))
 
with(sim.results, plot(Sim.targeted.q, Est.targeted.q, main = "Targeted catchability", xlab = "Simulated", ylab = "Estimated (+- 2 SD)", pch = 19, type = "n", las = 1, cex = 0.2))
for(i in 1:nb.sim) with(sim.results, text( Sim.targeted.q, Est.targeted.q, pch = i))
abline(0,1, lty = 3)
with(sim.results, segments(Sim.targeted.q, Est.targeted.q - 2 * sd.Est.targeted.q, Sim.targeted.q, Est.targeted.q + 2 * sd.Est.targeted.q))

with(sim.results, plot(Sim.vm.mean, Est.vm.mean, main = "Von mises distribution's mean", xlab = "Simulated", ylab = "Estimated (+- 2 SD)", pch = 19, type = "p", las = 1, cex = 0.2, xlim = c(-pi, pi), ylim = c(-pi,pi)))
abline(0,1, lty = 3)
with(sim.results, segments(Sim.vm.mean, Est.vm.mean - 2 * sd.Est.vm.mean, Sim.vm.mean, Est.vm.mean + 2 * sd.Est.vm.mean))

with(sim.results, plot(Sim.vm.sigma, Est.vm.sigma, main = "Von mises distribution's sigma", xlab = "Simulated", ylab = "Estimated (+- 2 SD)", pch = 19, type = "p", las = 1, xlim = c(0,20), ylim = c(0,20), cex = 0.2))
abline(0,1, lty = 3)
with(sim.results, segments(Sim.vm.sigma, Est.vm.sigma - 2 * sd.Est.vm.sigma, Sim.vm.sigma, Est.vm.sigma + 2 * sd.Est.vm.sigma))

with(sim.results, plot(1e-6 * Sim.B1, 1e-6 * Est.B1, main = "Biomass(1) and biomass(2)", xlab = "Simulated (millions)", ylab = "Estimated (+- 2 SD) ( millions )", pch = 19, type = "p", las = 1, cex = 0.2))
abline(0,1, lty = 3)
with(sim.results, segments(1e-6 * Sim.B1, 1e-6 * (Est.B1 - 2 * sd.Est.B1), 1e-6 * Sim.B1, 1e-6 * (Est.B1 + 2 * sd.Est.B1)))

with(sim.results, points(1e-6 * Sim.B2, 1e-6 * Est.B2, pch = 19, cex = 0.2))
with(sim.results, segments(1e-6 * Sim.B2, 1e-6 * (Est.B2 - 2 * sd.Est.B2), 1e-6 * Sim.B2, 1e-6 * (Est.B2 + 2 * sd.Est.B2)))

dev.print(device=postscript)
