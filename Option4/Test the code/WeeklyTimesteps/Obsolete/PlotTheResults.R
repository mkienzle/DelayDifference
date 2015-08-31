# CREATED ...
# MODIFIED 19 April 2013

##### Use the output of fitting the model to CPUE
##### Load data

# Parameters of the simulation
Sim.par <- read.table("Results/SimulatedParameters.txt", header = T, row.names = 1)

# Estimates from model fitting to CPUE
tmp <- read.table("Results/ResultsByFittingCPUE.txt", header = T, row.names = 1);

#postscript( file = "Results/NegativeCorrelationBetweenCatchabilityAndBiomassEstimatesByFittingCPUE.ps")
postscript( file = "Results/CompareOutputOfDifferentObjectiveFct.ps")
par(mfrow = c(2,2), mgp = c(2.5,1,0), mar = c(5.1, 6.1, 4.1, 2.1))

# Compare (ratio of) estimates of catchability against estimates of biomasses
plot( tmp$CatchabilityEstimate/Sim.par$Catchability, tmp$SlopeOfRegressionBetweenBiomasses, pch = 19, cex = 0.5, las = 1, xlab = "Ratio of estimated over simulated catchability \n Underestimated <---|----> Overestimated", ylab = "Slope of estimated biomass = f(simulated)\n Underestimated <----|---> Overestimated", xlim = c(-1, 4.0), ylim = c(0.5, 2), main = "Results from delay difference fitted to CPUE");

# Add errors bars on catchability estimates
segments( tmp$CatchabilityEstimate/Sim.par$Catchability - tmp$CatchabilityEstimateUncertainty/Sim.par$Catchability, tmp$SlopeOfRegressionBetweenBiomasses, tmp$CatchabilityEstimate/Sim.par$Catchability + tmp$CatchabilityEstimateUncertainty/Sim.par$Catchability, col = "lightgrey")

points(tmp$CatchabilityEstimate/Sim.par$Catchability, tmp$SlopeOfRegressionBetweenBiomasses, pch = 19, cex = 0.5)

points(1,1, pch = 3, col = "red"); abline(v=1, lty=3, col="red"); abline(h=1, lty = 3, col = "red")

# Compare estimates of catchability (+- 1 S.D.) against simulated
#with(tmp, plot(1:length(CatchabilityEstimate), CatchabilityEstimate, xlab = "simulation number", ylab = "Catchability", ylim = c(ifelse(min(CatchabilityEstimate - CatchabilityEstimateUncertainty) < 0, 1.5 * min(CatchabilityEstimate - CatchabilityEstimateUncertainty), 0.8 * min(CatchabilityEstimate - CatchabilityEstimateUncertainty)), 1.5 * max(CatchabilityEstimate + CatchabilityEstimateUncertainty)), las = 1 , pch = 19))
#with(tmp, segments(1:length(CatchabilityEstimate), CatchabilityEstimate - CatchabilityEstimateUncertainty, 1:length(CatchabilityEstimate), CatchabilityEstimate + CatchabilityEstimateUncertainty))
#with(Sim.par, points(1:length(Catchability), Catchability, pch = 19, col = "red", cex = 0.5))
#legend(1,20, pch = c(19,19), col = c("black", "red"), legend = c("Estimated (+- 1 s.d.)", "Simulated"))
#dev.off()

## ##### Use the output of fitting the model to catch
##### Load data

# Parameters of the simulation
Sim.par <- read.table("Results/SimulatedParameters.txt", header = T, row.names = 1)

# Estimates from model fitting to CPUE
tmp1 <- read.table("Results/ResultsByFittingSqrtOfCatch.txt", header = T, row.names = 1);

#postscript( file = "Results/NegativeCorrelationBetweenCatchabilityAndBiomassEstimatesByFittingCatch.ps")
#par(mfrow = c(2,2), mgp = c(2.5,1,0), mar = c(5.1, 6.1, 4.1, 2.1))

# Compare (ratio of) estimates of catchability against estimates of biomasses
plot( tmp1$CatchabilityEstimate/Sim.par$Catchability, tmp1$SlopeOfRegressionBetweenBiomasses, pch = 19, cex = 0.5, las = 1, xlab = "Ratio of estimated over simulated catchability \n Underestimated <---|----> Overestimated", ylab = "Slope of estimated biomass = f(simulated)\n Underestimated <----|---> Overestimated", xlim = c(-1, 4.0), ylim = c(0.5, 2), main = "Results from delay difference fitted to sqrt of Catch");

# Add errors bars on catchability estimates
segments( tmp1$CatchabilityEstimate/Sim.par$Catchability - tmp1$CatchabilityEstimateUncertainty/Sim.par$Catchability, tmp1$SlopeOfRegressionBetweenBiomasses, tmp1$CatchabilityEstimate/Sim.par$Catchability + tmp1$CatchabilityEstimateUncertainty/Sim.par$Catchability, col = "lightgrey")

points(tmp1$CatchabilityEstimate/Sim.par$Catchability, tmp1$SlopeOfRegressionBetweenBiomasses, pch = 19, cex = 0.5)

points(1,1, pch = 3, col = "red"); abline(v=1, lty=3, col="red"); abline(h=1, lty = 3, col = "red")

# Compare estimates of catchability (+- 1 S.D.) against simulated
#with(tmp1, plot(1:length(CatchabilityEstimate), CatchabilityEstimate, xlab = "simulation number", ylab = "Catchability", ylim = c(ifelse(min(CatchabilityEstimate - CatchabilityEstimateUncertainty) < 0, 1.5 * min(CatchabilityEstimate - CatchabilityEstimateUncertainty), 0.8 * min(CatchabilityEstimate - CatchabilityEstimateUncertainty)), 1.5 * max(CatchabilityEstimate + CatchabilityEstimateUncertainty)), las = 1 , pch = 19))
#with(tmp1, segments(1:length(CatchabilityEstimate), CatchabilityEstimate - CatchabilityEstimateUncertainty, 1:length(CatchabilityEstimate), CatchabilityEstimate + CatchabilityEstimateUncertainty))
#with(Sim.par, points(1:length(Catchability), Catchability, pch = 19, col = "red", cex = 0.5))
#legend(1,20, pch = c(19,19), col = c("black", "red"), legend = c("Estimated (+- 1 s.d.)", "Simulated"))
#dev.off()


##### Use the output of fitting the model to log of catch
##### Load data

# Parameters of the simulation
Sim.par <- read.table("Results/SimulatedParameters.txt", header = T, row.names = 1)

# Estimates from model fitting to CPUE
tmp2 <- read.table("Results/ResultsByFittingLogOfCatch.txt", header = T, row.names = 1);

#postscript( file = "Results/NegativeCorrelationBetweenCatchabilityAndBiomassEstimatesByFittingLogOfCatch.ps")
#par(mfrow = c(2,2), mgp = c(2.5,1,0), mar = c(5.1, 6.1, 4.1, 2.1))

# Compare (ratio of) estimates of catchability against estimates of biomasses
plot( tmp2$CatchabilityEstimate/Sim.par$Catchability, tmp2$SlopeOfRegressionBetweenBiomasses, pch = 19, cex = 0.5, las = 1, xlab = "Ratio of estimated over simulated catchability \n Underestimated <---|----> Overestimated", ylab = "Slope of estimated biomass = f(simulated)\n Underestimated <----|---> Overestimated", xlim = c(-1, 4.0), ylim = c(0.5, 2), main = "Results from delay difference fitted to log of catch");

# Add errors bars on catchability estimates
segments( tmp2$CatchabilityEstimate/Sim.par$Catchability - tmp2$CatchabilityEstimateUncertainty/Sim.par$Catchability, tmp2$SlopeOfRegressionBetweenBiomasses, tmp2$CatchabilityEstimate/Sim.par$Catchability + tmp2$CatchabilityEstimateUncertainty/Sim.par$Catchability, col = "lightgrey")

points(tmp2$CatchabilityEstimate/Sim.par$Catchability, tmp2$SlopeOfRegressionBetweenBiomasses, pch = 19, cex = 0.5)

points(1,1, pch = 3, col = "red"); abline(v=1, lty=3, col="red"); abline(h=1, lty = 3, col = "red")

## # Compare estimates of catchability (+- 1 S.D.) against simulated
## with(tmp2, plot(1:length(CatchabilityEstimate), CatchabilityEstimate, xlab = "simulation number", ylab = "Catchability", ylim = c(ifelse(min(CatchabilityEstimate - CatchabilityEstimateUncertainty) < 0, 1.5 * min(CatchabilityEstimate - CatchabilityEstimateUncertainty), 0.8 * min(CatchabilityEstimate - CatchabilityEstimateUncertainty)), 1.5 * max(CatchabilityEstimate + CatchabilityEstimateUncertainty)), las = 1 , pch = 19))
## with(tmp2, segments(1:length(CatchabilityEstimate), CatchabilityEstimate - CatchabilityEstimateUncertainty, 1:length(CatchabilityEstimate), CatchabilityEstimate + CatchabilityEstimateUncertainty))
## with(Sim.par, points(1:length(Catchability), Catchability, pch = 19, col = "red", cex = 0.5))
## legend(1,20, pch = c(19,19), col = c("black", "red"), legend = c("Estimated (+- 1 s.d.)", "Simulated"))
## dev.off()

## tmp2 <- read.table("Results/ResultsByFittingCatch.txt");

## postscript( file = "Results/NegativeCorrelationBetweenCatchabilityAndBiomassEstimatesByFittingCatch.ps")
## par(mfrow = c(2,2), mgp = c(2.5,1,0))

## # Compare (ratio of) estimates of catchability against estimates of biomasses
## with(tmp2, plot( V3/V2, V5, pch = 19, cex = 0.5, las = 1, xlab = "Ratio of estimated over simulated catchability \n Underestimated <---|----> Overestimated", ylab = "Slope of estimated biomass = f(simulated)\n Underestimated <----|---> Overestimated", xlim = c(-1, 4.0), ylim = c(0.5, 2)));
## with(tmp2, segments( V3/V2 - V4/V2, V5, V3/V2 + V4/V2, col = "lightgrey")) 

## with(tmp2, points(V3/V2, V5, pch = 19, cex = 0.5))

## points(1,1, pch = 3, col = "red"); abline(v=1, lty=3, col="red"); abline(h=1, lty = 3, col = "red")

## # Compare estimates of catchability (+- 1 S.D.) against simulated
## with(tmp2, plot(V1, V3, xlab = "simulation number", ylab = "Catchability", ylim = c(ifelse(min(V3-V4) < 0, 1.5 * min(V3-V4), 0.8 * min(V3-V4)), 1.5 * max(V3+V4)), las = 1 , pch = 19))
## with(tmp2, segments(V1, V3-V4, V1, V3+V4))
## with(tmp2, points(V1, V2, pch = 19, col = "red", cex = 0.5))
## legend(1,20, pch = c(19,19), col = c("black", "red"), legend = c("Estimated (+- 1 s.d.)", "Simulated"))

##### Use the output of fitting the model to catch
##### Load data

# Parameters of the simulation
Sim.par <- read.table("Results/SimulatedParameters.txt", header = T, row.names = 1)

# Estimates from model fitting to CPUE
tmp3 <- read.table("Results/ResultsByFittingCatch.txt", header = T, row.names = 1);

#postscript( file = "Results/NegativeCorrelationBetweenCatchabilityAndBiomassEstimatesByFittingLogOfCatch.ps")
#par(mfrow = c(2,2), mgp = c(2.5,1,0), mar = c(5.1, 6.1, 4.1, 2.1))

# Compare (ratio of) estimates of catchability against estimates of biomasses
plot( tmp3$CatchabilityEstimate/Sim.par$Catchability, tmp3$SlopeOfRegressionBetweenBiomasses, pch = 19, cex = 0.5, las = 1, xlab = "Ratio of estimated over simulated catchability \n Underestimated <---|----> Overestimated", ylab = "Slope of estimated biomass = f(simulated)\n Underestimated <----|---> Overestimated", xlim = c(-1, 4.0), ylim = c(0.5, 2), main = "Results from delay difference fitted to catch");

# Add errors bars on catchability estimates
segments( tmp3$CatchabilityEstimate/Sim.par$Catchability - tmp3$CatchabilityEstimateUncertainty/Sim.par$Catchability, tmp3$SlopeOfRegressionBetweenBiomasses, tmp3$CatchabilityEstimate/Sim.par$Catchability + tmp3$CatchabilityEstimateUncertainty/Sim.par$Catchability, col = "lightgrey")

points(tmp3$CatchabilityEstimate/Sim.par$Catchability, tmp3$SlopeOfRegressionBetweenBiomasses, pch = 19, cex = 0.5)

points(1,1, pch = 3, col = "red"); abline(v=1, lty=3, col="red"); abline(h=1, lty = 3, col = "red")


dev.off()
