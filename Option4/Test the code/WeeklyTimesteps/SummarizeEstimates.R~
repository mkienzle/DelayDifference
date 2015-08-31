# CREATED   10 June 2015
# MODIFIED  10 June 2015

# PURPOSE summary estimates from the delay difference model coded in C++

# Load the data
error.range <- seq(10, 40, 10)

for(i in 1:length(error.range)){

ifelse(i == 1, tmp <- read.csv(paste("Results/Tables/WhiteNoiseLevel", error.range[i], "SimulationsSummary.csv", sep = "")),
tmp <- rbind(tmp, read.csv(paste("Results/Tables/WhiteNoiseLevel", error.range[i], "SimulationsSummary.csv", sep = ""))))

}

##### Natural mortality

# discrepancy in units of standard deviation
tmp.NatMort <- with(tmp, aggregate(abs(Est.NatMort - Sim.NatMort)/sd.Est.NatMort, by = list(WhiteNoiseLevel), function(x) quantile(x, c(seq(0.1,0.9,0.1), 0.95, 0.975, 0.99, 1), na.rm = TRUE)))

library(xtable)
tmp.NatMort.xtable <- as.matrix(tmp.NatMort[,-1])
dimnames(tmp.NatMort.xtable)[[1]] <- paste("White noise level ", substr(tmp.NatMort[,1], 16, 20), "%")
tmp.NatMort.xtable <- xtable(tmp.NatMort.xtable, caption = "Quantiles of the distribution of standardized natural mortality estimates ($|\\hat{M} - M|/\\sigma_{\\hat{M}}$) for varying levels of random error (in rows) applied to simulated catch and effort.", label = "tab:SummaryOfDistributionOfNaturalMortality")
print(tmp.NatMort.xtable, file = "Results/Tables/SummaryOfDistributionOfNaturalMortality.tex", floating.environment = 'sidewaystable')

# discrepancy in percentage of true value
tmp.NatMort2 <- with(tmp, aggregate(abs(Est.NatMort - Sim.NatMort)/Sim.NatMort, by = list(WhiteNoiseLevel), function(x) quantile(x, c(seq(0.1,0.9,0.1), 0.95, 0.975, 0.99, 1), na.rm = TRUE)))

# Plot histogram
library(ggplot2)
p0 <- ggplot(tmp, aes(x=(Est.NatMort - Sim.NatMort)/sd.Est.NatMort)) + geom_histogram(aes(y=..density..), binwidth = 1)
p0 <- p0 + stat_function(fun = dnorm, aes(colour = 'Normal'))
p0 <- p0 + facet_wrap(~WhiteNoiseLevel)
print(p0)
postscript(file = "Results/Graphics/DistributionOfStandardizedNaturalMortEstimates.ps")
print(p0)
dev.off()


#### Targeted catchability

# Discrepancy in unit of standard deviation
tmp.TargCatch <- with(tmp, aggregate(abs(Est.targeted.q - Sim.targeted.q)/sd.Est.targeted.q, by = list(WhiteNoiseLevel), function(x) quantile(x, c(seq(0.1,0.9,0.1), 0.95, 0.975, 0.99, 1), na.rm = TRUE)))

# Discrepancy in percentage of the true value
tmp.TargCatch2 <- with(tmp, aggregate(abs(Est.targeted.q - Sim.targeted.q)/Sim.targeted.q, by = list(WhiteNoiseLevel), function(x) quantile(x, c(seq(0.1,0.9,0.1), 0.95, 0.975, 0.99, 1), na.rm = TRUE)))


##### Sigma of the fit (sd deviation of the normal errors around catch)
tmp.sigma <- with(tmp, aggregate(Sigma, by = list(WhiteNoiseLevel), function(x) quantile(x, c(seq(0.1,0.9,0.1), 0.95, 0.975, 0.99, 1), na.rm = TRUE)))

tmp.sigma.xtable <- as.matrix(tmp.sigma[,-1])
dimnames(tmp.sigma.xtable)[[1]] <- paste("White noise level ", substr(tmp.sigma[,1], 16, 20), "%")
tmp.sigma.xtable <- xtable(tmp.sigma.xtable, caption = "Quantiles of the distribution of the residuals of the fit ($\\sigma$) for varying levels of random error (in rows) applied to simulated catch and effort.", label = "tab:SummaryOfDistributionOfSigma")
print(tmp.sigma.xtable, file = "Results/Tables/SummaryOfDistributionOfSigma.tex", floating.environment = 'sidewaystable')
