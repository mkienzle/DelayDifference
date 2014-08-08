## Compare simulated biomasses with estimated biomasses
time.series1 <- read.csv("Results/EstimatedBiomassByFittingSqrtOfCatch.csv", header = F)[105:520,1];
time.series2 <- read.csv("Data/SimulatedBiomass.csv")[seq(2860-8*52+1, 2860),2]

linear.model <- lm(time.series1 ~ time.series2)

write(round(coef(linear.model)[2],2), file = "Results/ResultsByFittingSqrtOfCatch.txt", append = TRUE)


