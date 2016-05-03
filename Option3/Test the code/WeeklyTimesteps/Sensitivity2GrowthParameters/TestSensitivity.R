# CREATED  10 October 2014
# MODIFIED 10 October 2014

# PURPOSE test delay difference model response to variation in somatic growth parameters input (rho, wk, wk_1)
#                                                           and age at recruitment
# Range of variation: for von Bertalanffy parameters +- 10%, for age at recruitment between 19 and 25 weeks

# Load useful functions
source("UsefulFunctions.R")

# An object holding the results
v.Linf <- seq(0.9,1.2,0.05)
v.k <- seq(0.9,1.1,0.05)

results <- expand.grid( vb.Linf.factor = v.Linf, vb.k.factor = v.k, rec.age = seq(20,24))
results[,"ExitStatus"] <- NA
results[,"logL"] <- NA

# Generate a simulated dataset
source("SimulatePopDynamic2.R")

# Loop over, create the input file and fit the model
for(i in 1:nrow(results)){
    
    fixed.par <- GenerateSomaticGrowthParameters(results[i,"vb.Linf.factor"],  results[i,"vb.k.factor"], results[i,"rec.age"])
    write.inputfile(fixed.par)
    #fit <- system("DelayDifference_Option3 Data/SimData4.txt FixedWeeklyParameters.txt", intern=FALSE)
    fit <- system("DelayDifference_Option3 Data/SimData4.txt FixedWeeklyParameters.txt", intern=FALSE)
    results[i,"ExitStatus"] <- fit

    # Get outcome of the fit if successful
    if(!fit) results[i,"logL"] <- read.table(file = "Results/FitOutcome.txt", skip = 3, nrow =1, sep = ":")$V2

}

with(results, which(logL == min(logL, na.rm = T)))
with(subset(results, rec.age ==22), which(logL == min(logL, na.rm = T)))

a.df <- with(subset(results, rec.age ==22), tapply(logL, list(vb.Linf.factor, vb.k.factor), I))
par(mfrow=c(1,1))
contour(v.Linf, v.k, a.df - min(a.df, na.rm=TRUE), levels = 10^seq(-8,0), xlab = "VB Linf multiplier", ylab ="VB k multiplier")
print(a.df - min(a.df, na.rm = T))
