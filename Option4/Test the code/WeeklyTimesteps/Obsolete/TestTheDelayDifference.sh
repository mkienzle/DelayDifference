#!/bin/bash 

# CREATED   4 June  2012
# MODIFIED 23 April 2013

# PURPOSE to test the fit of the Delay Difference model in this directory to many set of simulated data

# NOTE the minimization algorithm returned different minimum when ran a second time starting with the first optimum !
#      ------> so a while loop was implemented to make  sure that it end up in the same place

# Plot the results in R
# OBSOLETE tmp <- read.table("Results/Results.txt"); with(tmp, plot( V3/V2, V4, pch = 19, cex = 0.5, las = 1, xlab = "Underestimated <--- Ratio of estimated over simulated catchability ----> Overestimated", ylab = "Underestimated <---- Slope of estimated biomass = f(simulated) ---> overestimated")); points(1,1, pch = 3, col = "red"); abline(v=1, lty=3, col="red"); abline(h=1, lty = 3, col = "red")

# Create file to hold results
echo "SimNumber Catchability" > Results/SimulatedParameters.txt
echo "SimNumber CatchabilityEstimate CatchabilityEstimateUncertainty SlopeOfRegressionBetweenBiomasses" > Results/ResultsByFittingCPUE.txt
echo "SimNumber CatchabilityEstimate CatchabilityEstimateUncertainty SlopeOfRegressionBetweenBiomasses" > Results/ResultsByFittingCatch.txt
echo "SimNumber CatchabilityEstimate CatchabilityEstimateUncertainty SlopeOfRegressionBetweenBiomasses" > Results/ResultsByFittingSqrtOfCatch.txt
echo "SimNumber CatchabilityEstimate CatchabilityEstimateUncertainty SlopeOfRegressionBetweenBiomasses" > Results/ResultsByFittingLogOfCatch.txt

for i in $(seq 1 100); 
do 

echo -n "$i " >> Results/SimulatedParameters.txt;
echo -n "$i " >> Results/ResultsByFittingCPUE.txt; 
echo -n "$i " >> Results/ResultsByFittingCatch.txt; 
echo -n "$i " >> Results/ResultsByFittingSqrtOfCatch.txt; 
echo -n "$i " >> Results/ResultsByFittingLogOfCatch.txt; 

# Simulate a set of data
R --vanilla < SimulatePopDynamic.R

### Fit the model 

# # to CPUE data only
# matlab -nodesktop -nosplash -r "EKPdatasets; EKPparameters; start_value = [ones(1,11)]; [mle,fval1,exitflag1] = fminsearch(@(par) EKPLogLikDelayDifference(par, cpue,1), start_value,  optimset('MaxFunEvals',1e5, 'MaxIter', 1e5, 'Display','iter')); while (sum(abs(mle-start_value)) > 1e-2), [mle,fval1,exitflag1] = fminsearch(@(par) EKPLogLikDelayDifference(par, cpue,1), start_value,  optimset('MaxFunEvals',1e5, 'MaxIter', 1e5, 'Display','iter')); start_value = mle; end; ConfidenceInterval2; csvwrite('Results/EstimatedBiomassByFittingCPUE.csv', Biomass); fid = fopen('Results/ResultsByFittingCPUE.txt', 'at'); fprintf(fid, '%2.2f %2.2f ', mle(1),error(1)); fclose(fid); exit"

# to square root of catch 
matlab -nodesktop -nosplash -r "EKPdatasets; EKPparameters; start_value = [ones(1,11)]; [mle,fval1,exitflag1] = fminsearch(@(par) EKPLogLikOfSqrtOfCatch(par, cpue,1), start_value,  optimset('MaxFunEvals',1e5, 'MaxIter', 1e5, 'Display','iter')); while (sum(abs(mle-start_value)) > 1e-2), [mle,fval1,exitflag1] = fminsearch(@(par) EKPLogLikOfSqrtOfCatch(par, cpue,1), start_value,  optimset('MaxFunEvals',1e5, 'MaxIter', 1e5, 'Display','iter')); start_value = mle; end; ConfidenceInterval2; csvwrite('Results/EstimatedBiomassByFittingSqrtOfCatch.csv', Biomass); fid = fopen('Results/ResultsByFittingSqrtOfCatch.txt', 'at'); fprintf(fid, '%2.2f %2.2f ', mle(1),error(1)); fclose(fid); exit"

# # to catch 
# matlab -nodesktop -nosplash -r "EKPdatasets; EKPparameters; start_value = [ones(1,11)]; [mle,fval1,exitflag1] = fminsearch(@(par) EKPLogLikOfCatch(par, cpue,1), start_value,  optimset('MaxFunEvals',1e5, 'MaxIter', 1e5, 'Display','iter')); while (sum(abs(mle-start_value)) > 1e-2), [mle,fval1,exitflag1] = fminsearch(@(par) EKPLogLikOfCatch(par, cpue,1), start_value,  optimset('MaxFunEvals',1e5, 'MaxIter', 1e5, 'Display','iter')); start_value = mle; end; ConfidenceInterval2; csvwrite('Results/EstimatedBiomassByFittingCatch.csv', Biomass); fid = fopen('Results/ResultsByFittingCatch.txt', 'at'); fprintf(fid, '%2.2f %2.2f ', mle(1),error(1)); fclose(fid); exit"

# # to log of catch 
# matlab -nodesktop -nosplash -r "EKPdatasets; EKPparameters; start_value = [ones(1,11)]; [mle,fval1,exitflag1] = fminsearch(@(par) EKPLogLikOfLogOfCatch(par, cpue,1), start_value,  optimset('MaxFunEvals',1e5, 'MaxIter', 1e5, 'Display','iter')); while (sum(abs(mle-start_value)) > 1e-2), [mle,fval1,exitflag1] = fminsearch(@(par) EKPLogLikOfLogOfCatch(par, cpue,1), start_value,  optimset('MaxFunEvals',1e5, 'MaxIter', 1e5, 'Display','iter')); start_value = mle; end; ConfidenceInterval2; csvwrite('Results/EstimatedBiomassByFittingLogOfCatch.csv', Biomass); fid = fopen('Results/ResultsByFittingLogOfCatch.txt', 'at'); fprintf(fid, '%2.2f %2.2f ', mle(1),error(1)); fclose(fid); exit"


# Save the Results
R --vanilla < CompareEstimate2Simulated.R;
done

R --vanilla < PlotTheResults.R