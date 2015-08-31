#!/bin/sh

# Store current directory
OriginalDir=`pwd`;

# Create file to hold results
echo "SimNumber Catchability" > Results/SimulatedParameters.txt
echo "SimNumber CatchabilityEstimate CatchabilityEstimateUncertainty SlopeOfRegressionBetweenBiomasses" > Results/Fitting.txt

for i in $(seq 1 1); 
do 

echo -n "$i " >> Results/SimulatedParameters.txt;
echo -n "$i " >> Results/Fitting.txt;

R --vanilla < SimulatePopDynamic.R

# to square root of catch 
matlab -nodesktop -nosplash -r "SimulatedDatasets; MbTigerPrawnParameters; start_value = [ones(1,1)]; [mle,fval,exitflag] = fminsearch(@(par) FitDDusingLogLikOfSqrtOfCatch(par, true), [ones(1,1)],  optimset('MaxFunEvals',1e5, 'MaxIter', 1e5, 'Display','iter')); while (sum(abs(mle-start_value)) > 1e-2),  [mle,fval,exitflag] = fminsearch(@(par) FitDDusingLogLikOfSqrtOfCatch(par, true), [ones(1,1)],  optimset('MaxFunEvals',1e5, 'MaxIter', 1e5, 'Display','iter'));  start_value = mle; end; csvwrite('Results/EstimatedBiomass.csv', Biomass); tmle = transpose([mle fval exitflag]); save('Results/MLE.txt', 'tmle','-ascii'); exit;"

# Estimate error
cd ErrorEstimate;

matlab -nodesktop -nosplash -r "tmpMLE=csvread('../Results/MLE.txt'); start_value = transpose(tmpMLE(1)); SimulatedDatasets; MbTigerPrawnParametersBasedOnPreviousRun; [mle,fval,exitflag] = fminsearch(@(par) FitDDusingLogLikOfSqrtOfCatch(par, true), start_value,  optimset('MaxFunEvals',1e5, 'MaxIter', 1e5, 'Display','iter')); ConfidenceInterval2; csvwrite('../Results/EstimatedBiomassFromErrorAnalysis.csv', Biomass); fid = fopen('../Results/Fitting.txt', 'at'); fprintf(fid, '%2.2f %2.2f ', mle(1),error(1)); fclose(fid); zz = [fval exitflag]; save('../Results/convergence.txt', 'zz', '-ascii'); exit;";

cd "$OriginalDir"
R --vanilla < CompareEstimate2Simulated.R
done