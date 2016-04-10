#!/bin/sh

R --vanilla < SimulatePopDynamic.R

# to square root of catch 
matlab -nodesktop -nosplash -r "SimulatedDatasets; MbTigerPrawnParameters; par = ones(1,2); DelayDifference; csvwrite('Results/EstimatedBiomassByFittingSqrtOfCatch.csv', Biomass); exit"

R --vanilla < CompareEstimate2Simulated.R
