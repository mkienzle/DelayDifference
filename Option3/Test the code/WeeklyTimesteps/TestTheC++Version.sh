#!/bin/bash

# CREATED   14 Nov 2013
# MODIFIED   9 Jun 2015

# PURPOSE the present test is to assess whether the weekly delay difference model using MINUIT and coded in C++ 
#         is capable of estimating several parameters INCLUDING NATURAL MORTALITY

for j in $(seq 10 10 40);
do
echo "##### I am working on white noise level $j #####"

for i in $(seq 1 1000); 
do 

# Inform users

echo "##### I am working on iteration $i #####"

# Simulate dataset using random natural mortality
Rscript SimulatePopDynamic2.R $j

# Fit the delay difference model to simulated data
echo "##### I am working on iteration $i #####"
DelayDifference_Option3 Data/SimData4.txt FixedWeeklyParameters.txt
cat Results/SimPar.txt

# Store the results
if [[ ! -e Results/WhiteNoiseLevel$j/$i ]]; then
            mkdir -p Results/WhiteNoiseLevel$j/$i
fi

cp Data/SimData4.txt Results/WhiteNoiseLevel$j/$i/.
cp Data/SimulatedBiomass.csv Results/WhiteNoiseLevel$j/$i/.

cp Results/FitOutcome.txt Results/WhiteNoiseLevel$j/$i/.
cp Results/ParameterEstimates.txt Results/WhiteNoiseLevel$j/$i/.
cp Results/EstimatedFisheriesQuantities.txt Results/WhiteNoiseLevel$j/$i/.
cp Results/SimPar.txt Results/WhiteNoiseLevel$j/$i/.

done
done # end of j loop over varying levels of white noise

# Automatically check for MINUIT failure to fit
for j in $(seq 10 10 40);
do
./FindFitFailure.sh WhiteNoiseLevel$j # cat FailedFit.txt
done

# Plot the results of simulation
for j in $(seq 10 10 40);
do
Rscript PlotTheC++Results.R WhiteNoiseLevel$j
done

