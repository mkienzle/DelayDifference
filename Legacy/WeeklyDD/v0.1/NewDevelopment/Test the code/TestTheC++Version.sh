#!/bin/bash

# CREATED   14 Nov 2013
# MODIFIED  14 Nov 2013

# PURPOSE assess the capacity of the weekly delay difference model using MINUIT and coded in C++ to estimate parameters used to generate simulated datasets

for i in $(seq 1 100); 
do 

# Inform users
echo "##### I am working on iteration $i #####"

# Simulate data
R --vanilla < SimulatePopDynamic.R

# Fit the delay difference model to simulated data
echo "##### I am working on iteration $i #####"
FitWeeklyDelayDifference Data/SimData4.txt

# Store the results
mkdir -p Results/C++Version/$i

cp Data/SimData4.txt Results/C++Version/$i/.
cp Data/SimulatedBiomass.csv Results/C++Version/$i/.

cp Results/FitOutcome.txt Results/C++Version/$i/.
cp Results/ParameterEstimates.txt Results/C++Version/$i/.
cp Results/EstimatedFisheriesQuantities.txt Results/C++Version/$i/.
cp Results/SimPar.txt Results/C++Version/$i/.

done



