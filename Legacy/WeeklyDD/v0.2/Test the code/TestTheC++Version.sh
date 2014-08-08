#!/bin/bash

# CREATED   14 Nov 2013
# MODIFIED  29 May 2014

# PURPOSE assess the capacity of the weekly delay difference model using MINUIT and coded in C++ to estimate parameters used to generate simulated datasets

for i in $(seq 1 100); 
do 

# Inform users
echo "##### I am working on iteration $i #####"

# Simulate data
R --vanilla < SimulatePopDynamic.R

# Fit the delay difference model to simulated data
echo "##### I am working on iteration $i #####"
FitWeeklyDelayDifference_v0.2 Data/SimData4.txt
cat Results/SimPar.txt

# Store the results
mkdir -p Results/$i

cp Data/SimData4.txt Results/$i/.
cp Data/SimulatedBiomass.csv Results/$i/.

cp Results/FitOutcome.txt Results/$i/.
cp Results/ParameterEstimates.txt Results/$i/.
cp Results/EstimatedFisheriesQuantities.txt Results/$i/.
cp Results/SimPar.txt Results/$i/.

done


