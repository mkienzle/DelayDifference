#!/bin/bash

# CREATED   14 Nov 2013
# MODIFIED   9 Dec 2013

# PURPOSE assess the capacity of the weekly delay difference model using MINUIT and coded in C++ to estimate parameters used to generate simulated datasets

for i in $(seq 1 100); 
do 

# Inform users
echo "##### I am working on iteration $i #####"

# Simulate data
R --vanilla < SimulatePopDynamic2.R

# Fit the delay difference model to simulated data
echo "##### I am working on iteration $i #####"
FitWeeklyDelayDifference2 Data/SimData4.txt

# Store the results
mkdir -p Results/C++Version2/$i

cp Data/SimData4.txt Results/C++Version2/$i/.
cp Data/SimulatedBiomass.csv Results/C++Version2/$i/.

cp Results/FitOutcome.txt Results/C++Version2/$i/.
cp Results/ParameterEstimates.txt Results/C++Version2/$i/.
cp Results/EstimatedFisheriesQuantities.txt Results/C++Version2/$i/.
cp Results/SimPar.txt Results/C++Version2/$i/.

done



