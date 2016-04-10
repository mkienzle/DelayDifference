#!/bin/bash

# CREATED   14 Nov 2013
# MODIFIED  29 May 2014

# PURPOSE the present test is to assess whether the weekly delay difference model using MINUIT and coded in C++ 
#         is capable of estimating several parameters INCLUDING NATURAL MORTALITY

for i in $(seq 1 100); 
do 

# Inform users
echo "##### I am working on iteration $i #####"

# Simulate dataset using random natural mortality
#R --slave --vanilla < SimulatePopDynamic2.R
R --slave --vanilla < SimulatePopDynamic.R

# Fit the delay difference model to simulated data
echo "##### I am working on iteration $i #####"
FitDelayDifference_v0.3 Data/SimData4.txt FixedMonthlyParameters.txt
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



