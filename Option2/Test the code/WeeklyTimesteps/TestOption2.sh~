#!/bin/bash

# CREATED   14 Nov 2013
# MODIFIED  26 May 2014

# PURPOSE assess the capacity of the delay difference model v0.2 [using MINUIT and coded in C++] to estimate parameters of a simulated population which data are available on weekly timestepd 

for i in $(seq 1 100);
do 

# Inform users
echo "##### I am working on iteration $i #####"

# Simulate data
R --slave --vanilla < SimulatePopDynamic.R

# Fit the delay difference model to simulated data
echo "##### I am working on iteration $i #####"
DelayDifference_v0.2 Data/SimData4.txt FixedWeeklyParameters.txt
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



