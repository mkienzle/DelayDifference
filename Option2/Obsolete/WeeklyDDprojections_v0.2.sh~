#!/bin/bash

# USAGE in ~/mystuff/Work/DEEDI/Moreton Bay Prawn Trawl fishery/Analysis/Scripts/DelayDifferenceModel/1989-2010/Weekly/model7
#          WhatIsMSY.sh Results/ParameterEstimates.txt Data/WeeklyPercentageOfSpawners.txt Results/LinearizedRickerCoef Data/AverageEffortPattern.txt Data/Availability.txt

#  CREATED  19 November 2013
#  MODIFIED 26 November 2013

# Get only the parameter estimates from the result file
cut -f2 -d"," $1 > /tmp/ParEstimates.txt

for i in `seq 0 100 15000`; 
do echo $i;
WeeklyDDprojections /tmp/ParEstimates.txt $2 $3 $4 $5 $i;
mkdir -p Results/Simulation/$i;
cp Results/Simulation/Projections.txt Results/Simulation/$i/.;
done
