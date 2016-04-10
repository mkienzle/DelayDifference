#!/bin/bash

#  CREATED  19 Nov 2013
#  MODIFIED  6 Apr 2016

# USAGE DD_Option2Projections.sh Results/DelayDifferenceModelParameters.csv Data/WeeklyPercentageOfSpawners.txt Results/LinearizedRickerCoef Data/AverageEffortPattern.txt Data/Availability.txt

for i in `seq 0 100 15000`; 
do echo $i;
DD_Option2Projections $1 $2 $3 $4 $5 $i;
mkdir -p Results/Simulation/$i;
cp Results/Simulation/Projections.txt Results/Simulation/$i/.;
done
