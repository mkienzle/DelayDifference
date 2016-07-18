#!/bin/bash

# USAGE DD_Option3Projections.sh Results/DelayDifferenceModelParameters.csv Data/WeeklyPercentageOfSpawners.txt Results/LinearizedRickerCoef Data/AverageEffortPattern.txt Data/Availability.txt


#  CREATED  19 November 2013
#  MODIFIED 17 December 2013

# Get only the parameter estimates from the result file
#cut -f2 -d"," $2 > /tmp/ParEstimates.txt

for i in `seq 0 100 20000`; 
do echo $i;
DD_Option3Projections $1 $2 $3 $4 $5 $i;
mkdir -p Results/Simulation/$i;
cp Results/Simulation/Projections.txt Results/Simulation/$i/.;
done
