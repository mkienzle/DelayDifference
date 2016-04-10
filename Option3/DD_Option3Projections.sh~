#!/bin/bash

# USAGE DDprojections_v0.3.sh FixedWeeklyParameters.txt Results/ParameterEstimates.txt Data/WeeklyPercentageOfSpawners.txt Results/LinearizedRickerCoef Data/AverageEffortPattern.txt Data/Availability.txt


#  CREATED  19 November 2013
#  MODIFIED 17 December 2013

# Get only the parameter estimates from the result file
cut -f2 -d"," $2 > /tmp/ParEstimates.txt

for i in `seq 0 100 20000`; 
do echo $i;
DD_v0.3Projections $1 /tmp/ParEstimates.txt $3 $4 $5 $6 $i;
mkdir -p Results/Simulation/$i;
cp Results/Simulation/Projections.txt Results/Simulation/$i/.;
done
