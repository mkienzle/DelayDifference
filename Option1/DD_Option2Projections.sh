#!/bin/bash

#  CREATED  19 November 2013
#  MODIFIED 29 April 2014

# USAGE DDprojections_Option2.sh FixedWeeklyParameters.txt Results/ParameterEstimates.txt Data/WeeklyPercentageOfSpawners.txt Results/LinearizedRickerCoef Data/AverageEffortPattern.txt Data/Availability.txt

# OBSOLETE USAGE in ~/mystuff/Work/DEEDI/Moreton Bay Prawn Trawl fishery/Analysis/Scripts/DelayDifferenceModel/1989-2010/Weekly/model7
# OBSOLETE         WhatIsMSY.sh Results/ParameterEstimates.txt Data/WeeklyPercentageOfSpawners.txt Results/LinearizedRickerCoef Data/AverageEffortPattern.txt Data/Availability.txt


# Get only the parameter estimates from the result file
cut -f2 -d"," $2 > /tmp/ParEstimates.txt

for i in `seq 0 100 15000`; 
do echo $i;
DD_Option2Projections $1 /tmp/ParEstimates.txt $3 $4 $5 $6 $i;
mkdir -p Results/Simulation/$i;
cp Results/Simulation/Projections.txt Results/Simulation/$i/.;
done
