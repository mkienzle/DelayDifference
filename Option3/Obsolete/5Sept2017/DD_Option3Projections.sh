#!/bin/bash

# USAGE DD_Option3Projections.sh Results/DelayDifferenceModelParameters.csv Data/WeeklyPercentageOfSpawners.txt Results/LinearizedRickerCoef Data/AverageEffortPattern.txt Data/Availability.txt


#  CREATED  19 November 2013
#  MODIFIED 17 December 2013

# Get only the parameter estimates from the result file
#cut -f2 -d"," $2 > /tmp/ParEstimates.txt

tmpfile="/tmp/TmpProj.txt"

if [ -f $tmpfile ]; then
   rm $tmpfile;
fi

   echo "catch,effort,SSB" > $tmpfile

for i in `seq 0 1000 30000`; 
do
    #echo $i;

   for j in `seq 1 10000`; do
   echo -ne $i"/"$j'\r'
   DD_Option3Projections $1 $2 $3 $4 $5 $i;
   tail -1 Results/Simulation/Projections_byYear.txt | cut -f3,4,6 -d, >> $tmpfile

   done
   
#mkdir -p Results/Simulation/$i;
#cp Results/Simulation/Projections.txt Results/Simulation/$i/.;
done

mv $tmpfile Results/Simulation/CatchAndEffortForMAY.csv
