#!/bin/bash

# USAGE ProjectionsForMAY.sh 


#  CREATED  19 November 2013
#  MODIFIED  1 September 2017

# Get only the parameter estimates from the result file
#cut -f2 -d"," $2 > /tmp/ParEstimates.txt

tmpfile="/tmp/TmpProj.txt"

if [ -f $tmpfile ]; then
   rm $tmpfile;
fi

   echo "Catch,Effort,SSB" > $tmpfile

for i in 0 `seq 500 1000 29500`; 
do
    #echo $i;

   for j in `seq 1 5000`; do
   echo -ne $i"/"$j'\r'
   50YearsProjectionsForMAY $i;
   tail -1 Results/Simulation/Projections_forMDP_byYear.csv | cut -f3,4,6 -d, >> $tmpfile

   done
   
#mkdir -p Results/Simulation/$i;
#cp Results/Simulation/Projections.txt Results/Simulation/$i/.;
done

mv $tmpfile Results/Simulation/CatchAndEffortForMAY.csv