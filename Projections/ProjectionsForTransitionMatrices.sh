#!/bin/bash

# USAGE ProjectionsForTransitionMatrices.sh

# PURPOSE run projections over 50 years multiple times (i times) for effort varying within specific ranges

#  CREATED  19 November 2013
#  MODIFIED 18 October 2017

### Define effort range
LowerEffortValues=(0 1 $(seq 1001 1000 29001))
UpperEffortValue=($(seq 0 1000 30000))

# and mid-range value
EffortValues=(0 $(seq 500 1000 29500))

### Perform i number of simulations over all range of effort
for j in `seq 0 $((${#EffortValues[@]} - 1))`;
	 do
	     for i in `seq 1 10000`; 
	     do
		echo -ne $j"/"$i'\r';
		50YearsProjectionsWithEffortBetween ${LowerEffortValues[$j]} ${UpperEffortValue[$j]} ;
		mkdir -p Results/Simulation/ForTransitionMatrices/${EffortValues[$j]};
		cp Results/Simulation/Projections_forMDP_byYear.csv Results/Simulation/ForTransitionMatrices/${EffortValues[$j]}/Projections_byYear$i.txt;
	     done
done

# Extract catch, effort and recruitment from year before last and recruitment in the last year
echo "Now extracting results from each simulation file"
for j in `seq 0 $((${#EffortValues[@]} - 1))`;
	 do

	 # A file to save all the last transition between recruitment
	 file="Results/Simulation/ForTransitionMatrices/${EffortValues[$j]}/FinalTransition.txt"
	 if [ -f $file ]; then
	     rm $file
	 fi

	 # create file with header
	 echo "Catch,Effort,St,St+1" > $file
	 echo -ne $j'\r';
	 
	 ls Results/Simulation/ForTransitionMatrices/${EffortValues[$j]}/Projections_byYear* | while read filename;
	 do
	     # Extract data from year before last and last year of simulation 
	     echo `tail -2 $filename | head -1 | cut -d, -f3,4,5`","`tail -1 $filename | cut -d, -f5` >> $file
	 done
done

# Combine all together

finalFile="Results/Simulation/ForTransitionMatrices/DataToBuildTransitionMatrices.csv"

	 if [ -f $finalFile ]; then
	     rm $finalFile
	 fi

	 echo "Catch,Effort,St,St+1" > $finalFile

	 find Results/Simulation/ForTransitionMatrices/ -name FinalTransition.txt | while read datafile; do

	     tail -n +2 $datafile >> $finalFile;
	 done
	 
