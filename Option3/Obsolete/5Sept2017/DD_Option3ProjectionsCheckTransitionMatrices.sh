#!/bin/bash

# USAGE DD_Option3ProjectionsCheckTransitionMatrices.sh Results/DelayDifferenceModelParameters.csv Data/WeeklyPercentageOfSpawners.txt Results/LinearizedRickerCoef Data/AverageEffortPattern.txt Data/Availability.txt

# PURPOSE run projections over x years multiple times for effort varying within a range

#  CREATED  19 November 2013
#  MODIFIED 22 August 2017

### Define effort range
LowerEffortValues=(0 1 1000 2000 3000 4000 5000 6000 7000 8000 9000 10000 11000 12000 13000 14000 15000 16000 17000 18000 19000 20000 21000 22000 23000 24000 25000 26000 27000 28000 29000)
UpperEffortValue=(0 1000 2000 3000 4000 5000 6000 7000 8000 9000 10000 11000 12000 13000 14000 15000 16000 17000 18000 19000 20000 21000 22000 23000 24000 25000 26000 27000 28000 29000 30000)

# and mid-range value
EffortValues=( 0 500 1500 2500 3500 4500 5500 6500 7500 8500 9500 10500 11500 12500 13500 14500 15500 16500 17500 18500 19500 20500 21500 22500 23500 24500 25500 26500 27500 28500 29500)

### Perform i number of simulations over all range of effort
for j in `seq 0 30`;
	 do
	     for i in `seq 1 2000`; 
	     do
		echo -ne $j"/"$i'\r';
		DD_Option3Projections_forMDPWithEffortLimits $1 $2 $3 $4 $5 50 ${LowerEffortValues[$j]} ${UpperEffortValue[$j]} ;
		mkdir -p Results/Simulation/Checks/${EffortValues[$j]};
		cp Results/Simulation/Projections_forMDP_byYearWithEffortLimits.txt Results/Simulation/Checks/${EffortValues[$j]}/Projections_byYear$i.txt;
	     done
done

# Extract recruitment in the last year
echo "Now extracting results from each simulation file"
for j in `seq 0 30`;
	 do
	 file="Results/Simulation/Checks/${EffortValues[$j]}/FinalRec.txt"
	 if [ -f $file ]; then
	     rm $file
	 fi

	 touch $file
	 echo -ne $j'\r';
	 
	 ls Results/Simulation/Checks/${EffortValues[$j]}/Projections_byYear* | while read filename;
	 do
	     tail -1 $filename | cut -d',' -f5 >> $file
	 done
done
