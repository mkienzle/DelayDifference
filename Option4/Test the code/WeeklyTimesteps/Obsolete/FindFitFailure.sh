#!/bin/bash

# CREATED   9 October 2014
# MODIFIED 10 June 2015

# PURPOSE automatically look for 

# USAGE ./FindFitFailure WhiteNoiseLevel10

args=("$@")

OUTPUT="Results/${args[0]}/FailedFit.txt"
rm -f $OUTPUT

maxSim=100;

# Look for failure to converge
for i in `seq 1 $maxSim`
 do 
    if grep -q "did not converge" Results/${args[0]}/$i/FitOutcome.txt
    then echo $i >> $OUTPUT
    fi
done

# Look parameter at upper/lower limit
for i in `seq 1 $maxSim`
 do 
    if grep -q "Parameter is at Lower limit" Results/${args[0]}/$i/FitOutcome.txt
    then echo $i >> $OUTPUT
    fi

    if grep -q "Parameter is at Upper limit" Results/${args[0]}/$i/FitOutcome.txt
    then echo $i >> $OUTPUT
    fi

done

# Find von mises mean's estimate close to pi or multiple
for i in `seq 1 $maxSim`
 do 
    if grep -q " -3.14" Results/${args[0]}/$i/FitOutcome.txt
    then echo $i >> $OUTPUT
    fi

    if grep -q " -6.28" Results/${args[0]}/$i/FitOutcome.txt
    then echo $i >> $OUTPUT
    fi

    if grep -q " 3.14" Results/${args[0]}/$i/FitOutcome.txt
    then echo $i >> $OUTPUT
    fi

    if grep -q " 6.28" Results/${args[0]}/$i/FitOutcome.txt
    then echo $i >> $OUTPUT
    fi

done

