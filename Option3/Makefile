# CREATED  23 April 2014
# MODIFIED  6 April 2016

# AUTHOR Marco.Kienzle@gmail.com

# PURPOSE build the weekly delay difference model

CFLAGS = -fopenmp -std=c++0x -O2 -I/usr/include -L/usr/local/lib -static -lMinuit2
#CFLAGS = -std=c++0x -O2 -I/usr/include $(shell root-config --glibs --cflags) -lMinuit2
#CFLAGS = -std=c++0x -O2 -I/usr/include -L/usr/local/lib/root -I/usr/local/include/root -lMinuit2

### Attempt to build a static binary (doesn't work!)
#CFLAGS = -std=c++0x -O2 -I/usr/include -L/usr/local/lib -L/usr/local/lib/root -I/usr/local/include/root -lRoot -lMinuit2 -static
#CFLAGS = -std=c++0x -O2 -I/usr/include -L/usr/local/lib/root -I/usr/local/include/root -lRoot -lMinuit2

##############################################################################################################################################################
### 1. A program to fit the delay difference model to data and estimate parameters
##############################################################################################################################################################

DelayDifference_Option3: FitDelayDifference_Option3.o LogLikelihoodFunction.o DelayDifference.o UsefulFunctions.o
			g++ -o DelayDifference_Option3 /usr/lib/prob.o FitDelayDifference_Option3.o LogLikelihoodFunction.o DelayDifference.o UsefulFunctions.o $(CFLAGS)

UsefulFunctions.o: ../lib_facilities.h ../UsefulFunctions.h ../UsefulFunctions.cpp
			g++ $(CFLAGS) -c ../UsefulFunctions.cpp

#vonMisesRecDist.o: ../vonMisesRecDist.cpp
#	g++ $(CFLAGS) -c ../vonMisesRecDist.cpp

LogLikelihoodFunction.o: LogLikelihoodFunction.cpp LogLikelihoodFunction.h DelayDifference.h 
	g++ $(CFLAGS) -c LogLikelihoodFunction.cpp

DelayDifference.o: DelayDifference.cpp
	g++ $(CFLAGS) -c DelayDifference.cpp

FitDelayDifference_Option3.o: FitDelayDifference_Option3.cpp DelayDifference.h LogLikelihoodFunction.h
	g++ $(CFLAGS) -c FitDelayDifference_Option3.cpp

##############################################################################################################################################################
### 2. A program to simulate a fishery with the delay difference model based on parameters estimated for this model 
##############################################################################################################################################################

DD_Option3Projections: /usr/lib/prob.o DDprojections2.o Projections2.o UsefulFunctions.o
			g++ -std=c++0x -g -o DD_Option3Projections -I/usr/include /usr/lib/prob.o DDprojections2.o Projections2.o UsefulFunctions.o

DDprojections2.o: DDprojections2.cpp
			g++ $(CFLAGS) -c DDprojections2.cpp

Projections2.o: Projections2.cpp
			g++ $(CFLAGS) -c Projections2.cpp

##############################################################################################################################################################
### make all
##############################################################################################################################################################

all: DelayDifference_Option3 DD_Option3Projections

##############################################################################################################################################################
### make clean
##############################################################################################################################################################

clean:
	rm -f *.o DelayDifference_Option3 DD_Option3Projections

############################################################################################################################################
### make test
###########################################################################################################################################
testDelayDifference:
	./DelayDifference_Option3 Data/SimData4.txt InputParameterDescription.csv
	#diff Results/FitOutcome.txt Results/FitOutcome-2016Apr6.txt
	diff Results/ParameterEstimates.txt Results/ParameterEstimates-2016Apr6.txt


