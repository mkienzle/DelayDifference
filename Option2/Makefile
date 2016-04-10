# CREATED  23 April 2014
# MODIFIED  5 April 2016

# AUTHOR Marco.Kienzle@gmail.com

# PURPOSE build the weekly delay difference model

#CFLAGS = -std=c++0x -O2 -I/usr/include $(shell root-config --glibs --cflags) -lMinuit2
#CFLAGS = -std=c++0x -O2 -I/usr/include -L/usr/local/root/lib/root -I/usr/local/root/include/root -lMinuit2

# Compile using static libraries for 32bit machine (NOT WORKING because of libMinuit2 not being 32bit)
# CFLAGS = -m32 -g -O2 -Wall -fopenmp -std=c++0x -O2 -I/usr/include -L/usr/local/lib -static -lMinuit2 # with lots of warning
# CFLAGS = -m32 -fopenmp -std=c++0x -O2 -I/usr/include -L/usr/local/lib -static -lMinuit2

# Compile using static libraries for 64bit machine
# CFLAGS = -m32 -g -O2 -Wall -fopenmp -std=c++0x -O2 -I/usr/include -L/usr/local/lib -static -lMinuit2 # with lots of warning
CFLAGS = -fopenmp -std=c++0x -O2 -I/usr/include -L/usr/local/lib -static -lMinuit2

############################################################################################################################################
### 1. A program to fit the delay difference model to data and estimate parameters
############################################################################################################################################

DelayDifference_Option2: FitDelayDifference.o LogLikelihoodFunction.o DelayDifference.o UsefulFunctions.o
			g++ -o DelayDifference_Option2 /usr/lib/prob.o FitDelayDifference.o LogLikelihoodFunction.o DelayDifference.o UsefulFunctions.o $(CFLAGS)

UsefulFunctions.o: ../lib_facilities.h ../UsefulFunctions.h ../UsefulFunctions.cpp
			g++ $(CFLAGS) -c ../UsefulFunctions.cpp

#vonMisesRecDist.o: ../vonMisesRecDist.cpp
#	g++ $(CFLAGS) -c ../vonMisesRecDist.cpp

LogLikelihoodFunction.o: LogLikelihoodFunction.cpp LogLikelihoodFunction.h DelayDifference.h 
	g++ $(CFLAGS) -c LogLikelihoodFunction.cpp

DelayDifference.o: DelayDifference.cpp
	g++ $(CFLAGS) -c DelayDifference.cpp

FitDelayDifference.o: FitDelayDifference.cpp DelayDifference.h LogLikelihoodFunction.h
	g++ $(CFLAGS) -c FitDelayDifference.cpp

############################################################################################################################################
### 2. A program to simulate a fishery with the delay difference model based on parameters estimated for this model 
############################################################################################################################################

DD_Option2Projections: /usr/lib/prob.o DDprojections.o Projections.o UsefulFunctions.o
			g++ -std=c++0x -g -o DD_Option2Projections -I/usr/include /usr/lib/prob.o DDprojections.o Projections.o UsefulFunctions.o

DDprojections.o: DDprojections.cpp
			g++ $(CFLAGS) -c DDprojections.cpp

Projections.o: Projections.cxx
			g++ $(CFLAGS) -c Projections.cxx

############################################################################################################################################
### make all
############################################################################################################################################

all: DelayDifference_Option2 DD_Option2Projections

############################################################################################################################################
### make clean
############################################################################################################################################

clean:
	rm -f *.o DelayDifference_Option2 DDprojections_Option2Projections

############################################################################################################################################
### make test
############################################################################################################################################

testDelayDifference:
	./DelayDifference_Option2 Data/SimData4.txt InputParameterDescription.csv
	#diff Results/FitOutcome.txt Results/FitOutcome-1Apr2016.txt
	diff Results/ParameterEstimates.txt Results/ParameterEstimates-2016Apr05.txt
