# CREATED  23 April 2014
# MODIFIED 19 July  2014

# AUTHOR Marco.Kienzle@gmail.com

# PURPOSE build the weekly delay difference model

CFLAGS = -std=c++0x -O2 -I/usr/include $(shell root-config --glibs --cflags) -lMinuit2

##############################################################################################################################################################
### 1. A program to fit the delay difference model to data and estimate parameters
##############################################################################################################################################################

DelayDifference_v0.2: FitDelayDifference.o LogLikelihoodFunction.o DelayDifference.o UsefulFunctions.o
			g++ -o DelayDifference_v0.2 /usr/lib/prob.o FitDelayDifference.o LogLikelihoodFunction.o DelayDifference.o UsefulFunctions.o $(CFLAGS)

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

##############################################################################################################################################################
### 2. A program to simulate a fishery with the delay difference model based on parameters estimated for this model 
##############################################################################################################################################################

DD_v0.2Projections: /usr/lib/prob.o DDprojections.o Projections.o UsefulFunctions.o
			g++ -std=c++0x -g -o DD_v0.2Projections -I/usr/include /usr/lib/prob.o DDprojections.o Projections.o UsefulFunctions.o

DDprojections.o: DDprojections.cpp
			g++ $(CFLAGS) -c DDprojections.cpp

Projections.o: Projections.cxx
			g++ $(CFLAGS) -c Projections.cxx

##############################################################################################################################################################
### make all
##############################################################################################################################################################

all: DelayDifference_v0.2 DD_v0.2Projections

##############################################################################################################################################################
### make clean
##############################################################################################################################################################

clean:
	rm -f *.o DelayDifference_v0.2 DDprojections_v02

##############################################################################################################################################################
### make test
##############################################################################################################################################################

test:
	./DelayDifference_v0.2 Data/SimData4.txt > /tmp/output
	tail --lines=+3 /tmp/output > /tmp/file1; tail --lines=+3 Results/v0.1.output > /tmp/file2
	diff /tmp/file1 /tmp/file2 # header differ by the version number
