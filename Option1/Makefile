# CREATED  23 April 2014
# MODIFIED 12 August 2014

# AUTHOR Marco.Kienzle@gmail.com

# PURPOSE build the weekly delay difference model

#CFLAGS = -std=c++0x -O2 -I/usr/include $(shell root-config --glibs --cflags) -lMinuit2

CFLAGS = -std=c++0x -O2 -I/usr/include -L/usr/local/root/lib/root -I/usr/local/root/include/root -lMinuit2

##############################################################################################################################################################
### 1. A program to fit the delay difference model to data and estimate parameters
##############################################################################################################################################################

DelayDifference_Option1: FitDelayDifference.o LogLikelihoodFunction.o DelayDifference.o UsefulFunctions.o
			g++ -o DelayDifference_Option1 /usr/lib/prob.o FitDelayDifference.o LogLikelihoodFunction.o DelayDifference.o UsefulFunctions.o $(CFLAGS)

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

DD_Option1Projections: /usr/lib/prob.o DDprojections.o Projections.o UsefulFunctions.o
			g++ -std=c++0x -g -o DD_Option1Projections -I/usr/include /usr/lib/prob.o DDprojections.o Projections.o UsefulFunctions.o

DDprojections.o: DDprojections.cpp
			g++ $(CFLAGS) -c DDprojections.cpp

Projections.o: Projections.cxx
			g++ $(CFLAGS) -c Projections.cxx

##############################################################################################################################################################
### make all
##############################################################################################################################################################

all: DelayDifference_Option1 DD_Option1Projections

##############################################################################################################################################################
### make clean
##############################################################################################################################################################

clean:
	rm -f *.o DelayDifference_Option1 DDprojections_Option1Projections

##############################################################################################################################################################
### make test
##############################################################################################################################################################

test:
	./DelayDifference_Option1 Data/SimData4.txt > /tmp/output
	tail --lines=+3 /tmp/output > /tmp/file1; tail --lines=+3 Results/v0.1.output > /tmp/file2
	diff /tmp/file1 /tmp/file2 # header differ by the version number
