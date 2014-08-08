# CREATED  23 April 2014
# MODIFIED 29 April 2014

# AUTHOR Marco.Kienzle@gmail.com

# PURPOSE build the weekly delay difference model

CFLAGS = -std=c++0x -O2 -I/usr/include $(shell root-config --glibs --cflags) -lMinuit2

##############################################################################################################################################################
### 1. A program to fit the delay difference model to data and estimate parameters
##############################################################################################################################################################

FitWeeklyDelayDifference_v0.2: vonMisesRecDist.o FitWeeklyDelayDifference.o LogLikelihoodFunction.o WeeklyDelayDifference.o UsefulFunctions.o
			g++ -o FitWeeklyDelayDifference_v0.2 /usr/lib/prob.o vonMisesRecDist.o FitWeeklyDelayDifference.o LogLikelihoodFunction.o WeeklyDelayDifference.o UsefulFunctions.o $(CFLAGS)

UsefulFunctions.o: std_lib_facilities.h UsefulFunctions.h
			g++ -c UsefulFunctions.cpp

vonMisesRecDist.o: vonMisesRecDist.cxx
	g++ $(CFLAGS) -c vonMisesRecDist.cxx

LogLikelihoodFunction.o: LogLikelihoodFunction.cxx LogLikelihoodFunction.h WeeklyDelayDifference.h 
	g++ $(CFLAGS) -c LogLikelihoodFunction.cxx

WeeklyDelayDifference.o: WeeklyDelayDifference.cxx
	g++ $(CFLAGS) -c WeeklyDelayDifference.cxx

FitWeeklyDelayDifference.o: FitWeeklyDelayDifference.cxx WeeklyDelayDifference.h LogLikelihoodFunction.h
	g++ $(CFLAGS) -c FitWeeklyDelayDifference.cxx

##############################################################################################################################################################
### 2. A program to simulate a fishery with the delay difference model based on parameters estimated for this model 
##############################################################################################################################################################

WeeklyDDprojections_v02: WeeklyDDprojections.cpp /usr/lib/prob.o Projections.cxx vonMisesRecDist.o WeeklyDDprojections.o UsefulFunctions.o
			g++ -std=c++0x -g -o WeeklyDDprojections_v0.2 -I/usr/include /usr/lib/prob.o Projections.cxx vonMisesRecDist.o WeeklyDDprojections.o UsefulFunctions.o

##############################################################################################################################################################
### make all
##############################################################################################################################################################

all: FitWeeklyDelayDifference_v0.2 WeeklyDDprojections_v02

##############################################################################################################################################################
### make clean
##############################################################################################################################################################

clean:
	rm -f *.o FitWeeklyDelayDifference_v0.2 WeeklyDDprojections_v02

##############################################################################################################################################################
### make test
##############################################################################################################################################################

test:
	./FitWeeklyDelayDifference_v0.2 Data/SimData4.txt > /tmp/output
	tail --lines=+3 /tmp/output > /tmp/file1; tail --lines=+3 Results/v0.1.output > /tmp/file2
	diff /tmp/file1 /tmp/file2 # header differ by the version number
