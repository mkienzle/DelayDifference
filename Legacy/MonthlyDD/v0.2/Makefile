# CREATED  23 April 2014
# MODIFIED  9 May 2014

# AUTHOR Marco.Kienzle@gmail.com

# PURPOSE build the weekly delay difference model

CFLAGS = -std=c++0x -O2 -I/usr/include $(shell root-config --glibs --cflags) -lMinuit2

##############################################################################################################################################################
### 1. A program to fit the delay difference model to data and estimate parameters
##############################################################################################################################################################

#g++ -std=c++0x -g -o FitMonthlyDelayDifference -I/usr/include /usr/lib/prob.o vonMisesRecDist.cxx FitMonthlyDelayDifference.cxx LogLikelihoodFunction.cxx MonthlyDelayDifference.cxx `root-config --glibs --cflags` -lMinuit2

FitMonthlyDelayDifference_v0.2: FitMonthlyDelayDifference_v0.2.cxx vonMisesRecDist.o FitMonthlyDelayDifference_v0.2.o LogLikelihoodFunction.o MonthlyDelayDifference.o UsefulFunctions.o
			g++ -o FitMonthlyDelayDifference_v0.2 /usr/lib/prob.o vonMisesRecDist.o FitMonthlyDelayDifference_v0.2.o LogLikelihoodFunction.o MonthlyDelayDifference.o UsefulFunctions.o $(CFLAGS)

UsefulFunctions.o: std_lib_facilities.h UsefulFunctions.h
			g++ -c UsefulFunctions.cpp

vonMisesRecDist.o: vonMisesRecDist.cxx
	g++ $(CFLAGS) -c vonMisesRecDist.cxx

LogLikelihoodFunction.o: LogLikelihoodFunction.cxx LogLikelihoodFunction.h MonthlyDelayDifference.h 
	g++ $(CFLAGS) -c LogLikelihoodFunction.cxx

MonthlyDelayDifference.o: MonthlyDelayDifference.cxx
	g++ $(CFLAGS) -c MonthlyDelayDifference.cxx

FitMonthlyDelayDifference_v0.2.o: FitMonthlyDelayDifference_v0.2.cxx MonthlyDelayDifference.h LogLikelihoodFunction.h
	g++ $(CFLAGS) -c FitMonthlyDelayDifference_v0.2.cxx

############################################################################################################################################################
clean:
	rm -f *.o FitMonthlyDelayDifference_v0.2 MonthlyDDprojections_v02


