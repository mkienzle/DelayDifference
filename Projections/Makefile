# CREATED  25 Aug 2017
# MODIFIED 18 Nov 2020

# AUTHOR Marco.Kienzle@gmail.com

# PURPOSE compile the programs to project a delay-difference model

CFLAGS = -pthread -std=c++11 -Wno-deprecated-declarations -m64

#############################################################################################################################################
### 1. 
#############################################################################################################################################
#50YearsProjectionsForMAY: DelayDifferenceModel.o DelayDifferenceModelProjection.o UsefulFunctions.o 50YearsProjectionsForMAY.cpp
#			g++ $(CFLAGS) -o 50YearsProjectionsForMAY 50YearsProjectionsForMAY.cpp /usr/lib/prob.o DelayDifferenceModel.o DelayDifferenceModelProjection.o UsefulFunctions.o 

50YearsProjectionsWithEffortBetween: DelayDifferenceModel.o DelayDifferenceModelProjection.o UsefulFunctions.o 50YearsProjectionsWithEffortBetween.cpp
			g++ $(CFLAGS) -o 50YearsProjectionsWithEffortBetween 50YearsProjectionsWithEffortBetween.cpp /usr/lib/prob.o DelayDifferenceModel.o DelayDifferenceModelProjection.o UsefulFunctions.o 

MDPprojectionType1: DelayDifferenceModel.o DelayDifferenceModelProjection.o UsefulFunctions.o MDPprojectionType1.cpp 
			g++ $(CFLAGS) -o MDPprojectionType1 MDPprojectionType1.cpp /usr/lib/prob.o DelayDifferenceModel.o DelayDifferenceModelProjection.o UsefulFunctions.o 

Tutorial: DelayDifferenceModel.o DelayDifferenceModelProjection.o UsefulFunctions.o Tutorial.cpp 
			g++ $(CFLAGS) -o Tutorial Tutorial.cpp /usr/lib/prob.o DelayDifferenceModel.o DelayDifferenceModelProjection.o UsefulFunctions.o $(LIBS)

Squid_150YearsProjectionWithEffort: DelayDifferenceModel.o DelayDifferenceModelProjection.o UsefulFunctions.o Squid_150YearsProjectionWithEffort.cpp
			g++ $(CFLAGS) -o Squid_150YearsProjectionWithEffort Squid_150YearsProjectionWithEffort.cpp /usr/lib/prob.o DelayDifferenceModel.o DelayDifferenceModelProjection.o UsefulFunctions.o $(LIBS)
UsefulFunctions.o: ../lib_facilities2.h ../UsefulFunctions.h ../UsefulFunctions.cpp
			g++ $(CFLAGS) -c ../UsefulFunctions.cpp

DelayDifferenceModel.o: DelayDifferenceModel.h DelayDifferenceModel.cpp ../UsefulFunctions.cpp ../UsefulClasses.h
	g++ $(CFLAGS) -c DelayDifferenceModel.cpp

DelayDifferenceModelProjection.o: DelayDifferenceModelProjection.h DelayDifferenceModelProjection.cpp ../UsefulFunctions.cpp ../UsefulClasses.h
	g++ $(CFLAGS) -c DelayDifferenceModelProjection.cpp

#############################################################################################################################################
### make all
#############################################################################################################################################
all: 50YearsProjectionsWithEffortBetween MDPprojectionType1 Tutorial Squid_150YearsProjectionWithEffort
#all: 50YearsProjectionsWithEffortBetween 50YearsProjectionsForMAY MDPprojectionType1 Tutorial 
#all: DelayDifference_Option3 DD_Option3Projections DD_Option3Projections_forMDP DD_Option3Projections_forMDPNumber2 DD_Option3Projections_fortest DD_Option3Projections_forMDPWithEffortLimits

#############################################################################################################################################
### make clean
#############################################################################################################################################

clean:
	rm -f *.o 




