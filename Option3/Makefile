# CREATED  23 April 2014
# MODIFIED 31 January 2017

# AUTHOR Marco.Kienzle@gmail.com

# PURPOSE build the weekly delay difference model

#CFLAGS = -fopenmp -std=c++0x -O2 -I/usr/include -L/usr/local/lib -static -lMinuit2
#CFLAGS = -std=c++0x -O2 -I/usr/include $(shell root-config --glibs --cflags) -lMinuit2
#CFLAGS = -std=c++0x -O2 -I/usr/include -L/usr/local/lib/root -I/usr/local/include/root -lMinuit2

### Attempt to build a static binary (doesn't work!)
#CFLAGS = -std=c++0x -O2 -I/usr/include -L/usr/local/lib -L/usr/local/lib/root -I/usr/local/include/root -lRoot -lMinuit2 -static
#CFLAGS = -std=c++0x -O2 -I/usr/include -L/usr/local/lib/root -I/usr/local/include/root -lRoot -lMinuit2

CFLAGS = $(shell root-config --cflags) 
#CFLAGS = -pthread -std=c++11 -Wno-deprecated-declarations -m64 -I/usr/local/include
LIBS = $(shell root-config --glibs) -lMinuit2

#############################################################################################################################################
### 1. A program to fit the delay difference model to data and estimate parameters
#############################################################################################################################################

DelayDifference_Option3: FitDelayDifference_Option3.o LogLikelihoodFunction.o DelayDifference.o UsefulFunctions.o
			g++ $(CFLAGS) -o DelayDifference_Option3 /usr/lib/prob.o FitDelayDifference_Option3.o LogLikelihoodFunction.o DelayDifference.o UsefulFunctions.o $(LIBS)

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

# # Projections in this directory are obsolete, replaced by the code in ../Projections
# #############################################################################################################################################
# ### 2. A program to simulate a fishery with the delay difference model based on parameters estimated for this model 
# #############################################################################################################################################

# DD_Option3Projections: /usr/lib/prob.o DDprojections2.o Projections2.o UsefulFunctions.o
# 			g++ -std=c++0x -g -o DD_Option3Projections -I/usr/include /usr/lib/prob.o DDprojections2.o Projections2.o UsefulFunctions.o

# DDprojections2.o: DDprojections2.cpp
# 		g++ $(CFLAGS) -c DDprojections2.cpp

# Projections2.o: Projections2.cpp
# 		g++ $(CFLAGS) -c Projections2.cpp

# #############################################################################################################################################
# ### 3. A program to simulate a fishery with the delay difference model based on parameters estimated for this model for MDP (Markov Decision Processes)
# #############################################################################################################################################

# DD_Option3Projections_forMDP: /usr/lib/prob.o DDprojections2_forMDP.o Projections2_forMDP.o UsefulFunctions.o
# 			g++ -std=c++0x -g -o DD_Option3Projections_forMDP -I/usr/include /usr/lib/prob.o DDprojections2_forMDP.o Projections2_forMDP.o UsefulFunctions.o

# DDprojections2_forMDP.o: DDprojections2_forMDP.cpp
# 			g++ $(CFLAGS) -c DDprojections2_forMDP.cpp

# Projections2_forMDP.o: Projections2_forMDP.cpp
# 			g++ $(CFLAGS) -c Projections2_forMDP.cpp

# #############################################################################################################################################
# ### 4. A program to simulate a fishery with the delay difference model based on parameters estimated for this model for MDP (Markov Decision Processes) that takes effort limit as input
# #############################################################################################################################################

# DD_Option3Projections_forMDPWithEffortLimits: /usr/lib/prob.o DDprojections2_forMDPWithEffortLimits.o Projections2_forMDPWithEffortLimits.o UsefulFunctions.o
# 			g++ -std=c++0x -g -o DD_Option3Projections_forMDPWithEffortLimits -I/usr/include /usr/lib/prob.o DDprojections2_forMDPWithEffortLimits.o Projections2_forMDPWithEffortLimits.o UsefulFunctions.o

# DDprojections2_forMDPWithEffortLimits.o: DDprojections2_forMDPWithEffortLimits.cpp
# 			g++ $(CFLAGS) -c DDprojections2_forMDPWithEffortLimits.cpp

# Projections2_forMDPWithEffortLimits.o: Projections2_forMDPWithEffortLimits.cpp
# 			g++ $(CFLAGS) -c Projections2_forMDPWithEffortLimits.cpp

# #############################################################################################################################################
# ### 5. A program to simulate a fishery with the delay difference model based on parameters estimated for this model for MDP (Markov Decision Processes) after the idea of Jerzy Filar and Iadine Chades
# #############################################################################################################################################

# DD_Option3Projections_forMDPNumber2: /usr/lib/prob.o DDprojections2_forMDPNumber2.o Projections2_forMDPNumber2.o UsefulFunctions.o
# 			g++ -std=c++0x -g -o DD_Option3Projections_forMDPNumber2 -I/usr/include /usr/lib/prob.o DDprojections2_forMDPNumber2.o Projections2_forMDPNumber2.o UsefulFunctions.o

# DDprojections2_forMDPNumber2.o: DDprojections2_forMDPNumber2.cpp
# 			g++ $(CFLAGS) -c DDprojections2_forMDPNumber2.cpp

# Projections2_forMDPNumber2.o: Projections2_forMDPNumber2.cpp
# 			g++ $(CFLAGS) -c Projections2_forMDPNumber2.cpp

# #############################################################################################################################################
# ### 5. A program to simulate a fishery with the delay difference model based on parameters estimated for this model for testing
# #############################################################################################################################################

# DD_Option3Projections_fortest: /usr/lib/prob.o DDprojections2_fortest.o Projections2_fortest.o UsefulFunctions.o
# 			g++ -std=c++0x -g -o DD_Option3Projections_fortest -I/usr/include /usr/lib/prob.o DDprojections2_fortest.o Projections2_fortest.o UsefulFunctions.o

# DDprojections2_fortest.o: DDprojections2_fortest.cpp
# 			g++ $(CFLAGS) -c DDprojections2_fortest.cpp

# Projections2_fortest.o: Projections2_fortest.cpp
# 			g++ $(CFLAGS) -c Projections2_fortest.cpp

#############################################################################################################################################
### make all
#############################################################################################################################################

#all: DelayDifference_Option3 DD_Option3Projections DD_Option3Projections_forMDP DD_Option3Projections_forMDPNumber2 DD_Option3Projections_fortest DD_Option3Projections_forMDPWithEffortLimits
all: DelayDifference_Option3

#############################################################################################################################################
### make clean
#############################################################################################################################################

clean:
	rm -f *.o DelayDifference_Option3

############################################################################################################################################
### make test
###########################################################################################################################################
testDelayDifference:
	./DelayDifference_Option3 Data/SimData4.txt InputParameterDescription.csv
	#diff Results/FitOutcome.txt Results/FitOutcome-2016Apr6.txt
	diff Results/ParameterEstimates.txt Results/ParameterEstimates-2016Apr6.txt


