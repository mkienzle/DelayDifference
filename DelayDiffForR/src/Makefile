CFLAGS = -std=c++0x -O2 -I/usr/include $(shell root-config --libs --cflags --ldflags) -lMinuit2
#CFLAGS = -fopenmp -std=c++0x -O2 -I/usr/include -L/usr/local/lib -static -lMinuit2
PKG_CXXFLAGS=-I/usr/local/lib/R/include -DNDEBUG  -I/usr/local/include   -I/home/mkienzle/R/x86_64-unknown-linux-gnu-library/3.2/Rcpp/include -fpic -g -O2
#CFLAGS = -fopenmp -std=c++0x -O2 -I/usr/include -L/usr/local/lib -static -lMinuit2

#FitLW: ObjectiveFCN.o FitLW.o
#	g++ -o FitLW ObjectiveFCN.o FitLW.o $(CFLAGS)
#FitLW.o: FitLW.cpp lib_facilities.h
#	g++ $(CFLAGS) -c FitLW.cpp

DelayDiffForR.so: FitLW.o ObjectiveFCN.o
	 g++ -shared -o DelayDiffForR.so FitLW.o ObjectiveFCN.o $(shell ls /home/mkienzle/mystuff/Software/root-6.06.06/math/minuit2/src/*.o) $(CFLAGS) 
FitLW.o: FitLW.cpp ObjectiveFCN.o
	g++ $(PKG_CXXFLAGS) $(CFLAGS) -o FitLW.o -c FitLW.cpp
ObjectiveFCN.o: ObjectiveFCN.cpp ObjectiveFCN.h
	g++ $(PKG_CXXFLAGS) $(CFLAGS) -c ObjectiveFCN.cpp
