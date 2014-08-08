//g++ -c test.cpp
// g++ -shared -Wl,-soname -o test.so test.o -lc


#include <iostream>

void test(double &par){

  std::cout << par << "\n";
}
