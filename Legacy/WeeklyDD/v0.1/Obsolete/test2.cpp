#include <stdio.h>
#include <stdlib.h>
#include <iostream>
 
int main(int argc, char *argv[]) {
    double n = strtod(argv[1],NULL);
 
    printf("The number: %f\n",n);
 
    std::cout << "And the second arg is " << argv[2] << "\n";
    std::cout << "And the second arg is " << argv[3] << "\n";
    std::cout << "And the second arg is " << argv[4] << "\n";
    std::cout << "And the second arg is " << argv[5] << "\n";
    std::cout << "And the second arg is " << argv[6] << "\n";
    std::cout << "And the second arg is " << argv[7] << "\n";
    std::cout << "And the second arg is " << argv[8] << "\n";


    return 0;
}
