// CREATED   4 February 2014
// MODIFIED  4 February 2014

// PURPOSE given 2 parameters (mu and kappa), display the density distribution in each month

// COMPILE g++ -std=c++0x -g -o DisplayVonMisesMontlyDistribution -I/usr/include /usr/lib/prob.o vonMisesRecDist.cxx DisplayVonMisesMontlyDistribution.cxx

# include <math.h>
# include <vector>
# include <iostream>
# define NMPY 12

std::vector<double> vonMisesRecDist(double a, double b);

int main(){

double mu=0.0, kappa=0.0;

std::cout << "Program to display the von mises distribution\n";

std::cout << "Enter parameter mu \n";
 std::cin >> mu;

std::cout << "Enter parameter kappa \n";
 std::cin >> kappa;

  // Calculate the proportion of recruitment in each month
  std::vector<double> RecDist(NMPY, 0.0);
  RecDist = vonMisesRecDist(mu, kappa);

  for(unsigned int counter = 0; counter < NMPY; ++counter)
    std::cout << RecDist.at(counter) << ",";

  std::cout << "\n";


 return 0;}
