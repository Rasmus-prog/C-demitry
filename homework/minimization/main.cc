#include "sfuns.h"

#include <iomanip>
#include <iostream>
#include <string>
#include <vector>


int main() {
    std::cout << std::setprecision(12);

    reportSolution("Rosenbrock", rosenbrock, std::vector<double>{-1.2, 1.0});
    reportSolution("Himmelblau", himmelblau, std::vector<double>{3.1, 2.1});
    
    return 0;
}
