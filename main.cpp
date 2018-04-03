#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <fstream>

#include "GridPricer.h"

int main() {
//    std::vector<double>
//            x{0.0, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0},
//            y{0.00544859, 0.00697146, 0.00788897, 0.01306847, 0.02447046, 0.03070799, 0.03350761, 0.03243098, 0.02976297,
//                0.02966949, 0.02928756, 0.02442422};

//    TestOption opt;
    auto res = PricerHelper::SolveTridiagonalAxBd({{1, 2, 1}, {3, 4, 5}, {5, 6, 7}}, {}, {8, 26, 38});
    for(auto &x : res )
        std::cout << x << " ";
    std::cout << std::endl;

//    std::cout << std::left << std::showpoint << std::setprecision(3) << std::fixed;
//    std::cout << std::setw(n) << "BS: " << bs << std::endl;
//    std::cout << std::setw(n) << "Grid: " << price << std::endl;
//    std::cout << std::setw(n) << "error: " << price / bs * 100 - 100 << '%' << std::endl;

    return 0;
}