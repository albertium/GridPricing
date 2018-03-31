#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <fstream>

#include "GridPricer.h"
#include "Options.h"

int main() {
    std::vector<double>
            x{0.0, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0},
            y{0.00544859, 0.00697146, 0.00788897, 0.01306847, 0.02447046, 0.03070799, 0.03350761, 0.03243098, 0.02976297,
                0.02966949, 0.02928756, 0.02442422};

    CubicSpline spline(x, y);
    auto vec = PricerHelper::GetLinSpace(0, 30, 200);

    std::ofstream file;
    double dx = 0.00001;
    file.open("/Users/wonalbe/Downloads/curve.csv");
    for(auto &x : vec)
        file << spline(x, 1) << ',' << (spline(x + dx) - spline(x - dx)) / 2 / dx << std::endl;
    file.close();
//    TestOption opt;
//    AmericanOption opt(41, 40, 0.04, 0.02, 0.35, 0.75);
//    auto price = opt.GetPrice();
//    std::cout << price << std::endl << price / 4.083817051176386 << std::endl;
//    auto price = opt.GetPrice(), bs = PricerHelper::GetBlackScholesPrice(100, 100, 0.01, 0.02, 0.2, 1);
//    int n = 8;
//    std::cout << std::left << std::showpoint << std::setprecision(3) << std::fixed;
//    std::cout << std::setw(n) << "BS: " << bs << std::endl;
//    std::cout << std::setw(n) << "Grid: " << price << std::endl;
//    std::cout << std::setw(n) << "error: " << price / bs * 100 - 100 << '%' << std::endl;

    return 0;
}