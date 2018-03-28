#include <iostream>
#include <iomanip>
#include <vector>

#include "GridPricer.h"
#include "Options.h"

int main() {
//    TestOption opt;
    AmericanOption opt(41, 40, 0.04, 0.02, 0.35, 0.75);
    auto price = opt.GetPrice();
    std::cout << price << std::endl << price / 4.083817051176386 << std::endl;
//    auto price = opt.GetPrice(), bs = PricerHelper::GetBlackScholesPrice(100, 100, 0.01, 0.02, 0.2, 1);
//    int n = 8;
//    std::cout << std::left << std::showpoint << std::setprecision(3) << std::fixed;
//    std::cout << std::setw(n) << "BS: " << bs << std::endl;
//    std::cout << std::setw(n) << "Grid: " << price << std::endl;
//    std::cout << std::setw(n) << "error: " << price / bs * 100 - 100 << '%' << std::endl;

    return 0;
}