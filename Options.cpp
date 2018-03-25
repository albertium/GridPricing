//
// Created by Albert Wong on 3/23/18.
//

#include "Options.h"
#include "GridPricer.h"

#include <algorithm>


EuropeanOption::EuropeanOption(double S,
                               double K,
                               double r,
                               double sig,
                               double T,
                               OptionType type
) : m_S(S), m_K(K), m_r(r), m_sig(sig), m_T(T){

    // define grid
    auto prices = PricerHelper::GetLinSpace(0, 2 * m_S, 200);
    auto dates = PricerHelper::GetLinSpace(0, m_T, 50);
    std::unique_ptr<DiffusionPDE> diffusion(new BlackScholesPDE(r, sig));
    auto upperBound = [K, r](const double &S, const double &t) {return S - K * exp(-r * t); };
    auto lowerBound = [](const double &S, const double &t){return 0; };
    GridPricer pricer(prices, dates, upperBound, lowerBound, std::move(diffusion));

    // define payoff
    std::function<double(double)> payoff;
    if(type == OptionType::Call)
        payoff = Curve(K, 0.0, 0.0, 1.0);
    else if(type == OptionType::Put)
        payoff = Curve(K, 0.0, -1.0,0.0);
    else
        throw std::runtime_error("unrecognized option type");

    // calculate price
    pricer.SetPayout(payoff);
    pricer.StepBack(0);
}


TestOption::TestOption() {
    std::vector<double> prices = PricerHelper::GetLinSpace(-2, 2, 8);
    std::vector<double> times = PricerHelper::GetLinSpace(0, 1, 256);
    std::unique_ptr<HeatPDE> diffusion(new HeatPDE());
    auto upperBound = [](const double &x, const double &t) { return exp(3 - t); };
    auto lowerBound = [](const double &x, const double &t) { return exp(-1 - t); };
    auto payout = [](const double &x) { return exp(x); };
    GridPricer pricer(prices, times, upperBound, lowerBound, std::move(diffusion));

    std::vector<double> exact(prices.size());
    std::transform(prices.begin(), prices.end(), exact.begin(), [](const double &x) { return exp(x + 1); } );

    // test explicit finite difference
    pricer.SetPayout(payout);
    pricer.StepBack(0);
    std::cout << "Test explicit finite difference: "
              << PricerHelper::CalculateRMSE(pricer.GetValues(), exact) << std::endl;

    // test implicit finite difference
    pricer.SetPayout(payout);
    pricer.StepBack(0, Implicit);
    std::cout << "Test implicit finite difference: "
              << PricerHelper::CalculateRMSE(pricer.GetValues(), exact) << std::endl;

    // test tridiagonal solver
    PricerHelper PH(3);
    auto result = PH.SolveTridiagonal({0, 1, 2}, {1, 2, 1}, {1, 3, 0}, {3, 14, 7});
    std::cout << "Test tridiagonal: " << PH.CalculateRMSE(result, {1, 2, 3});
}

double TestOption::GetPrice() {
    return 0;
}