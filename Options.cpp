//
// Created by Albert Wong on 3/23/18.
//

#include "Options.h"


EuropeanOption::EuropeanOption(double S,
                               double K,
                               double r,
                               double q,
                               double sig,
                               double T,
                               OptionType type
) : m_S(S), m_K(K), m_r(r), m_q(q), m_sig(sig), m_T(T){

    // define grid
    const int N = 1000;
    auto prices = PricerHelper::GetLinSpace(0, 2 * m_S, N);
    auto dates = PricerHelper::GetLinSpace(0, m_T, 1000);
    std::unique_ptr<DiffusionPDE> diffusion(new BlackScholesPDE(r, q, sig));
    auto upperBound = [K, r, T](const double &S, const double &t) {return S - K * exp(-r* (T - t)); };
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
    pricer.StepBack(0, CrankNicolson);
    m_price = pricer.GetValue(static_cast<size_t>(N / 2));
}


AmericanOption::AmericanOption(double S, double K, double r, double q, double sig, double T) {
    const int N = 1000, steps = 1000;
    auto prices = PricerHelper::GetLinSpace(0, 2 * S, N);
    auto times = PricerHelper::GetLinSpace(0, T, steps);
    std::unique_ptr<DiffusionPDE> diffusion(new BlackScholesPDE(r, q, sig));
    auto upperBound = [](const double &S, const double &t) { return 0; };
    auto lowerBound = [K](const double &S, const double &t) { return K; };
    GridPricer pricer(prices, times, upperBound, lowerBound, std::move(diffusion));

    // payoff
    auto payoff = [K](const double &S) { return std::max(K - S, 0.0); };
    auto mask = [K](const double &u, const double& S) { return std::max(u, K - S); };

    // stepping
    pricer.SetPayout(payoff);
    pricer.StepBack(0, CrankNicolson, mask);
    m_price = pricer.GetValue(static_cast<size_t>(N / 2));
}


TestOption::TestOption() {
    std::vector<double> prices = PricerHelper::GetLinSpace(-2, 2, 8);
    std::vector<double> times = PricerHelper::GetLinSpace(0, 1, 8);
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

    // test Crank Nicolson finite difference
    pricer.SetPayout(payout);
    pricer.StepBack(0, CrankNicolson);
    std::cout << "Test Crank Nicolson finite difference: "
              << PricerHelper::CalculateRMSE(pricer.GetValues(), exact) << std::endl;

    // test tridiagonal solver
    auto result = PricerHelper::SolveTridiagonal({0, 1, 2}, {1, 2, 1}, {1, 3, 0}, {3, 14, 7});
    std::cout << "Test tridiagonal: " << PricerHelper::CalculateRMSE(result, {1, 2, 3}) << std::endl;

    // test SOR solver
    result = PricerHelper::SolveSOR({1, 1}, {3, 3}, {1, 1}, {10, 15}, {1, 0, 0, 4});
    std::cout << "Test SOR: " << PricerHelper::CalculateRMSE(result, {1, 2, 3, 4}) << std::endl;
}