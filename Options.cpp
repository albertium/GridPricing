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
        payoff = Spline::MakeLinear(K, 0.0, 0.0, 1.0);
    else if(type == OptionType::Put)
        payoff = Spline::MakeLinear(K, 0.0, -1.0, 0.0);
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
              << PricerHelper::CalculateRMSE(pricer.GetValues(), exact) << " (0.0784)" <<  std::endl;

    // test implicit finite difference
    pricer.SetPayout(payout);
    pricer.StepBack(0, Implicit);
    std::cout << "Test implicit finite difference: "
              << PricerHelper::CalculateRMSE(pricer.GetValues(), exact) << " (0.1517)" << std::endl;

    // test Crank Nicolson finite difference
    pricer.SetPayout(payout);
    pricer.StepBack(0, CrankNicolson);
    std::cout << "Test Crank Nicolson finite difference: "
              << PricerHelper::CalculateRMSE(pricer.GetValues(), exact) << " (0.0407)" <<  std::endl;

    // test American option
    AmericanOption opt(41, 40, 0.04, 0.02, 0.35, 0.75);
    auto price = opt.GetPrice();
    std::cout << "Test American option: " << abs(price / 4.083817051176386 - 1) << std::endl;

    // test linear spline
    Spline payoff = Spline::MakeLinear(5, 0, -1, 2);
    std::vector<double> vec{2, 3, 4, 5, 6, 7};
    for (auto &x : vec)
        x = payoff(x);
    std::cout << "Test linear piecewise: " << PricerHelper::CalculateRMSE(vec, {3, 2, 1, 0, 2, 4}) << std::endl;

    // test cubic spline
    auto xs = PricerHelper::GetLinSpace(-3, 3, 10);
    std::vector<double> ys(xs.size()), vec2;
    std::transform(xs.begin(), xs.end(), ys.begin(), [](double x) { return x * x * x; });
    Spline cubic = Spline::MakeCubic(xs, ys);
    vec.resize(xs.size());
    std::transform(xs.begin(), xs.end(), vec.begin(), cubic);
    std::cout << "Test cubic: " << PricerHelper::CalculateRMSE(vec, ys) << std::endl;

    // test cubic derivatives
    Spline deriv = cubic.GetDerivative(), deriv2 = cubic.GetDerivative(2);
    xs = PricerHelper::GetLinSpace(-4, 4, 50);
    double dx = 0.000001;
    vec.resize(50);
    vec2.resize(50);
    std::transform(xs.begin(), xs.end(), vec.begin(), deriv);
    std::transform(xs.begin(), xs.end(), vec2.begin(),
                   [cubic, dx](double x ) { return (cubic(x + dx) - cubic(x - dx)) / 2 / dx; });
    std::cout << "Test cubic 1st derivative: " << PricerHelper::CalculateRMSE(vec, vec2) << std::endl;

    std::transform(xs.begin(), xs.end(), vec.begin(), deriv2);
    std::transform(xs.begin(), xs.end(), vec2.begin(),
                   [cubic, dx](double x ) { return (cubic(x + dx) - 2*cubic(x) + cubic(x - dx)) / dx / dx; });
    std::cout << "Test cubic 2nd derivative: " << PricerHelper::CalculateRMSE(vec, vec2) << std::endl;

    // test Ax=Bd
    auto res = PricerHelper::SolveTridiagonalAxBd({}, {{1, 2, 3}, {3, 4, 5}, {5, 6, 7}}, {1, 2, 3}, {}, {1, 2});
    std::cout << "Test Ax=Bd (1): " << PricerHelper::CalculateRMSE(res, {9, 26, 42}) << std::endl;
    res = PricerHelper::SolveTridiagonalAxBd({}, {{1, 2, 3}, {3, 4, 5}, {5, 6, 7}}, {1, 2, 3});
    std::cout << "Test Ax=Bd (2): " << PricerHelper::CalculateRMSE(res, {14, 26, 38}) << std::endl;
    res = PricerHelper::SolveTridiagonalAxBd({{1, 2, 3}, {3, 4, 5}, {5, 6, 7}}, {}, {9, 26, 42}, {1, 2});
    std::cout << "Test Ax=Bd (3): " << PricerHelper::CalculateRMSE(res, {1, 2, 3}) << std::endl;
    res = PricerHelper::SolveTridiagonalAxBd({{1, 2, 1}, {3, 4, 5}, {5, 6, 7}}, {}, {8, 26, 38});
    std::cout << "Test Ax=Bd (4): " << PricerHelper::CalculateRMSE(res, {1, 2, 3}) << std::endl;


    // test tridiagonal solver
    auto result = PricerHelper::SolveTridiagonal({0, 1, 2}, {1, 2, 1}, {1, 3, 0}, {3, 14, 7});
    std::cout << "Test tridiagonal: " << PricerHelper::CalculateRMSE(result, {1, 2, 3}) << std::endl;

    // test SOR solver
    result = PricerHelper::SolveSOR({1, 1}, {3, 3}, {1, 1}, {10, 15}, {1, 0, 0, 4});
    std::cout << "Test SOR: " << PricerHelper::CalculateRMSE(result, {1, 2, 3, 4}) << std::endl;
}