//
// Created by Albert Wong on 3/23/18.
//

#ifndef HULLWHITE_PRICERHELPER_H
#define HULLWHITE_PRICERHELPER_H

#include <vector>
#include <iostream>
#include <cassert>
#include <cmath>

#include "Types.h"


class PricerHelper {
private:
    std::vector<double> m_pu_p, m_d_p, m_temp;

public:
    explicit PricerHelper(size_t size) : m_pu_p(size), m_d_p(size), m_temp(size) {};
    virtual ~PricerHelper() = default;

    // Linear algebra
    static std::vector<double> GetLinSpace(const double &start, const double &end, const size_t &steps);

    std::vector<double> SolveTridiagonal(
            const std::vector<double>& pd,
            const std::vector<double>& pm,
            const std::vector<double>& pu,
            const std::vector<double>& d,
            size_t offset = 0
    );  // will pad offset 0s in front and in the end of the returned vector

    static std::vector<double> SolveSOR(
            const std::vector<double>& pd,
            const std::vector<double>& pm,
            const std::vector<double>& pu,
            const std::vector<double>& d,
            std::vector<double> guess,
            const std::function<double(double)> &filter = {},
            const double &omega = 0.5,
            const double &tol = 0.0001
    );

    // Formula
    inline static double NormalCDF(double x);
    static double GetBlackScholesPrice(double S, double K, double r, double q, double sig, double T);

    // Misc
    template<typename T>
    static void PrintVector(const std::vector<T>& toPrint);

    template<typename T>
    static double CalculateRMSE(const std::vector<T> &x, const std::vector<T> &y);
};


double PricerHelper::NormalCDF(double x) {
    return erfc(-x * M_SQRT1_2) / 2;
}

template<typename T>
void PricerHelper::PrintVector(const std::vector<T>& toPrint) {
    for(auto &x : toPrint)
        std::cout << x << ' ';
    std::cout << std::endl;
}


template<typename T>
double PricerHelper::CalculateRMSE(const std::vector<T> &x, const std::vector<T> &y) {
    assert(x.size() == y.size());
    T error = 0;
    for(size_t i = 0; i < x.size(); ++i)
        error += pow(x[i] - y[i], 2);
    return sqrt(error / x.size());
}

#endif //HULLWHITE_PRICERHELPER_H
