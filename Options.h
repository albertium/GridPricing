//
// Created by Albert Wong on 3/23/18.
//

#ifndef HULLWHITE_OPTIONS_H
#define HULLWHITE_OPTIONS_H

#include "PricerHelper.h"
#include "DiffusionPDE.h"
#include "GridPricer.h"
#include "Curve.h"

#include <cmath>
#include <algorithm>


class Options {
protected:
    double m_price;

public:
    explicit Options() : m_price(0.0) {};
    virtual ~Options() = default;
    inline double GetPrice() const { return m_price; };
};


class EuropeanOption : public Options {
private:
    double m_S, m_K, m_r, m_q, m_sig, m_T;

public:
    explicit EuropeanOption(double S, double K, double r, double q, double sig, double T, OptionType type);
    ~EuropeanOption() override = default;
};


class AmericanOption : public Options {
public:
    explicit AmericanOption(double S, double K, double r, double q, double sig, double T);  // put option
    ~AmericanOption() override = default;
};


class TestOption : public Options {
public:
    explicit TestOption();
    ~TestOption() override = default;
};

#endif //HULLWHITE_OPTIONS_H
