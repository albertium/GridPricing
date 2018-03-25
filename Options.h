//
// Created by Albert Wong on 3/23/18.
//

#ifndef HULLWHITE_OPTIONS_H
#define HULLWHITE_OPTIONS_H

#include "PricerHelper.h"
#include "DiffusionPDE.h"
#include "Curve.h"

#include <cmath>


enum OptionType {
    Call,
    Put
};


class Options {
public:
    explicit Options() = default;
    virtual ~Options() = default;
    virtual double GetPrice() = 0;
};


class EuropeanOption : public Options {
private:
    double m_S, m_K, m_r, m_sig, m_T;
public:
    explicit EuropeanOption(double S, double K, double r, double sig, double T, OptionType type=OptionType::Call);
    ~EuropeanOption() override = default;
    double GetPrice() override { return 0.0; };
};


class TestOption : public Options {
public:
    explicit TestOption();
    ~TestOption() override = default;
    double GetPrice() override;
};

#endif //HULLWHITE_OPTIONS_H
