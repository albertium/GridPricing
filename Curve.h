//
// Created by Albert Wong on 3/19/18.
//

#ifndef HULLWHITE_CURVE_H
#define HULLWHITE_CURVE_H

#include <vector>
#include <algorithm>

#include "PricerHelper.h"


class Curve {
public:
    Curve() = default;
    virtual ~Curve() = default;
    virtual double operator()(const double &x) const = 0;
};


class Spline : public Curve {
private:
    std::vector<std::vector<double>> m_coefs;
    std::vector<double> m_xs;

public:
    Spline(std::vector<double> xs, std::vector<std::vector<double>> coefs);
    ~Spline() override = default;
    double operator()(const double &x) const override;
    Spline GetDerivative(size_t order = 1) const;
    static Spline MakeLinear(std::vector<double> xs, std::vector<double> ys, double left_grad = 0, double right_grad = 0);
    static Spline MakeLinear(double x, double y, double left_grad, double right_grad);
    static Spline MakeCubic(std::vector<double> xs, std::vector<double> ys);
};


#endif //HULLWHITE_CURVE_H
