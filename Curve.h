//
// Created by Albert Wong on 3/19/18.
//

#ifndef HULLWHITE_CURVE_H
#define HULLWHITE_CURVE_H

#include <vector>


class Curve {
private:
    std::vector<double> xValues_, yValues_;
    double leftGrad_, rightGrad_;

public:
    Curve(
            std::vector<double> xValues,
            std::vector<double> yValues,
            double leftGrad=0.0,
            double rightGrad=0.0
    ) : xValues_(std::move(xValues)), yValues_(std::move(yValues)), leftGrad_(leftGrad), rightGrad_(rightGrad) {};

    Curve(
            double xValue,
            double yValue,
            double leftGrad=0.0,
            double rightGrad=0.0
    ) : Curve(std::vector<double>(1, xValue), std::vector<double>(1, yValue), leftGrad, rightGrad) {};

    virtual ~Curve() = default;

    virtual double operator()(double x) const;
};


#endif //HULLWHITE_CURVE_H
