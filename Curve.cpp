//
// Created by Albert Wong on 3/19/18.
//

#include <vector>
#include "Curve.h"

double Curve::operator()(double x) const {
    int idx = 0;
    while(xValues_[idx] < x && idx < xValues_.size())
        idx++;

    if(idx == 0)
        return yValues_[0] + leftGrad_ * (x - xValues_[0]);

    if(idx >= xValues_.size())
        return yValues_.back() + rightGrad_ * (x - xValues_.back());

    return yValues_[idx - 1] +
            (yValues_[idx] - yValues_[idx - 1]) / (xValues_[idx] - xValues_[idx - 1]) * (x - xValues_[idx - 1]);
}