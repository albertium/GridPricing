//
// Created by Albert Wong on 3/22/18.
//

#include "DiffusionPDE.h"


std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> DiffusionPDE::GetParameters(
        const std::vector<double> &xRanges,
        const double &t1,
        const double &t2,
        FdType method,
        double offset
) const {
    std::vector<double> pd(xRanges.size() - 2, 0.0), pm(xRanges.size() - 2, 0.0), pu(xRanges.size() - 2, 0.0);
    double dt = t2 - t1, sign;
    if (method == Explicit)
        sign = 1;
    else if (method == Implicit)
        sign = -1;
    else
        throw std::runtime_error("method not supported by GetParameters");

    for(size_t j = 1; j < xRanges.size() -1; ++j) {
        const double &Xd = xRanges[j - 1], &Xm = xRanges[j], &Xu = xRanges[j + 1];
        double dX = Xu - Xd, mu = GetDrift(Xm, t1), sig2 = GetDiffusion(Xm, t1);
        pd[j - 1] = sign * dt * (2 * sig2 / (Xm - Xd) - mu) / dX;
        pm[j - 1] = offset - sign * dt * (GetInterest(Xm, t1) + 2 * sig2 / (Xm - Xd) / (Xu - Xm));
        pu[j - 1] = sign * dt * (2 * sig2 / (Xu - Xm) + mu) / dX;
    }
    return std::make_tuple(std::move(pd), std::move(pm), std::move(pu));
}