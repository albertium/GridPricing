//
// Created by Albert Wong on 3/22/18.
//

#include "DiffusionPDE.h"

std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> DiffusionPDE::GetExplicitParameters(
        const std::vector<double> &xRanges,
        const double &t1,  // t(i)
        const double &t2   // t(i+1)
) const {
    std::vector<double> pd(xRanges.size() - 2, 0.0), pm(xRanges.size() - 2, 0.0), pu(xRanges.size() - 2, 0.0);
    double dt = t2 - t1;
    for(size_t j = 1; j < xRanges.size() - 1; ++j) {
        const double &Xu = xRanges[j + 1], &Xm = xRanges[j], &Xd = xRanges[j - 1];
        double dX = Xu - Xd;
        pu[j - 1] = 2 * dt * (GetDiffusion(Xm, t1) / (Xm - Xd) / dX - GetDrift(Xm, t1) / 2 / dX);
        pm[j - 1] = 1 - 2 * dt * (GetDiffusion(Xm, t1) / (Xu - Xm) / (Xm - Xd) - GetInterest(Xm, t1));
        pd[j - 1] = 2 * dt * (GetDiffusion(Xm, t1) / (Xu - Xm) / dX + GetDrift(Xm, t1) / 2 / dX);
    }
    return std::make_tuple(std::move(pd), std::move(pm), std::move(pu));
}


std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> DiffusionPDE::GetImplicitParameters(
        const std::vector<double> &xRanges,
        const double &t1,
        const double &t2) const {
    // still use the convention of starting at 0 instead of 1 because we don't need to do any adjustment later
    // when using SolveTridiagonal
    std::vector<double> pd(xRanges.size() - 2, 0.0), pm(xRanges.size() - 2, 0.0), pu(xRanges.size() - 2, 0.0);
    double dt = t2 - t1;
    for(size_t j = 1; j < xRanges.size() -1; ++j) {
        const double &Xd = xRanges[j - 1], &Xm = xRanges[j], &Xu = xRanges[j + 1];
        double dX = Xu - Xd, mu = GetDrift(Xm, t1), sig2 = GetDiffusion(Xm, t1);
        pd[j - 1] = dt * (-2 * sig2 / (Xm - Xd) + mu) / dX;
        pm[j - 1] = 1 + dt * (GetInterest(Xm, t1) + 2 * sig2 / (Xm - Xd) / (Xu - Xm));
//        pm[j - 1] = 1;
        pu[j - 1] = -dt * (2 * sig2 / (Xu - Xm) + mu) / dX;
    }
    return std::make_tuple(std::move(pd), std::move(pm), std::move(pu));
}