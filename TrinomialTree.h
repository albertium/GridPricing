//
// Created by Albert Wong on 3/19/18.
//

#ifndef HULLWHITE_TRINOMIALTREE_H
#define HULLWHITE_TRINOMIALTREE_H

#include <vector>
#include <functional>
#include "Curve.h"

class TrinomialTree {
private:
    Curve zeroRates_;
    std::function<double(double)> payoff_;
    std::vector<std::vector<double>> grid_, pu_, pm_, pd_;
    std::vector<double> zcBonds_;
    double dR_, dt_, T_, a_, sig_;
    int N_, jmax_, jmin_;

    void BuildTree();

public:
    explicit TrinomialTree(
            Curve& zeroRates,
            double T,
            int N
    ) : zeroRates_(std::move(zeroRates)), T_(T), N_(N) {}

    void SetParameters(double a, double sig);
};


#endif //HULLWHITE_TRINOMIALTREE_H
