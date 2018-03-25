//
// Created by Albert Wong on 3/19/18.
//

#include "TrinomialTree.h"
#include <cmath>

void TrinomialTree::SetParameters(double a, double sig) {
    a_ = a;
    sig_ = sig;
    BuildTree();
}

void TrinomialTree::BuildTree() {
    dt_ = T_ / N_;
    dR_ = sig_ * sqrt(3 * dt_);
    jmax_ = static_cast<int>(ceil(0.184 / a_));
    jmin_ = -jmax_;

    // initialize grid and probabilities
    for(size_t i = 0; i < N_; ++i) {
        grid_.emplace_back(N_, 0.0);
        pu_.emplace_back(N_, 0.0);
        pm_.emplace_back(N_, 0.0);
        pd_.emplace_back(N_, 0.0);
    }

    // stage 1

}