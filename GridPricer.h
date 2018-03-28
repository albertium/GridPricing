//
// Created by Albert Wong on 3/21/18.
//

#ifndef HULLWHITE_GRIDPRICER_H
#define HULLWHITE_GRIDPRICER_H

#include <vector>
#include <functional>
#include <iostream>
#include <algorithm>
#include <exception>
#include <memory>

#include "DiffusionPDE.h"
#include "PricerHelper.h"
#include "Types.h"


class GridPricer {
private:
    std::vector<double> m_xRanges, m_timeRanges, m_states;  // m_states represents the contract value
    std::function<double(const double&, const double&)> m_upperBound, m_lowerBound;
    std::unique_ptr<DiffusionPDE> m_diffusionPDE;
    int m_timeAnchor;

public:
    explicit GridPricer(
            std::vector<double> xRanges,
            std::vector<double> timeRanges,
            std::function<double(const double&, const double&)> upperBound,
            std::function<double(const double&, const double&)> lowerBound,
            std::unique_ptr<DiffusionPDE> diffusionPDE
    );

    virtual ~GridPricer() = default;

    void StepBack(double timePoint=0, FdType method=Explicit, std::function<double(const double&, const double&)> mask = {});
    void SetPayout(std::function<double(double)> payoff);
    double GetValue(size_t pos) const;
    inline const std::vector<double>& GetValues() const { return m_states; };
    void PrintStates() const;
};


#endif //HULLWHITE_GRIDPRICER_H
