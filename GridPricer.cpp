//
// Created by Albert Wong on 3/21/18.
//

#include "GridPricer.h"


GridPricer::GridPricer(
        std::vector<double> xRanges,
        std::vector<double> timeRanges,
        std::function<double(const double&, const double&)> upperBound,
        std::function<double(const double&, const double&)> lowerBound,
        std::unique_ptr<DiffusionPDE> diffusionPDE
): m_xRanges(std::move(xRanges)),
   m_timeRanges(std::move(timeRanges)),
   m_upperBound(std::move(upperBound)),
   m_lowerBound(std::move(lowerBound)),
   m_diffusionPDE(std::move(diffusionPDE)),
   m_timeAnchor(-1) {
    m_states.resize(m_xRanges.size(), 0.0);
}


void GridPricer::SetPayout(std::function<double(double)> payoff) {
    std::transform(m_xRanges.begin(), m_xRanges.end(), m_states.begin(), std::move(payoff));
    m_timeAnchor = static_cast<int>(m_timeRanges.size() - 2);
}


void GridPricer::StepBack(double timePoint, FdType method, std::function<double(const double&, const double&)> mask) {
    if (m_timeAnchor < 0)
        throw std::runtime_error("SetPayout not called yet");

    // TODO: for hull white, when parameters are time-dependent, need to test if should use t0 or t1 to generate the parameters
    for (; m_timeAnchor >= 0 && m_timeRanges[m_timeAnchor] >= timePoint; --m_timeAnchor) {
        // set up variables
        std::vector<double> pd, pm, pu;
        double &t0 = m_timeRanges[m_timeAnchor], &t1 = m_timeRanges[m_timeAnchor + 1];
        double lb = m_lowerBound(m_xRanges.front(), t0), ub = m_upperBound(m_xRanges.back(), t0);

        // finite difference
        if (method == Explicit || method == CrankNicolson) {
            std::vector<double> temp(m_xRanges.size());
            if (method == Explicit)
                std::tie(pd, pm, pu) = m_diffusionPDE->GetParameters(m_xRanges, t0, t1, Explicit);
            else  // crank nicolson parameter is half of that of explicit
                std::tie(pd, pm, pu) = m_diffusionPDE->GetParameters(m_xRanges, t0, t1, Explicit, 2);

            for(size_t j = 1; j < m_xRanges.size() - 1; ++j)
                temp[j] = pd[j - 1] * m_states[j - 1] + pm[j - 1] * m_states[j] + pu[j - 1] * m_states[j + 1];
            m_states.swap(temp);
        }

        if (method == Implicit || method == CrankNicolson || method == SOR) {
            if (method == CrankNicolson)  // crank nicolson parameter is half of that of implicit
                std::tie(pd, pm, pu) = m_diffusionPDE->GetParameters(m_xRanges, t0, t1, Implicit, 2);
            else
                std::tie(pd, pm, pu) = m_diffusionPDE->GetParameters(m_xRanges, t0, t1, Implicit);

            if (method == SOR) {
                std::vector<double> guess = m_states;
                m_states[1] -= pd[0] * lb;
                m_states.rbegin()[1] -= pu.back() * ub;
                PricerHelper::PrintVector(m_states);
                m_states = PricerHelper::SolveSOR(pd, pm, pu, m_states, guess);
            } else {
                m_states[1] -= pd[0] * lb;
                m_states.rbegin()[1] -= pu.back() * ub;
                m_states = PricerHelper::SolveTridiagonal(pd, pm, pu, m_states, 1);
            }
        }

        if(mask)
            for(auto it1=m_states.begin(), it2=m_xRanges.begin(); it1 != m_states.end(); it1++, it2++)
                *it1 = mask(*it1, *it2);

        m_states.front() = lb;
        m_states.back() = ub;
    }
}


double GridPricer::GetValue(size_t pos) const {
    assert(pos < m_states.size());
    return m_states[pos];
}


void GridPricer::PrintStates() const {
    for(auto x: m_states)
        std::cout << x << ' ';
    std::cout << std::endl;
}