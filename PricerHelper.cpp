//
// Created by Albert Wong on 3/23/18.
//

#include "PricerHelper.h"

std::vector<double> PricerHelper::GetLinSpace(const double &start, const double &end, const size_t &steps) {
    std::vector<double> result(steps + 1);
    result[0] = start;
    double delta = (end - start) / steps, curr = start;
    for(auto iter = result.begin() + 1; iter != result.end(); ++iter)
        *iter = curr += delta;
    return result;
}

std::vector<double> PricerHelper::SolveTridiagonal(const std::vector<double> &pd, const std::vector<double> &pm,
                                                   const std::vector<double> &pu, const std::vector<double> &d,
                                                   size_t offset
) {
    assert(pd.size() == pm.size());
    assert(pm.size() == pu.size());

    std::vector<double> result(pm.size() + 2*offset);

    // forward
    m_pu_p[0] = pu[0] / pm[0];
    m_d_p[0] = d[offset] / pm[0];
    for(size_t i = 1; i < pm.size(); ++i) {
        double den = pm[i] - pd[i] * m_pu_p[i - 1];
        m_pu_p[i] = pu[i] / den;
        m_d_p[i] = (d[offset + i] - pd[i] * m_d_p[i - 1]) / den;
    }

    // backward
    result.rbegin()[offset] = m_d_p.back();
    for(auto i = static_cast<int>(pm.size()) - 2; i >= 0; --i) {
        result[i + offset] = m_d_p[i] - m_pu_p[i] * result[i + 1 + offset];
    }

    return result;
}