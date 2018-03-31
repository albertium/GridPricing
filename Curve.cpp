//
// Created by Albert Wong on 3/19/18.
//

#include <vector>
#include "Curve.h"

double Curve::operator()(const double &x, int order) const {
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


CubicSpline::CubicSpline(
        std::vector<double> xValues,
        std::vector<double> yValues
) : m_a(xValues.size()), m_b(xValues.size() - 1), m_d(xValues.size() - 1)
{
    assert(xValues.size() == yValues.size());
    // spline is of the form y = a + bx + cx^2 + dx^3
    const size_t N = yValues.size();
    std::vector<double> h(N), C(N);  // C as in curvature
    std::vector<double> d(N), m(N), u(N);  // tridiagonal

    for(size_t i = 0; i < N; i++) {
        m_a[i] = yValues[i];
        h[i] = xValues[i + 1] - xValues[i];
    }

    for(size_t i = 1; i < N; i++) {
        // set A as in AX=b
        u[i] = h[i];
        m[i] = 2 * (h[i - 1] + h[i]);
        d[i] = h[i - 1];

        // set b as in AX=b
        C[i] = 3 * ((m_a[i + 1] - m_a[i]) / h[i] - (m_a[i] - m_a[i - 1]) / h[i - 1]);
    }
    u.front() = 0;
    m.front() = m.back() = 1;
    d.back() = 0;
    C.front() = C.back() = 0;

    m_c = PricerHelper::SolveTridiagonal(d, m, u, C);
    for(size_t i = 0; i < N; i++) {
        m_b[i] = (m_a[i + 1] - m_a[i]) / h[i] - h[i] * (m_c[i + 1] + 2*m_c[i]) / 3;
        m_d[i] = (m_c[i + 1] - m_c[i]) / 3 / h[i];
    }
    m_Xs.swap(xValues);
}

double CubicSpline::operator()(const double &x, int order) const {
    auto idx = std::max(std::distance(m_Xs.cbegin(), std::upper_bound(m_Xs.cbegin(), m_Xs.cend(), x)) - 1, 0l);
    double dx = x - m_Xs[idx];
    if (order == 0) {
        return m_a[idx] + m_b[idx] * dx + m_c[idx] * dx * dx + m_d[idx] * dx * dx * dx;
    } else if (order == 1) {
        return m_b[idx] + 2 * m_c[idx] * dx + 3 * m_d[idx] * dx * dx;
    }
    return 0.0;
}