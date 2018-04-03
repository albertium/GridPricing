//
// Created by Albert Wong on 3/19/18.
//

#include "Curve.h"

Spline::Spline(std::vector<double> xs, std::vector<std::vector<double>> coefs) {
    size_t N = coefs[0].size();
    for (auto &x : coefs)
        if (x.size() != N)
            throw std::runtime_error("non-matching coefficients size");

    if (xs.size() != coefs[0].size() - 1)
        throw std::runtime_error("non-matching number of segments");

    m_xs.swap(xs);
    m_coefs.swap(coefs);
}


double Spline::operator()(const double &x) const {
    long idx = std::distance(m_xs.cbegin(), std::upper_bound(m_xs.cbegin(), m_xs.cend(), x));
    double dx = x - m_xs[std::max(idx - 1, 0L)];
    double basis = 1, result = 0.0;
    for(const auto &coef : m_coefs) {
        result += coef[idx] * basis;
        basis *= dx;
    }
    return result;
}


Spline Spline::GetDerivative(size_t order) const {
    assert(order <= m_coefs.size());
    assert(order > 0);

    const size_t N = m_coefs[0].size();
    std::vector<std::vector<double>> newCoefs;

    for (size_t i = order; i < m_coefs.size(); i++) {
        newCoefs.emplace_back(N);
        double mult = 1;
        for (int j = 1; j <= order; j++)
            mult *= (i - order + j);
        for(size_t j = 0; j < N; j++)
            newCoefs.back()[j] = mult * m_coefs[i][j];
    }

    return Spline(m_xs, newCoefs);
}


Spline Spline::MakeLinear(std::vector<double> xs, std::vector<double> ys, double left_grad, double right_grad) {
    assert(xs.size() == ys.size());
    size_t N = xs.size();
    std::vector<double> a_tmp(N + 1), b_tmp(N + 1);

    // left boundary spline
    a_tmp[0] = ys[0];
    b_tmp[0] = left_grad;

    // right boundary spline
    a_tmp.back() = ys.back();
    b_tmp.back() = right_grad;

    for (size_t i = 1; i < N; i++) {
        b_tmp[i] = (ys[i] - ys[i - 1]) / (xs[i] - xs[i - 1]);
        a_tmp[i] = ys[i] - b_tmp[i] * xs[i];
    }

    return Spline(std::move(xs), {std::move(a_tmp), std::move(b_tmp)});
}


Spline Spline::MakeLinear(double x, double y, double left_grad, double right_grad) {
    return Spline::MakeLinear(std::vector<double>{x}, std::vector<double>{y}, left_grad, right_grad);
}


Spline Spline::MakeCubic(std::vector<double> xs, std::vector<double> ys) {
    assert(xs.size() == ys.size());

    // spline is of the form y = a + bx + cx^2 + dx^3
    const size_t N = ys.size();
    std::vector<double> a_tmp(N + 1), b_tmp(N + 1), c_tmp(N + 1, 0), d_tmp(N + 1);  // coefficients of cubic spline
    std::vector<double> h(N), C(N);  // C as in curvature
    std::vector<double> d(N), m(N), u(N);  // tridiagonal

    for(size_t i = 0; i < N; i++) {
        a_tmp[i + 1] = ys[i];
        h[i] = xs[i + 1] - xs[i];
    }

    for(size_t i = 1; i < N; i++) {
        // set A as in AX=b
        u[i] = h[i];
        m[i] = 2 * (h[i - 1] + h[i]);
        d[i] = h[i - 1];

        // set b as in AX=b
        C[i] = 3 * ((a_tmp[i + 2] - a_tmp[i + 1]) / h[i] - (a_tmp[i + 1] - a_tmp[i]) / h[i - 1]);
    }
    u.front() = 0;
    m.front() = m.back() = 1;
    d.back() = 0;
    C.front() = C.back() = 0;

    std::vector<double> tmp = PricerHelper::SolveTridiagonal(d, m, u, C);
    std::move(tmp.begin(), tmp.end(), c_tmp.begin() + 1);
    for(size_t i = 1; i <= N; i++) {
        b_tmp[i] = (a_tmp[i + 1] - a_tmp[i]) / h[i - 1] - h[i - 1] * (c_tmp[i + 1] + 2 * c_tmp[i]) / 3;
        d_tmp[i] = (c_tmp[i + 1] - c_tmp[i]) / 3 / h[i - 1];
    }

    // add left boundary spline
    a_tmp[0] = ys[0];
    b_tmp[0] = b_tmp[1];
    c_tmp[0] = d_tmp[0] = 0;

    // add right boundary spline
    double dx = xs.back() - xs.rbegin()[1];
    a_tmp.back() = ys.back();
    b_tmp.back() = b_tmp.rbegin()[1] + 2 * c_tmp.rbegin()[1] * dx + 3 * d_tmp.rbegin()[1] * dx * dx;
    c_tmp.back() = d_tmp.back() = 0;

    return Spline(std::move(xs), {std::move(a_tmp), std::move(b_tmp), std::move(c_tmp), std::move(d_tmp)});
}