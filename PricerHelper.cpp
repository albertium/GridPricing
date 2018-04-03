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
    // todo: add assert and probably remove offset
    assert(pd.size() == pm.size());
    assert(pm.size() == pu.size());

    std::vector<double> result(pm.size() + 2 * offset);

    // forward
    std::vector<double> pu_p(pm.size()), d_p(pm.size());
    pu_p[0] = pu[0] / pm[0];
    d_p[0] = d[offset] / pm[0];
    for(size_t i = 1; i < pm.size(); ++i) {
        double den = pm[i] - pd[i] * pu_p[i - 1];
        pu_p[i] = pu[i] / den;
        d_p[i] = (d[offset + i] - pd[i] * d_p[i - 1]) / den;
    }

    // backward
    result.rbegin()[offset] = d_p.back();
    for(auto i = static_cast<int>(pm.size()) - 2; i >= 0; --i) {
        result[i + offset] = d_p[i] - pu_p[i] * result[i + 1 + offset];
    }

    return result;
}


std::vector<double> PricerHelper::SolveTridiagonal(const std::vector<std::vector<double>> A,
                                                   const std::vector<double> d) {
    /*
     * Thomas algorithm
     * A[0][0] and A[N][N] are not used
     */
    assert(A.size() == d.size());

    const size_t N = d.size();
    std::vector<double> cTmp(N), dTmp(N);

    // walk forward
    cTmp[0] = A[0][2] / A[0][1];
    dTmp[0] = d[0] / A[0][1];
    for (size_t i = 1; i < N; i++) {
        double den = A[i][1] - A[i][0] * cTmp[i - 1];
        cTmp[i] = A[i][2] / den;
        dTmp[i] = (d[i] - A[i][0] * dTmp[i - 1]) / den;
    }

    // walk backward
    for (int i = static_cast<int>(N - 2); i >= 0; i--) {
        dTmp[i] = dTmp[i] - cTmp[i] * dTmp[i + 1];
    }

    return dTmp;
}

std::vector<double> PricerHelper::SolveTridiagonalAxBd(std::vector<std::vector<double>> A,
                                                       std::vector<std::vector<double>> B,
                                                       std::vector<double> d,
                                                       std::vector<double> xBounds,
                                                       std::vector<double> dBounds) {
    /*
     * each of A or B should be of size 3
     * If xBounds is provided, then we assume a(0, 0) * xBound + a(0, 1) * x(0, 0) + a(0, 2) * x(0, 1)
     * If xBounds is not provided, then assume a(0, 0) * x(0, 0) + a(0, 1) * x(0, 1) + a(0, 2) * x(0, 2)
     * This is the same for B and d
     *
     * To solve for regular Ax=Bd, we need to set both xBounds and dBounds to {0, 0}
     */

    assert(d.size() > 2);
    assert(A.size() == d.size() || A.empty());
    assert(B.size() == d.size() || B.empty());
    assert(!A.empty() || !B.empty());
    assert(xBounds.empty() || xBounds.size() == 2);
    assert(dBounds.empty() || dBounds.size() == 2);

    std::vector<double> result(d.size());
    const size_t N = d.size();

    // solve for D = Bd first
    if (!B.empty()) {
        if (dBounds.empty()) {
            result[0] = B[0][0] * d[0] + B[0][1] * d[1] + B[0][2] * d[2];
            result.back() = B.back()[0] * d.rbegin()[2] + B.back()[1] * d.rbegin()[1] + B.back()[2] * d.back();
        } else {
            result[0] = B[0][0] * dBounds[0] + B[0][1] * d[0] + B[0][2] * d[1];
            result.back() = B.back()[0] * d.rbegin()[1] + B.back()[1] * d.back() + B.back()[2] * dBounds[1];
        }

        for (size_t i = 1; i < N - 1; i++)
            result[i] = B[i][0] * d[i - 1] + B[i][1] * d[i] + B[i][2] * d[i + 1];
    } else
        result.swap(d);

    // solve for Ax = D now
    if (!A.empty()) {
        if (xBounds.empty()) {
            double ratio = A[0][2] / A[1][2];
            A[0][2] = A[0][1] - A[1][1] * ratio;
            A[0][1] = A[0][0] - A[1][0] * ratio;
            result[0] -= result[1] * ratio;

            ratio = A[N - 1][N - 3] / A[N - 2][N - 3];
            A.back().rbegin()[2] = A.back().rbegin()[1] - A.rbegin()[1].rbegin()[1] * ratio;
            A.back().rbegin()[1] = A.back().back() - A.rbegin()[1].back() * ratio;
            result[N - 1] -= result[N - 2] * ratio;
        } else {
            result[0] -= A[0][0] * xBounds[0];
            result.back() -= A.back().back() * xBounds[1];
        }

        result = PricerHelper::SolveTridiagonal(A, result);
    }

    return result;
}


std::vector<double> PricerHelper::SolveSOR(const std::vector<double> &pd, const std::vector<double> &pm,
                                           const std::vector<double> &pu, const std::vector<double> &d,
                                           std::vector<double> guess, const std::function<double(double)> &filter,
                                           const double &omega, const double &tol) {
    assert(pd.size() == pm.size());
    assert(pm.size() == pu.size());
    assert(pm.size() + 2 == guess.size());

    double err;
    size_t count = 0;

    do {
        err = 0.0;
        if(filter) {
            for(size_t i = 0; i < pm.size(); ++i)
                guess[i + 1] = filter((1 - omega) * guess[i + 1]
                                       + omega * (d[i] - pd[i] * guess[i] - pu[i] * guess[i + 2]) / pm[i]);
        } else {
            for(size_t i = 0; i < pm.size(); ++i) {
                double r = omega / pm[i] * (d[i] - pd[i] * guess[i] - pm[i] * guess[i + 1] - pu[i] * guess[i + 2]);
                guess[i + 1] = guess[i + 1] + r;
                err += abs(r);
            }
        }
        count++;
    } while (err > tol * tol && count < 100);

    return guess;
}

double PricerHelper::GetBlackScholesPrice(double S, double K, double r, double q, double sig, double T) {
    double sig2 = sig * sig;
    double sqrt_T = sqrt(T);
    double d1 = (log(S / K) + (r - q + sig2 / 2) * T) / sig / sqrt_T;
    double d2 = d1 - sig * sqrt_T;
    return S * exp(-q*T) * NormalCDF(d1) - K * exp(-r*T) * NormalCDF(d2);
}