//
// Created by Albert Wong on 3/22/18.
//

#ifndef HULLWHITE_DIFFUSIONPDE_H
#define HULLWHITE_DIFFUSIONPDE_H

#include <vector>
#include <tuple>

class DiffusionPDE {
public:
    explicit DiffusionPDE() = default;
    virtual ~DiffusionPDE() = default;

    virtual inline double GetDrift(const double &x, const double &t) const = 0;
    virtual inline double GetDiffusion(const double &x, const double &t) const = 0;
    virtual inline double GetInterest(const double &x, const double &t) const = 0;

    std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> GetExplicitParameters(
            const std::vector<double> &xRanges,
            const double &t1,
            const double &t2
    ) const;

    std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> GetImplicitParameters(
            const std::vector<double> &xRanges,
            const double &t1,
            const double &t2
    ) const;
};


class BlackScholesPDE : public DiffusionPDE {
private:
    double m_r, m_sig2;

public:
    explicit BlackScholesPDE(const double &r, const double &sig) : m_r(r), m_sig2(sig * sig) {};
    ~BlackScholesPDE() override = default;

private:
    inline double GetDrift(const double &X, const double &t) const override { return m_r * X; };
    inline double GetDiffusion(const double &X, const double &t) const override { return m_sig2 * X * X; };
    inline double GetInterest(const double &X, const double &t) const override { return m_r; };
};


class HeatPDE : public DiffusionPDE {
public:
    explicit HeatPDE() = default;
    ~HeatPDE() override = default;

private:
    inline double GetDrift(const double &x, const double &t) const override { return 0.0; };
    inline double GetDiffusion(const double &x, const double &t) const override { return 1.0; };
    inline double GetInterest(const double &x, const double &t) const override { return 0.0; };
};


#endif //HULLWHITE_DIFFUSIONPDE_H
