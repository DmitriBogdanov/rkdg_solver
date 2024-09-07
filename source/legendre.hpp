#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <functional>
#include <vector>


struct Polynomial {
    // coefs = [ c0 c1 ... cN ]
    // define polynomial P(x) = c0 + c1 * x + ... + cN * x^N
    std::vector<double> coefs;

    Polynomial() = default;
    Polynomial(const std::vector<double>& coefs) : coefs(coefs) {}
    Polynomial(std::vector<double>&& coefs) : coefs(std::move(coefs)) {}

    Polynomial diff() const {
        if (this->coefs.size() <= 1) return Polynomial();

        std::vector<double> diff_coefs(this->coefs.size() - 1);
        for (std::size_t k = 0; k < diff_coefs.size(); ++k) diff_coefs[k] = this->coefs[k + 1] * (k + 1.);

        return diff_coefs;
    }

    double evaluate_at(double x) const {
        double value = 0.;
        double x_pwr = 1.;

        for (std::size_t k = 0; k < this->coefs.size(); ++k) {
            value += this->coefs[k] * x_pwr;
            x_pwr *= x;
        }

        return value;
    }
};

inline std::vector<Polynomial> diff_polynomials(const std::vector<Polynomial>& polys) {
    std::vector<Polynomial> res;
    res.reserve(polys.size());

    for (const auto& poly : polys) res.push_back(poly.diff());

    return res;
}


// Manually create table of legendre polynomials
// This is a dumb way of getting 'pre-evaluated' legendre polynomials without having to compute their
// coefs in general case at every call. Saves us lots of performance and we don't really need large orders
// and derivative in practice.

inline const std::vector<Polynomial> legendre_diff0 = {
    Polynomial({1.}),                                  // M = 0
    Polynomial({0., 1.}),                              // M = 1
    Polynomial({-1. / 2., 0., 3. / 2.}),               // M = 2
    Polynomial({0., -3. / 2., 0., 5. / 2.}),           // M = 3
    Polynomial({3. / 8., 0., -30. / 8., 0., 35. / 8.}) // M = 4
};

inline const std::vector<Polynomial> legendre_diff1 = diff_polynomials(legendre_diff0);
inline const std::vector<Polynomial> legendre_diff2 = diff_polynomials(legendre_diff1);
inline const std::vector<Polynomial> legendre_diff3 = diff_polynomials(legendre_diff2);
inline const std::vector<Polynomial> legendre_diff4 = diff_polynomials(legendre_diff3);

inline const std::vector<std::reference_wrapper<const std::vector<Polynomial>>> legendre_diffs = {
    legendre_diff0, legendre_diff1, legendre_diff2, legendre_diff3, legendre_diff4};

inline double evaluate_legendre_derivative(std::size_t l, std::size_t derivative, double x) {
    assert(l < legendre_diff0.size());          // Max supported polynomial order
    assert(derivative < legendre_diffs.size()); // Max supported derivative

    return legendre_diffs[derivative].get()[l].evaluate_at(x);
}