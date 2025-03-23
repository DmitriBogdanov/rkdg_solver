#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <execution>
#include <functional>
#include <future>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <type_traits>
#include <vector>



// Polynomial P_N(x) = c0 + c1 * x + ... + cN * x^N
// defined by an array: coefs = [ c0 c1 ... cN ]
template <typename T, std::enable_if_t<std::is_arithmetic_v<T>, bool> = true>
struct Polynomial {
private:
    std::vector<double> coefs;

    void _check_validity() const { assert(!this->coefs.empty()); }

public:
    Polynomial() : coefs(1) {} // Trivial polynomial ~ a single const that equals to 0

    // - Constructors -
    explicit Polynomial(const std::vector<double>& init_coefs) : coefs(init_coefs) {
        if (init_coefs.empty()) throw std::invalid_argument("Attempted polynomial construction from an empty array.");
    }

    // - Getters -
    std::size_t order() const {
        this->_check_validity();

        return this->coefs.size() - 1;
    }

    bool is_trivial() const {
        this->_check_validity();

        return this->order() == 1 && this->coefs.front() == T();
    }

    // Computation
    T operator()(T x) const {
        this->_check_validity();

        double sum     = coefs.front();
        double x_power = 1.;

        for (std::size_t k = 1; k < this->coefs.size(); ++k) {
            x_power *= x;
            sum += this->coefs[k] * x_power;
        }

        return sum;
    }

    Polynomial derivative() const {
        this->_check_validity();

        if (this->order() == 0) return Polynomial();

        std::vector<double> diff_coefs(this->coefs.size() - 1);
        for (std::size_t k = 0; k < diff_coefs.size(); ++k) diff_coefs[k] = this->coefs[k + 1] * (k + 1.);

        return Polynomial(diff_coefs);
    }

    std::string stringify() const {
        this->_check_validity();

        std::stringstream ss;

        ss << "P_" << this->order() << "(x) = " << this->coefs.front();

        for (std::size_t k = 1; k < this->coefs.size(); ++k) {
            if (this->coefs[k] == T()) continue;
            ss << coefs[k] << " * x^" << k << (k == this->coefs.size() - 1 ? "" : "  +  ");
        }

        return ss.str();
    }
};

using Poly = Polynomial<double>;

inline std::vector<Poly> diff_polynomials(const std::vector<Poly>& polys) {
    std::vector<Poly> res;
    res.reserve(polys.size());

    for (const auto& poly : polys) res.push_back(poly.derivative());

    return res;
}


// Manually create table of legendre polynomials
// This is a dumb way of getting 'pre-evaluated' legendre polynomials without having to compute their
// coefs in general case at every call. Saves us lots of performance and we don't really need large orders
// and derivative in practice.

inline const std::vector<Poly> legendre_diff0 = {
    Poly(std::vector<double>({1.})),                                  // M = 0
    Poly(std::vector<double>({0., 1.})),                              // M = 1
    Poly(std::vector<double>({-1. / 2., 0., 3. / 2.})),               // M = 2
    Poly(std::vector<double>({0., -3. / 2., 0., 5. / 2.})),           // M = 3
    Poly(std::vector<double>({3. / 8., 0., -30. / 8., 0., 35. / 8.})) // M = 4
};

inline const std::vector<Poly> legendre_diff1 = diff_polynomials(legendre_diff0);
inline const std::vector<Poly> legendre_diff2 = diff_polynomials(legendre_diff1);
inline const std::vector<Poly> legendre_diff3 = diff_polynomials(legendre_diff2);
inline const std::vector<Poly> legendre_diff4 = diff_polynomials(legendre_diff3);

inline const std::vector<std::vector<Poly>> legendre_diffs = {legendre_diff0, legendre_diff1, legendre_diff2,
                                                              legendre_diff3, legendre_diff4};

inline double evaluate_legendre_derivative(std::size_t l, std::size_t derivative, double x) {
    assert(l < legendre_diff0.size());          // Max supported polynomial order
    assert(derivative < legendre_diffs.size()); // Max supported derivative

    return legendre_diffs[derivative][l](x);
}



// Just testing out std::thread parallelization
template <typename IdxType = std::size_t, typename FuncType,
          typename std::enable_if_t<std::is_convertible_v<FuncType, std::function<void(IdxType)>>, bool> = true>
inline void parallel_for(IdxType start, IdxType end, FuncType func, std::size_t thread_count_override = 0) {
    // Execute sequentially if single thread was manually specified
    if (thread_count_override == 1) {
        for (IdxType j = start; j < end; ++j) func(j);
        return;
    }

    const IdxType workload = end - start;

    // Get a hint with a number of threads supported by hardware
    // (hint == 0) => Couldn't get info, so we assume a default number (8)
    constexpr IdxType thread_count_fallback = 4;
    const IdxType     thread_count_hint     = std::thread::hardware_concurrency();
    const IdxType     thread_count          = (thread_count_override) ? thread_count_override
                                              : (thread_count_hint)   ? thread_count_hint
                                                                      : thread_count_fallback;

    // Split work into batches
    const IdxType batch_size      = workload / thread_count;
    const IdxType batch_remainder = workload % thread_count;

    const auto execute_batch = [&](IdxType batch_start, IdxType batch_end) {
        for (IdxType j = batch_start; j < batch_end; ++j) func(j);
    };

    // Execute on multiple threads
    std::vector<std::thread> thread_pool(thread_count);

    for (std::size_t i = 0; i < thread_count; ++i) {
        IdxType start  = i * batch_size;
        thread_pool[i] = std::thread(execute_batch, start, start + batch_size);
    }

    // Handle remainder work
    const IdxType leftover_work_start = thread_count * batch_size;
    execute_batch(leftover_work_start, leftover_work_start + batch_remainder);

    // Wait for the other thread to finish their task
    std::for_each(thread_pool.begin(), thread_pool.end(), std::mem_fn(&std::thread::join));
}