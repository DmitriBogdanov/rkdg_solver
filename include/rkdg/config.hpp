#pragma once

#include "firstparty/UTL.hpp"
#include <filesystem>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>

// Integration params
const std::string solution_dir = "temp/rkdg_solution/";

// constexpr double t1            = 0;                 // 't' interval start
// constexpr double t2            = 2.0;               // 't' interval end
// constexpr double X_L           = 0;                 // 'x' interval start
// constexpr double X_R           = utl::math::PI_TWO; // 'x' interval end
// constexpr size_t N             = 80;                // spacial grid N (j=0,N => number of grid points = N + 1)
// constexpr size_t time_steps    = (t2 - t1) / 5e-4;  // time steps
// constexpr size_t M             = 1;                 // basis polynomial order
// constexpr size_t callback_freq = time_steps / 100;  // how often are time layers saved

constexpr double t1            = 0;                // 't' interval start
constexpr double t2            = 2.0;              // 't' interval end
constexpr double X_L           = -20;              // 'x' interval start
constexpr double X_R           = 45;               // 'x' interval end
constexpr size_t N             = 40;               // spacial grid N (j=0,N => number of grid points = N + 1)
constexpr size_t time_steps    = (t2 - t1) / 5e-6; // time steps
constexpr size_t M             = 1;                // basis polynomial order
constexpr size_t callback_freq = time_steps / 100; // how often are time layers saved


// TEMP:
// Basis function
inline double phi(unsigned int l, double x) { return std::legendre(l, x); }

template <class T, std::enable_if_t<std::is_integral_v<T>, bool> = true>
std::size_t digits_in_integer(T value, bool count_minus = true) {
    std::size_t digits = 0;
    if constexpr (std::is_signed_v<T>)
        if (count_minus && value < 0) digits = 1;
    while (value) {
        value /= 10;
        digits++;
    }
    return digits;
}

inline void ensure_clear_directory(std::filesystem::path path) {
    std::filesystem::remove_all(path);
    std::filesystem::create_directory(path);
}