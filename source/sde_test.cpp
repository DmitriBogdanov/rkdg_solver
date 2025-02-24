#include "firstparty/UTL.hpp"
#include "firstparty/gsse.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>

#include "sde_test/solver_wrappers.hpp"



// =====================
// --- Test problems ---
// =====================

namespace ode_test_problem {

using gsse::Vector;

constexpr double U = 5.;

inline const auto X0 = Vector{
    {10., 10.}
};

inline Vector A([[maybe_unused]] double t, [[maybe_unused]] const Vector& X) {
    return Vector{
        {U - X(0), X(0)}
    };
}

inline Vector B([[maybe_unused]] double t, [[maybe_unused]] const Vector& X) {
    return Vector{
        {0, 0}
    };
}

} // namespace ode_test_problem

namespace sde_test_problem {

using gsse::Vector;
using utl::math::sqr;

constexpr double alpha = -1.; // -1
constexpr double beta  = 1.;  // 1

inline const auto X0 = Vector{{0.}}; // 0

inline Vector A([[maybe_unused]] double t, [[maybe_unused]] const Vector& X) {
    return Vector{{-(alpha + sqr(beta) * X(0)) * (1. - sqr(X(0)))}};
}

inline Vector B([[maybe_unused]] double t, [[maybe_unused]] const Vector& X) {
    return Vector{{beta * (1 - sqr(X(0)))}};
}

inline Vector solution(double t, double Wt) {
    const double x0          = X0(0);
    const double nominator   = (1. + x0) * std::exp(-2. * alpha * t + 2. * beta * Wt) + x0 - 1.;
    const double denominator = (1. + x0) * std::exp(-2. * alpha * t + 2. * beta * Wt) - x0 + 1.;
    return Vector{{nominator / denominator}};
}

} // namespace sde_test_problem



int main() {
    using namespace utl;
    using namespace sde_test_problem;

    // --- Declare test problem parameters ---
    // ---------------------------------------

    const Vector     y0            = sde_test_problem::X0;
    const double     t1            = 0;
    const double     t2            = 1;
    const double     tau           = 1e-2;
    const gsse::uint callback_freq = (t2 - t1) / tau / 100; // 10 callbacks
    const gsse::uint seed          = 17;

    constexpr auto filepath_process           = "temp/sde_test_params/process.json";
    constexpr auto filepath_params            = "temp/sde_test_params/params.json";
    constexpr auto filepath_euler_maruyama    = "temp/sde_test_solution/euler_maruyama.json";
    constexpr auto filepath_modified_milstein = "temp/sde_test_solution/modified_milstein.json";

    json::Node solution_json = json::Array{};

    const auto callback = [&](const gsse::sde::State& state) {
        json::Node time_layer;
        time_layer["iteration"] = state.iteration;
        time_layer["t"]         = state.t;
        time_layer["y"]         = state.y0;
        solution_json.get_array().push_back(std::move(time_layer));
    };

    // Export params
    json::Node params_json;
    params_json["y1"]            = y0;
    params_json["t1"]            = t1;
    params_json["t2"]            = t2;
    params_json["tau"]           = tau;
    params_json["callback_freq"] = callback_freq;
    params_json["seed"]          = seed;
    std::filesystem::create_directories(std::filesystem::path(filepath_params).parent_path());
    params_json.to_file(filepath_params);

    // --- Export underlying Wiener's process & analythical solution ---
    // -----------------------------------------------------------------

    log::println("Exporting Wiener's process & analythical solution...");
    timer::start();

    const auto process_callback = [&](const gsse::sde::State& state) {
        json::Node time_layer;
        time_layer["iteration"] = state.iteration;
        time_layer["t"]         = state.t;
        time_layer["W"]         = state.y0;
        time_layer["sol"]       = sde_test_problem::solution(state.t, state.y0[0]);
        solution_json.get_array().push_back(std::move(time_layer));
    };

    const Vector W0                         = Vector{{0.}}; // Wiener's process always starts at zero
    const auto   trivial_deterministic_part = [](double, const Vector&) { return Vector{{0}}; };
    const auto   trivial_stochastic_part    = [](double, const Vector&) { return Vector{{1}}; };

    solution_json = json::Array{};
    gsse::sde::sde_solve(gsse::sde::Method::EULER_MARUYAMA, trivial_deterministic_part, trivial_stochastic_part, W0, t1,
                         t2, tau, process_callback, callback_freq, seed);
    solution_json.to_file(filepath_process);

    log::println(timer::elapsed_string_fullform(), '\n');

    // --- Solve test SDE using different methods ---
    // ---------------------------------------------

    // - Euler-Maruyama -
    log::println("Solving test SDE using Euler-Maruyama's method...");
    timer::start();

    const auto test_deterministic_part = sde_test_problem::A;
    const auto test_stochastic_part    = sde_test_problem::B;

    solution_json = json::Array{};
    gsse::sde::sde_solve(gsse::sde::Method::EULER_MARUYAMA, test_deterministic_part, test_stochastic_part, y0, t1, t2,
                         tau, callback, callback_freq, seed);
    solution_json.to_file(filepath_euler_maruyama);

    log::println(timer::elapsed_string_fullform(), '\n');

    // - Modified Milstein -
    log::println("Solving test SDE using Modified Miltein's method...");
    timer::start();

    solution_json = json::Array{};
    gsse::sde::sde_solve(gsse::sde::Method::MODIFIED_MILSTEIN, test_deterministic_part, test_stochastic_part, y0, t1,
                         t2, tau, callback, callback_freq, seed);
    solution_json.to_file(filepath_modified_milstein);

    log::println(timer::elapsed_string_fullform(), '\n');
}