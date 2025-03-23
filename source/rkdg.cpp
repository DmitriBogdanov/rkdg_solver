#include "rkdg/rkdg_convdiff.hpp"
#include "rkdg/rkdg_generic.hpp"
#include "rkdg/rkdg_pdf.hpp"
#include <cmath>
#include <cstddef>
#include <math.h>

// ===================
// --- PDF Problem ---
// ===================

namespace pdf_problem {
constexpr double C_E  = 1.;
constexpr double E    = 1.;
constexpr double wq_n = 2.;

constexpr double tau_Theta           = 4.;
constexpr double T_Theta             = 4.;
constexpr double avg_theta_g_squared = 2.;

constexpr auto W_Q = [](double Theta) -> double {
    const double numenator   = std::exp(Theta / (1. + Theta / E));
    const double demoninator = 1. + std::pow((C_E / E) * (1. + Theta / E), wq_n);
    return numenator / demoninator;
};

constexpr auto tau_Q = [](double Theta) -> double {
    // NOTE: This is likely incorrect
    constexpr double eps            = 0.5e-8;
    const double     W_Q_derivative = (W_Q(Theta + eps) - W_Q(Theta - eps)) / (eps * 2.);
    return (1. / W_Q_derivative);
};

constexpr auto f_Theta = [](double Theta, double t) -> double {
    // How do we compute 'tau_Q(Theta)' ?
    const double term_1 = 1. + (tau_Theta / T_Theta) - (tau_Theta / tau_Q(Theta));
    const double term_2 = 1. - std::exp(-t / tau_Theta * (term_1));
    return (1. / term_1) * term_2;
};

constexpr auto U = [](double Theta, [[maybe_unused]] double t) -> double { return -Theta / tau_Theta + W_Q(Theta); };

constexpr auto D = [](double Theta, double t) -> double {
    return 1. / tau_Theta * f_Theta(Theta, t) * avg_theta_g_squared;
};
} // namespace pdf_problem

double gaussian_distribution(double mean, double variance, double x) {
    using namespace utl;
    const double two_sigma_squared = 2. * math::sqr(variance);
    const double coef              = 1. / std::sqrt(math::constants::pi * two_sigma_squared);
    const double term              = std::exp(-math::sqr(x - mean) / two_sigma_squared);
    return coef * term;
    // computes N(mean, variance) PDF at 'x'
}

void solve_pdf() {

    constexpr std::size_t time_steps = 100;

    rkdg::IntegratorParams params;

    params.x_L              = -10;
    params.x_R              = 10;
    params.t1               = 0;
    params.t2               = 1;
    params.N                = 40;
    params.tau              = (t1 - t2) / time_steps;
    params.polynomial_order = 1;
    params.callback_freq    = time_steps / 100;

    constexpr double mean_temperature = 10; // initial normal distribution
    constexpr double variance         = 0.5;

    rkdg::CDProblem problem;
    problem.c    = pdf_problem::U;
    problem.a    = pdf_problem::D;
    problem.mu_0 = [](double x) { return gaussian_distribution(mean_temperature, variance, x); };
    problem.nu_L = [](double) { return 0.; };
    problem.nu_R = [](double) { return 0.; };
    
    rkdg::solve_convection_diffusion(params, problem, "temp/rkdg_solution/");
}

int main() {
    using namespace utl;
    
    log::println("Start: ", timer::datetime_string());
    
    // Solve
    timer::start();
    try {
        solve_pdf();
    } catch (const std::exception& except) { std::cerr << "ERROR: Caught exception!\n" << except.what() << "\n"; }
    const auto elapsed = timer::elapsed_string_fullform();
    
    log::println();
    log::println("Done!");
    log::println("End:    ", timer::datetime_string());
    log::println("Elapsed time: ", elapsed);
}