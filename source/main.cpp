#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <ctime>
#include <exception>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "firstparty/gsse.hpp"
#include "firstparty/proto_utils.hpp"

#include "legendre.hpp"



// Integration params
const std::string solution_dir = "temp/solution/";

constexpr double t1            = 0;                 // 't' interval start
constexpr double t2            = 2.0;               // 't' interval end
constexpr double X_L           = 0;                 // 'x' interval start
constexpr double X_R           = utl::math::PI_TWO; // 'x' interval end
constexpr size_t N             = 50;                // spacial grid N (j=0,N => number of grid points = N + 1)
constexpr size_t time_steps    = (t2 - t1) / 5e-4;  // time steps
constexpr size_t M             = 2;                 // basis polynomial order
constexpr size_t callback_freq = time_steps / 100;  // how often are time layers saved

// Problem params
// u_t + c u_x - a u_xx = 0
constexpr double c    = 1;                                                            // coef 1
constexpr double a    = 1;                                                            // coef 2
constexpr auto   mu_0 = [](double x) { return std::sin(x); };                         // initial condition
constexpr auto   nu_L = [](double t) { return -std::exp(-a * t) * std::sin(c * t); }; // left BC
constexpr auto   nu_R = [](double t) { return nu_L(t); };                             // right BC

// Basis function
inline double phi(unsigned int l, double x) { return std::legendre(l, x); }

template <class IntegerType>
std::size_t digits_in_integer(IntegerType value) {
    std::size_t digits = 0;
    if (value < 0) digits = 1; // remove this line if '-' counts as a digit
    while (value) {
        value /= 10;
        digits++;
    }
    return digits;
}

// System builder
void rkdg_solve_convection_diffusion() {
    using namespace utl;

    const size_t N_j = N;
    const size_t N_m = M + 1;

    // Clear directory
    std::error_code err;
    std::filesystem::remove_all(solution_dir, err);
    if (err) {
        std::cerr << "ERROR: Error erasing file -> " << err.message() << std::endl;
        std::exit(1);
    }
    std::filesystem::create_directory(solution_dir);

    // Build spacial grid {x_j} (0 <= j <= N)
    const std::vector<double> x = math::linspace(X_L, X_R, math::Intervals(N));

    // Compute grid steps {h_j} (0 <= j <= N-1)
    std::vector<double> h(N);
    for (size_t j = 0; j <= N - 1; ++j) h[j] = x[j + 1] - x[j];

    // Compute coefs {R_nm} ( 0 <= n, m <= M)
    const std::function<double(int, int)> recursively_compute_R = [&recursively_compute_R](int n, int m) -> double {
        assert(-1 <= n && n <= static_cast<int>(M));
        assert(-1 <= m && m <= static_cast<int>(M));

        if (n == -1) return 0;
        if (n == 0) return 0;
        if (n < m + 1) return 0;

        return 2. * math::kronecker_delta(n - 1, m) + recursively_compute_R(n - 2, m); // #expression (3.2)

        // Recursive lambda for computing {R_nm}.
        //
        // std::function<...> here is necessary due to recursion.
        // 'int' arguments are necessary due to "fake" value at "n = -1" that is needed for recursion
        //
        // This is NOT a performant way of doing it (we compute same values of R_nm repeatedly), but since
        // realistically 'M' isn't going to be larger that ~10, it doesn't really matter compared to main computation
    };

    gsse::Matrix R(M + 1, M + 1);
    // mvl::Matrix<double> R(M + 1, M + 1);

    for (size_t n = 0; n <= M; ++n)
        for (size_t m = 0; m <= M; ++m) R(n, m) = recursively_compute_R(static_cast<int>(n), static_cast<int>(m));

    // NOTE:
    // Seems to solve issues with partial differentiation blowing up due to index "misalignment".
    // However, it's a bit unclear where I messed up the order of variables.
    gsse::Matrix R_transposed = R.transpose();
    R                         = R_transposed;

    // Compute {alpha_jm(0)}
    gsse::Matrix alpha0(N_j, N_m);

    for (size_t j = 0; j < N_j; ++j) {
        for (size_t m = 0; m < N_m; ++m) {
            const double front_coef          = (2. * m + 1) * 0.5;
            const auto   integrated_function = [&](double xi) -> double {
                return mu_0(x[j] + 0.5 * h[j] * (xi + 1.)) * phi(m, xi);
            };
            const double integral = math::integrate_trapezoidal(integrated_function, -1., 1., math::Intervals(800));
            alpha0(j, m)          = front_coef * integral; // #expression (6.2)
            // NOTE: In preprint this expression is computed analythically for mu_0(x) = sin(x),
            // however I couldn't derive formula for arbitrary M so numerical integration it is.
            // Also makes method more generic since we no longer rely on specific mu_0(x).
        }
    }

    gsse::Vector y0(alpha0.reshaped()); // ODE solver wants vector, not a matrix

    // Declare F(t, y)
    gsse::VectorTimeFunction F = [&](double t, const gsse::Vector& y) -> gsse::Vector {
        // Reshape 1D data into {alpha_jm(t)} matrix
        // (working with matrices directly allows us to avoid all 1D <-> 2D index converion issues)
        gsse::Matrix alpha = y.reshaped(N_j, N_m);

        // 1) Using {alpha_jm(t)} compute {u_j^-}, {u_j^+}
        gsse::Vector u_minus(N_j);
        gsse::Vector u_plus(N_j);

        u_minus[0] = nu_L(t); // #expression (6.4)

        for (size_t j = 1; j <= N - 1; ++j) {
            double sum_1 = 0;
            double sum_2 = 0;

            for (size_t n = 0; n <= M; ++n) sum_1 += alpha(j, n) * math::power_of_minus_one(n);
            for (size_t n = 0; n <= M; ++n) sum_2 += alpha(j - 1, n);

            u_minus[j] = 0.5 * (sum_1 + sum_2); // #expression (5.7)
        }

        u_plus[N - 1] = nu_R(t); // #expression (6.4)

        for (size_t j = 0; j <= N - 2; ++j) {
            u_plus[j] = u_minus[j + 1]; // #expression (5.7)
        }

        // 2) Knowing {u_j^-}, {u_j^+} we can compute {beta_jm(t)}
        gsse::Matrix beta(N_j, N_m);

        for (size_t j = 0; j < N_j; ++j) {
            for (size_t m = 0; m < N_m; ++m) {
                const double coef_1 = c;
                const double coef_2 = a * (2. * m + 1) / h[j];

                const double u_diff = u_plus[j] - u_minus[j] * math::power_of_minus_one(m);

                double sum = 0;
                for (size_t n = 0; n <= M; ++n) sum += alpha(j, n) * R(n, m);

                beta(j, m) = coef_1 * alpha(j, m) - coef_2 * (u_diff - sum); // #expression (5.6)
            }
        }

        // 3) Knowing {beta_jm(t)} we can compute {f_j^-}, {f_j^+}
        gsse::Vector f_minus(N_j);
        gsse::Vector f_plus(N_j);

        {
            double sum_1 = 0;
            for (size_t n = 0; n <= M; ++n) sum_1 += beta(0, n) * math::power_of_minus_one(n);

            f_minus[0] = sum_1; // #expression (6.3)
        }

        for (size_t j = 1; j <= N - 1; ++j) {
            double sum_1 = 0;
            double sum_2 = 0;

            for (size_t n = 0; n <= M; ++n) sum_1 += beta(j, n) * math::power_of_minus_one(n);
            for (size_t n = 0; n <= M; ++n) sum_2 += beta(j - 1, n);

            f_minus[j] = 0.5 * (sum_1 + sum_2); // #expression (4.4)
        }

        {
            double sum_2 = 0;
            for (size_t n = 0; n <= M; ++n) sum_2 += beta(N - 1, n);

            f_plus[N - 1] = sum_2; // #expression (6.3)
        }

        for (size_t j = 0; j <= N - 2; ++j) {
            f_plus[j] = f_minus[j + 1]; // #expression (4.4)
        }

        // 4) Knowing {f_j^-}, {f_j^+} we can compute RHS F(t, y)
        gsse::Matrix d_alpha_d_t(N_j, N_m);

        for (size_t j = 0; j < N_j; ++j) {
            for (size_t m = 0; m < N_m; ++m) {
                const double front_coef = -(2. * m + 1.) / h[j];

                const double f_diff = f_plus[j] - f_minus[j] * math::power_of_minus_one(m);

                double sum = 0;
                for (size_t n = 0; n <= M; ++n) sum += beta(j, n) * R(n, m);

                d_alpha_d_t(j, m) = front_coef * (f_diff - sum); // #expression (4.3)
            }
        }

        // > Start
        // UTL_LOG("\nValues of {alpha_jm(t)}:\n", alpha.format(gsse::FORMAT_DEFAULT), "\n\n");

        // > 1)
        // UTL_LOG("\nValues of {u_j^-}:\n", u_minus.format(gsse::FORMAT_INLINE), "\n\n");
        // UTL_LOG("\nValues of {u_j^+}:\n", u_plus.format(gsse::FORMAT_INLINE), "\n\n");

        // > 2)
        // UTL_LOG("\nValues of {beta_jm(t)}:\n", beta.format(gsse::FORMAT_DEFAULT), "\n\n");

        // > 3)
        // UTL_LOG("\nValues of {f_j^-}:\n", f_minus.format(gsse::FORMAT_INLINE), "\n\n");
        // UTL_LOG("\nValues of {f_j^+}:\n", f_plus.format(gsse::FORMAT_INLINE), "\n\n");

        // > 4)
        // UTL_LOG("\nValues of {d_alpha_d_t(t)}:\n", d_alpha_d_t.format(gsse::FORMAT_DEFAULT), "\n\n");

        // Reshape matrix of derivatives into a 1D vector
        gsse::Vector result(d_alpha_d_t.reshaped());

        return result;
    };

    // Callback for saving results
    progressbar::Percentage progressbar;

    const gsse::ode::CallbackFunction callback = [&](gsse::ode::State& state) -> void {
        // UTL_LOG("Callback: saving results");

        progressbar.set_progress(static_cast<double>(state.iteration) / time_steps);

        // Restore solution u(x) ~ {u_j} using decomposition coefs
        std::vector<double> u(x.size());

        gsse::Matrix alpha = state.y0.reshaped(N_j, N_m);

        u[0] = nu_L(state.t); // #expression (6.4)
        for (size_t j = 1; j <= N - 1; ++j) {
            double sum_1 = 0;
            double sum_2 = 0;

            for (size_t n = 0; n <= M; ++n) sum_1 += alpha(j, n) * math::power_of_minus_one(n);
            for (size_t n = 0; n <= M; ++n) sum_2 += alpha(j - 1, n);

            u[j] = 0.5 * (sum_1 + sum_2); // #expression (5.7)
        }
        u[N] = nu_R(state.t); // #expression (6.4)

        config::export_json(solution_dir + stre::pad_with_zeroes(state.iteration, digits_in_integer(time_steps)) +
                                ".json",
                            config::entry("t", state.t), config::entry("alpha", state.y0), config::entry("u", u));
    };

    // Forward to ODE solver
    progressbar.start();

    gsse::ode::ode_solve(gsse::ode::Method::RK4,       // method
                         F,                            // system RHS
                         y0,                           // system Initial Condition
                         t1, t2 + 1e-6,                // time interval, NOTE: +1e-6 ensures we don't miss last "t"
                         (t2 - t1) / (time_steps - 1), // time step
                         callback,                     // callback for saving results
                         callback_freq                 // how often to save results
    );

    progressbar.finish();

    // Export params
    config::export_json(solution_dir + "params.json",                            //
                        config::entry("date", timer::datetime_string().data()),  //
                        config::entry("problem", "Convection-diffusion test 1"), //
                        config::entry("problem_coef_c", c),                      //
                        config::entry("problem_coef_a", a),                      //
                        config::entry("t1", t1),                                 //
                        config::entry("t2", t2),                                 //
                        config::entry("X_L", X_L),                               //
                        config::entry("X_R", X_R),                               //
                        config::entry("x", x),                                   //
                        config::entry("N", N),                                   //
                        config::entry("time_steps", time_steps),                 //
                        config::entry("callback_freq", callback_freq));
}


int main() {
    using namespace utl;
    std::cout << "\n" << timer::datetime_string() << "\n";
    timer::start();

    // Solve
    std::ofstream log_file("temp/log.txt");
    UTL_LOG_SET_OUTPUT(voidstream::vout);

    try {
        rkdg_solve_convection_diffusion();
    } catch (const std::exception& except) { std::cerr << "ERROR: Caught exception!\n" << except.what() << "\n"; }

    // std::cout << std::ifstream("temp/log.txt").rdbuf();

    gsse::Matrix mat = gsse::make_matrix({
        {1, 2, 3},
        {4, 5, 6},
    });

    std::cout << "\nDone!\n"
              << "Elapsed time: " << timer::elapsed_string_fullform() << "\n"
              << timer::datetime_string() << "\n\n";
    
    return 0;
}