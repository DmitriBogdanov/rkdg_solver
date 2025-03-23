// TODO: Beautify this header some more and turn it into a full-on project

#pragma once

#include <cmath>
#include <functional>
#include <initializer_list>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>

#include "thirdparty/Eigen/Dense"


// GSSE = "Generic System Solvers for Eigen"
namespace gsse {


// --- Typedefs ---
// ----------------
using uint = unsigned int;
using sint = int;
using idx  = Eigen::Index; // ptrdiff_t

using Scalar = double; // switching these typedefs float equivalents will this whole header to operate on floats
using Vector = Eigen::VectorXd;
using Matrix = Eigen::MatrixXd;

using ScalarFunction     = std::function<Scalar(Scalar)>;
using ScalarTimeFunction = std::function<Scalar(Scalar, Scalar)>;
using VectorFunction     = std::function<Vector(const Vector&)>;
using VectorTimeFunction = std::function<Vector(Scalar, const Vector&)>;


// --- Formats ---
// ---------------
inline const Eigen::IOFormat FORMAT_DEFAULT(6, 0, ", ", "\n", "[", "]");
inline const Eigen::IOFormat FORMAT_INLINE(6, 0, ", ", ", ", "", "", "{", "}");
inline const Eigen::IOFormat FORMAT_NONE(6, 0, " ", " ", "", "", "", "");


// --- Matrix initializers ---
// ---------------------------
inline Matrix make_matrix(std::initializer_list<std::initializer_list<Scalar>> matrix) {
    // Deduce size
    idx rows = matrix.size();
    idx cols = (*matrix.begin()).size();

    // Check that column sizes aren't messed up
    for (const auto& row : matrix)
        if (cols != static_cast<idx>(row.size())) std::cerr << "WARNING: make_matrix(): Incompatible column sizes.\n";

    // Random access indexing for std::initializer_list<> 'matrix'
    auto take_value_at = [&](idx i, idx j) -> Scalar { return (matrix.begin()[i]).begin()[j]; };

    // Init and fill with values
    Matrix res = Matrix::Zero(rows, cols);

    for (idx i = 0; i < rows; ++i)
        for (idx j = 0; j < cols; ++j) res(i, j) = take_value_at(i, j);

    return res;
}

inline Vector make_vector(std::initializer_list<Scalar> vector) {
    // Deduce size
    idx rows = vector.size();

    // Random access indexing for std::initializer_list<> 'vector'
    auto take_value_at = [&](idx i) -> Scalar { return vector.begin()[i]; };

    // Init and fill with values
    Vector res = Vector::Zero(rows);

    for (idx i = 0; i < rows; ++i) res(i) = take_value_at(i);

    return res;
}


// --- Converters ---
// ------------------

// Useful for when we need to feed matrix to an API that expects 'array of arrays'
inline std::vector<std::vector<Scalar>> matrix_to_nested_std_vector(const Matrix& matrix) {
    // Resize std::vector<std::vector<>> to fit a matrix
    std::vector<std::vector<Scalar>> res(matrix.rows());
    for (auto& res_row : res) res_row.resize(matrix.cols());

    // Copy matrix
    for (idx i = 0; i < matrix.rows(); ++i)
        for (idx j = 0; j < matrix.cols(); ++j) res[i][j] = matrix(i, j);

    return res;
}

// Usefult when we need to feed Eigen matrix/vector to an API that expects a more 'standard' array
inline std::vector<Scalar> matrix_to_std_vector(const Matrix& matrix) {
    // Since matrix doesn't come with 1D views, we use access to underlying 1D array
    auto raw_data = matrix.data();
    auto raw_size = matrix.size();

    // Copy array to vector
    std::vector<Scalar> res(raw_size);
    for (idx i = 0; i < raw_size; ++i) res[i] = raw_data[i];


    return res;
}


// --- Algebraic solvers ---
// -------------------------
namespace algebraic {

// Computes Jacobian of F(y) at y = X using numeric differentiation
inline Matrix jacobian(VectorFunction F, Vector X) {
    constexpr Scalar DIFF_EPS         = 2e-8;
    constexpr Scalar DIFF_INVERSE_EPS = 1. / DIFF_EPS;

    const auto N = X.size();

    Matrix J = Matrix::Zero(N, N);

    Vector dXj = Vector::Zero(N);

    for (int j = 0; j < N; ++j) {
        // Differentiate F with respect to x_j
        dXj(j) += DIFF_EPS;
        const auto dFj = (F(X + 0.5 * dXj) - F(X - 0.5 * dXj)) * DIFF_INVERSE_EPS;
        dXj(j)         = 0.;

        // As a result we get an entire column of jacobian
        J.col(j) = dFj;
    }

    return J;
}

// Solve system of algebraic equations F(y) = 0 using y0 as an initial guess
inline Vector nonlinear_solve(VectorFunction F, Vector y0, double precision, uint maxIterations = 100) {
    Vector y = y0;

    // Fill approximations
    for (uint iterations = 0; iterations < maxIterations; ++iterations) {
        y = y0 - jacobian(F, y0).inverse() * F(y0); // Newthon's method iteration

        if ((y - y0).norm() < precision) break;

        y0 = y;
    }

    return y;
}

} // namespace algebraic


// --- SDE solvers ---
// -------------------
namespace sde {

class NormalDistrbutionGenerator {
private:
    std::random_device         random_device;
    std::mt19937               random_generator;
    std::normal_distribution<> random_distribution;

public:
    inline NormalDistrbutionGenerator() : random_device(), random_generator(random_device()), random_distribution() {}

    inline void   seed(unsigned int seed) { this->random_generator.seed(seed); }
    inline Scalar get_next() { return random_distribution(random_generator); }
};

enum class Method {
    EULER_MARUYAMA,   // explicit, weak order = 1, strong order = 1/2
    MODIFIED_MILSTEIN // explicit, weak order = 1, strong order = 1
};

struct State {
    const Scalar& t;          // current time
    const uint&   iteration;  // current iteration
    const Vector& y0;         // current x
    const Vector& y_computed; // 'next' x
};

using CallbackFunction = std::function<void(State&)>;

// - Solver implementations -

// Simplest method, equivalent of Euler's method for ODEs but for SDEs
inline void _step_euler_maruyama(VectorTimeFunction& A, VectorTimeFunction& B, Scalar t, Scalar tau,
                                 NormalDistrbutionGenerator& generator, Vector& y0, Vector& y) {
    const auto sqrt_of_tau = std::sqrt(tau);

    Vector dW(y0.rows());
    for (auto& elem : dW) elem = generator.get_next() * sqrt_of_tau;

    y = y0 + tau * A(t, y0) + (dW.array() * B(t, y0).array()).matrix();
}

// Method from "[P.E. Kloeden, E. Platen] Numerical solution of stochastic differential equations" (page 374)
// Unnamed modification of the Milstein method,
// written for Ito process A = A(X), B = B(X)
// Has a strong order of 1 and a weak order of 1
// For ODE's equvalent to explicit Euler method
inline void _step_modified_milstein(VectorTimeFunction& A, VectorTimeFunction& B, Scalar t, Scalar tau,
                                    NormalDistrbutionGenerator& generator, Vector& y0, Vector& y) {
    const auto sqrt_of_tau = std::sqrt(tau);

    Vector dW(y0.rows());
    for (auto& elem : dW) elem = generator.get_next() * sqrt_of_tau;

    Vector tau_column(y0.rows());
    for (auto& elem : tau_column) elem = tau;

    const Vector Y = y0 + tau * A(t, y0) + sqrt_of_tau * B(t, y0);

    y = y0 + tau * A(t, y0) + (dW.array() * B(t, y0).array()).matrix() +
        (0.5 * (B(t, Y) - B(t, y0)).array() * (dW.array() * dW.array() - tau_column.array()) / sqrt_of_tau).matrix();
}

// Solve stochastic diff. equation (SDE) with Ito form:
// dX(t) = A(t, X(t)) dt +  B(t, X(t)) dW
// with dW = eta(t) * sqrt(t) and eta(t) denotes white noise
inline void sde_solve(Method             method,           // solver
                      VectorTimeFunction A,                // determenistic part
                      VectorTimeFunction B,                // stochastic part
                      Vector             y0,               // initial condition
                      Scalar t1, Scalar t2,                // time interval
                      Scalar           tau,                // time step
                      CallbackFunction callback,           // callback function
                      uint             callback_frequency, // callback frequency (in iterations)
                      uint             random_seed         // random generator seed
) {
    // 'x0' passed by value cause we use it as an internal variable

    NormalDistrbutionGenerator generator;
    generator.seed(random_seed);

    // Iteration
    Vector y         = y0;
    Scalar t         = t1;
    uint   iteration = 0;

    State state{t, iteration, y0, y};

    for (t = t1; t < t2; t += tau, ++iteration) {
        // Get 'x' from 'x0' using appropriate solver
        y0 = y;

        switch (method) {
        case Method::EULER_MARUYAMA: _step_euler_maruyama(A, B, t, tau, generator, y0, y); break;
        case Method::MODIFIED_MILSTEIN: _step_modified_milstein(A, B, t, tau, generator, y0, y); break;
        default: return;
        }

        // Check for solution divergence
        for (idx k = 0; k < y.size(); ++k) {
            // Stop and throw if any values in 'y' are 'Inf' or 'NaN'
            if (std::isfinite(y[k])) continue;
            std::stringstream ss;
            ss << "Inf or NaN encountered during SDE solution at { iteration = " << iteration << ", t = " << t
               << " } on index { k = " << k << " / " << y.size() << " }";
            throw std::domain_error(ss.str());
        }

        // Callback
        if (iteration % callback_frequency == 0) callback(state);
    }
}

} // namespace sde


// --- ODE solvers ---
// -------------------
namespace ode {

enum class Method {
    EXPLICIT_EULER,  // explicit, convergence order = 1
    IMPLICIT_EULER,  // implicit, convergence order = 1
    SYMMETRIC_EULER, // implicit, convergence order = 1, conservative
    RK2,             // explicit, convergence order = 2
    RK4              // explicit, convergence order = 4
};

struct State {
    const Scalar& t;          // current time
    const uint&   iteration;  // current iteration
    const Vector& y0;         // current x
    const Vector& y_computed; // 'next' x
};

using CallbackFunction = std::function<void(State&)>;

// - Solver implementations -

inline void _step_explicit_euler(VectorTimeFunction& F, Scalar t, Scalar tau, Vector& y0, Vector& y) {
    y = y0 + tau * F(t, y0);
}

inline void _step_implicit_euler(VectorTimeFunction& F, Scalar t, Scalar tau, Vector& y0, Vector& y) {

    VectorFunction implicit_equation = [&](const Vector& yn) -> Vector { return yn - y0 - tau * F(t + tau, yn); };

    y = algebraic::nonlinear_solve(implicit_equation, y0, 1e-8); // NOTE: Harcoded solver precision, kinda bad
}

inline void _step_symmetric_euler(VectorTimeFunction& F, Scalar t, Scalar tau, Vector& y0, Vector& y) {

    VectorFunction implicit_equation = [&](const Vector& yn) -> Vector {
        return yn - y0 - tau * 0.5 * (F(t + tau, yn) + F(t, y0));
    };

    y = algebraic::nonlinear_solve(implicit_equation, y0, 1e-8); // NOTE: Harcoded solver precision, kinda bad
}

inline void _step_rk2(VectorTimeFunction& F, Scalar t, Scalar tau, Vector& y0, Vector& y) {
    const Vector k1 = F(t, y0);
    const Vector k2 = F(t + tau, y0 + tau * k1);

    y = y0 + tau / 2. * (k1 + k2);
}

inline void _step_rk4(VectorTimeFunction& F, Scalar t, Scalar tau, Vector& y0, Vector& y) {
    const Vector k1 = F(t, y0);
    const Vector k2 = F(t + 0.5 * tau, y0 + 0.5 * tau * k1);
    const Vector k3 = F(t + 0.5 * tau, y0 + 0.5 * tau * k2);
    const Vector k4 = F(t + tau, y0 + tau * k3);

    y = y0 + tau / 6. * (k1 + 2. * k2 + 2. * k3 + k4);
}

// Solve ordinary diff. equation (ODE) with form:
// dX(t)/dt = F(t, X(t))
inline void ode_solve(Method             method,          // solver
                      VectorTimeFunction F,               // determenistic part
                      Vector             y0,              // initial condition
                      Scalar t1, Scalar t2,               // time interval
                      Scalar           tau,               // time step
                      CallbackFunction callback,          // callback function
                      uint             callback_frequency // callback frequency (in iterations)
) {
    // 'x0' passed by value cause we use it as an internal variable

    // Iteration
    Vector y         = y0;
    Scalar t         = t1;
    uint   iteration = 0;

    State state{t, iteration, y0, y};

    for (t = t1; t < t2; t += tau, ++iteration) {
        // Get 'x' from 'x0' using appropriate solver
        y0 = y;

        switch (method) {
        case Method::EXPLICIT_EULER: _step_explicit_euler(F, t, tau, y0, y); break;
        case Method::IMPLICIT_EULER: _step_implicit_euler(F, t, tau, y0, y); break;
        case Method::SYMMETRIC_EULER: _step_symmetric_euler(F, t, tau, y0, y); break;
        case Method::RK2: _step_rk2(F, t, tau, y0, y); break;
        case Method::RK4: _step_rk4(F, t, tau, y0, y); break;
        // TODO: Other methods
        default: return;
        }

        // Check for solution divergence
        for (idx k = 0; k < y.size(); ++k) {
            // Stop and throw if any values in 'y' are 'Inf' or 'NaN'
            if (std::isfinite(y[k])) continue;
            std::stringstream ss;
            ss << "Inf or NaN encountered during ODE solution at { iteration = " << iteration << ", t = " << t
               << " } on index { k = " << k << " / " << y.size() << " }";
            throw std::domain_error(ss.str());
        }

        // Callback
        if (iteration % callback_frequency == 0) callback(state);
    }
}

} // namespace ode

} // namespace gsse