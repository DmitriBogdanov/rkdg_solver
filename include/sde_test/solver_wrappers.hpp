#pragma once

#include "firstparty/gsse.hpp"

#include <fstream>
#include <random>
#include <string>

using gsse::Vector;
using gsse::VectorTimeFunction;



// Wiener's process
inline void wieners_process(double T, double tau, std::string savepath) {
    std::ofstream file_t(savepath + "[t].txt");
    std::ofstream file_y(savepath + "[X].txt");

    std::mt19937                   gen(0);
    std::uniform_real_distribution dist;

    const auto sqrt_of_tau = std::sqrt(tau);

    double W = 0;

    for (double t = 0; t < T; t += tau) {
        // Save results
        file_t << t << '\n';
        file_y << W << '\n';

        // Get new value
        const double dW = dist(gen) * sqrt_of_tau;
        W += dW;
    }
}

// Explicit Euler's method
inline void sodesolve_euler_maruyama(double T, double tau, Vector y0, VectorTimeFunction A, VectorTimeFunction B,
                                     std::string savepath) {
    std::ofstream file_t(savepath + "[t].txt");
    std::ofstream file_y(savepath + "[X].txt");

    Vector y = y0;

    std::mt19937                   gen(0);
    std::uniform_real_distribution dist;

    const auto sqrt_of_tau = std::sqrt(tau);

    for (double t = 0; t < T; t += tau) {
        // Save results
        file_t << t << '\n';
        for (const auto& coord : y) file_y << coord << ' ';
        file_y << '\n';

        // Get new values
        y0 = y;

        const double dW = dist(gen) * sqrt_of_tau;

        /// Euler-Maruyama
        // Has a strong order of 0.5 and a weak order of 1
        // For ODE's equvalent to explicit Euler method
        y = y0 + tau * A(t, y0) + dW * B(t, y0);
    }
}


// // Runge-Kutta's method
// void sodesolve_RK(double T, double tau, Vector y0, VectorTimeFunction A, VectorTimeFunction B, std::string savepath)
// { 	std::ofstream file_t(savepath + "[t].txt"); 	std::ofstream file_y(savepath + "[X].txt");

// 	NormalDistrGenerator generator;
// 	generator.seed(0);

// 	const auto sqrt_of_tau = sqrt(tau);

// 	// Iteration
// 	Vector y = y0;

// 	double W = 0;

// 	for (double t = 0; t < T; t += tau) {
// 		// Save results
// 		file_t << t << '\n';
// 		for (const auto &coord : y) file_y << coord << ' ';
// 		file_y << '\n';

// 		// Get new values
// 		y0 = y;

// 		const double dW = generator.get_next() * sqrt_of_tau;

// 		/// Runge-Kutta from "[P.E. Kloeden, E. Platen] Numerical solution of stochastic differential equations" (page
// 374)
// 		/// written for Ito process A = A(X), B = B(X)
// 		// Has a strong order of 1 and a weak order of 1
// 		// For ODE's equvalent to explicit Euler method
// 		const Vector Y = y0 + tau * A(t, y0) + sqrt_of_tau * B(t, y0);

// 		y = y0 + tau * A(t, y0) + dW * B(t, y0) + 0.5 * (B(t, Y) - B(t, y0)) * (dW * dW - tau) / sqrt_of_tau;

// 		/// "Intended" Runge-Kutta
// 		//constexpr double
// 		//	R21 = 2. / 3.,
// 		//	R31 = -1.,
// 		//	R32 = 1.; // other components are zero

// 		//constexpr double
// 		//	r2 = 3. / 4.,
// 		//	r3 = 1. / 4.; // other components are zero

// 		//const Vector x1 =
// 		//	y0;


// 		//const Vector x2 =
// 		//	y0
// 		//	+ tau * (R21 * A(t, x1))
// 		//	+ dW  * (R21 * B(t, x1));


// 		//const Vector x3 =
// 		//	y0
// 		//	+ tau * (R31 * A(t, x1) + R32 * A(t, x2))
// 		//	+ dW  * (R31 * B(t, x1) + R32 * B(t, x2));

// 		//y =
// 		//	y0
// 		//	+ tau * (r2 * A(t, x2) + r3 * A(t, x3))
// 		//	+ dW  * (r2 * B(t, x2) + r3 * B(t, x3));

// 		/// Wikipedia "Runge-Kutta as Improved Euler"
// 		/*const Vector k1 = tau * A(t, y0) + dW * B(t, y0);

// 		const Vector k2 = tau * A(t + tau, y0 + k1) + dW * B(t + tau, y0 + k1);

// 		y = y0 + 0.5 * (k1 + k2);*/

// 		/// Buther's table (3) from a book
// 		//constexpr double
// 		//	R21 = 1. / 2.,
// 		//	R31 = -1.,
// 		//	R32 = 2.; // other components are zero

// 		//constexpr double
// 		//	r1 = 1. / 6.,
// 		//	r2 = 4. / 6.,
// 		//	r3 = 1. / 6.; // other components are zero

// 		//const double t1 =
// 		//	t;

// 		//const Vector x1 =
// 		//	y0;

// 		//const double t2 =
// 		//	t
// 		//	+ tau * 0.5;

// 		//const Vector x2 =
// 		//	y0
// 		//	+ tau * (R21 * A(t1, x1))
// 		//	+ dW  * (R21 * B(t1, x1));

// 		//const double t3 =
// 		//	t
// 		//	+ tau * 1.0;

// 		//const Vector x3 =
// 		//	y0
// 		//	+ tau * (R31 * A(t1, x1) + R32 * A(t2, x2))
// 		//	+ dW  * (R31 * B(t1, x1) + R32 * B(t2, x2));

// 		//y =
// 		//	y0
// 		//	+ tau * (r1 * A(t1, x1) + r2 * A(t, x2) + r3 * A(t, x3))
// 		//	+ dW  * (r1 * B(t1, x1) + r2 * B(t, x2) + r3 * B(t, x3));


// 		/// Adaptation of a regular RK2
// 		/*const auto k1 = A(t, y0) + B(t, y0) * (W - dW);
// 		const auto k2 = A(t + tau, y0 + tau * k1) + B(t + tau, y0 + tau * k1) * W;

// 		y = y0 + tau * (0.5 * k1 + 0.5 * k2);*/

// 		/*const Vector a1 = A(t, y0);
// 		const Vector b1 = (W - dW) * B(t, y0);
// 		const Vector k1 = a1 + b1;

// 		const Vector a2 = A(t + tau, y0 + tau * k1);
// 		const Vector b2 = W * B(t + tau, y0 + tau * k1);
// 		const Vector k2 = a2 + b2;

// 		std::cout
// 			<< "t = " << t << "\n"
// 			<< "k1 = " << (k1).format(INLINE) << "\n"
// 			<< "k2 = " << (k2).format(INLINE) << "\n"
// 			<< "a1 = " << (a1).format(INLINE) << "\n"
// 			<< "a2 = " << (a2).format(INLINE) << "\n"
// 			<< "b1 = " << (b1).format(INLINE) << "\n"
// 			<< "b2 = " << (b2).format(INLINE) << "\n"
// 			<< "A(t, y0) = " << (A(t, y0)).format(INLINE) << "\n"
// 			<< "A(t + tau, y0 + tau * k1) = " << (A(t + tau, y0 + tau * k1)).format(INLINE) << "\n"
// 			<< "B(t, y0) * (W - dW) = " << (B(t, y0) * (W - dW)).format(INLINE) << "\n"
// 			<< "B(t + tau, y0 + tau * k1) * W = " << (B(t + tau, y0 + tau * k1) * W).format(INLINE) << "\n";

// 		y = y0 + tau * (0.5 * k1 + 0.5 * k2);*/
// 	}
// }