#include "rkdg/rkdg_convdiff.hpp"
#include "rkdg/rkdg_generic.hpp"



int main() {
    using namespace utl;
    std::cout << "\n" << timer::datetime_string() << "\n";
    timer::start();

    // Solve
    try {
        rkdg_convdiff::solve_eq();
    } catch (const std::exception& except) { std::cerr << "ERROR: Caught exception!\n" << except.what() << "\n"; }

    // std::cout << std::ifstream("temp/log.txt").rdbuf();

    gsse::Matrix mat = gsse::make_matrix({
        {1, 2, 3},
        {4, 5, 6},
    });

    std::cout << "\nDone!\n"
              << "Elapsed time: " << timer::elapsed_string_fullform() << "\n"
              << timer::datetime_string() << "\n\n";
}