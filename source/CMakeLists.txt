# Macro for defining a test with proper doctest args
# Note: Adding '--force-colors' to 'add_test()' allows doctest coloring in terminal but messes up CTest logs formatting
macro(add_solver_executable filename)
    add_executable(${filename} ${filename}.cpp)
    target_compile_features(${filename} PRIVATE cxx_std_17)
    target_compile_options(${filename} PRIVATE -O2 -Wall -Wextra -Wpedantic -Werror -fmax-errors=0)
    # Add for more rigorous checks:           -Wall -Wextra -Wpedantic -Werror
    # Add for more debug info when profiling: -g -no-pie
    # Add when using OpenMP:                  -fopenmp
endmacro()

add_solver_executable(prng_test)
add_solver_executable(rkdg)
add_solver_executable(sde_test)