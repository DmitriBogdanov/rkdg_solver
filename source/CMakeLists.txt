add_executable(run main.cpp)

target_compile_features(run PRIVATE cxx_std_17)
target_compile_options(run PRIVATE -O2 -Wall -Wextra -Wpedantic -Werror)
# Add for more rigorous checks:           -Wall -Wextra -Wpedantic -Werror
# Add for more debug info when profiling: -g -no-pie