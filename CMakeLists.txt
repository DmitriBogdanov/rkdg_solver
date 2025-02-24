cmake_minimum_required(VERSION 3.9.1)

project(
    RKDG_SOLVER
    VERSION 1.0.1
    DESCRIPTION "Partial codebase for a thesis in stochastic modeling."
    HOMEPAGE_URL "https://github.com/DmitriBogdanov/rkdg_solver"
)

# Source
include_directories(include)

# Executables & Tests
add_subdirectory(source)
add_subdirectory(tests)