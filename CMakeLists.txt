cmake_minimum_required(VERSION 3.9.1)

project(project_name)

# Source
include_directories(source)
add_subdirectory(source)

# Tests
add_subdirectory(tests)