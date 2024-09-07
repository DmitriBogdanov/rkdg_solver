# ---------------------------------- CONTENTS -----------------------------------
#   
#   This script contains shortcuts for building, running
#   and testing the project. All action keywords can be
#   chained which causes them to be executed one after another.
#
#   See "docs/guide_building_project.md" for the whole building guide.
#   
# ------------------------------------ GUIDE ------------------------------------
#   
#   ---- Usage format ----
#     > bash actions.sh [ACTIONS]
#   
#   ------ Actions -------
#     clear   - Clears "build/" folder
#     config  - Configures CMake with appropriate args                                [requires installed cmake]
#     build   - Builds the project                     (requires successful 'config')
#     run     - Runs main executable                   (requires successful 'build' )
#     profile - Runs main executable with profiler     (requires successful 'build' ) [requires installed callgrind & kcachegrind]
#     test    - Runs CTest tests                       (requires successful 'build' )
#   
#   --- Usage examples ---
#     # Example 1: clear existing build & temp files -> configure cmake -> cmake build -> run main executable
#       > bash actions.sh clear config build run
#
#     # Example 2: cmake build -> run main executable with callgrind profiler and visualize dump with kcachegrind
#       > bash actions.sh build profile
#   
# -------------------------------------------------------------------------------

# -----------------------
# ---- Configuration ----
# -----------------------
directory_source="source/"
directory_build="build/"
directory_tests="${directory_build}tests/"
directory_temp="temp/"

path_executable="${directory_build}source/run"

compiler="g++"                                               # Selected compiler
test_flags="--rerun-failed --output-on-failure --timeout 60" # CTest flags
build_jobs="6"                                               # CMake build jobs

# -----------------------
# ------ Functions ------
# -----------------------
check_command_exists() {
    if ! command -v $1 &> /dev/null
    then
        echo "Command [ $1 ] could not be found."
        exit 1
    fi
}

clear_files() {
    if [ -d "$directory_build" ]; then
        rm --recursive $directory_build
        echo "Cleared directory [ $directory_build ]."
    else
        echo "Directory [ $directory_build ] is clear."
    fi
}

cmake_config() {
    check_command_exists "cmake"
    check_command_exists "$compiler"
    cmake -D CMAKE_CXX_COMPILER=$compiler -B $directory_build -S .
}

cmake_build() {
    check_command_exists "cmake"
    cmake --build $directory_build --parallel $build_jobs
}

cmake_test() {
    check_command_exists "ctest"
    cd $directory_tests
    ctest $test_flags
    cd ..
}

executable_run() {
    ./$path_executable
}

executable_profile() {
    check_command_exists "valgrind"
    check_command_exists "kcachegrind"
    valgrind --tool=callgrind --callgrind-out-file="${directory_temp}callgrind.latest" ./$path_executable
    kcachegrind "./${directory_temp}callgrind.latest"
}

# -----------------------
# --- Action selector ---
# -----------------------
valid_command=false

for var in "$@"
do
    valid_command=false
    
    if [ "$var" = "clear" ]; then
        echo "# Action: Clear Files"
        clear_files
        valid_command=true
    fi

    if [ "$var" = "config" ]; then
        echo "# Action: CMake Configure"
        cmake_config
        valid_command=true
    fi

    if [ "$var" = "build" ]; then
        echo "# Action: CMake Build"
        cmake_build
        valid_command=true
    fi
    
    if [ "$var" = "test" ]; then
        echo "# Action: CMake Test"
        cmake_test
        valid_command=true
    fi

    if [ "$var" = "run" ]; then
        echo "# Action: run"
        executable_run
        valid_command=true
    fi
    
    if [ "$var" = "profile" ]; then
        echo "# Action: profile"
        executable_profile
        valid_command=true
    fi
    
    if [ $valid_command = false ]; then
        echo "# Error: Invalid action name -> ${var}"
        break
    fi

done