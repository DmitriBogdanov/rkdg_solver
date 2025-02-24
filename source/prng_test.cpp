#include "prng_test/prng_test_cases.hpp"



int main() {
    using namespace utl;
    namespace gens = utl::random::generators;

    // --- Grayscale test ---
    // ----------------------

    std::filesystem::remove_all(prng_grayscale_directory);
    std::filesystem::create_directories(prng_grayscale_directory);

    // 'std' PRNGs
    prng_test_grayscale("minstd_rand", std::minstd_rand{});
    prng_test_grayscale("mt19937", std::mt19937{});

    // 'bad' PRNGs
    const auto compiler_name =
        stre::replace_all_occurences(stre::to_lower(std::string(predef::compiler_name)), " ", "_") + "_";
    prng_test_grayscale(compiler_name + "rand", CRandWrapper{}); // rand() is platfrom-dependent, this is important
    prng_test_grayscale("RANDU", RANDU{});

    // 'utl' PRNGs
    prng_test_grayscale("RomuMono16", gens::RomuMono16{});
    prng_test_grayscale("RomuTrio32", gens::RomuTrio32{});
    prng_test_grayscale("SplitMix32", gens::SplitMix32{});
    prng_test_grayscale("Xoshiro128PP", gens::Xoshiro128PP{});
    prng_test_grayscale("RomuDuoJr64", gens::RomuDuoJr64{});
    prng_test_grayscale("SplitMix64", gens::SplitMix64{});
    prng_test_grayscale("Xoshiro256PP", gens::Xoshiro256PP{});
    prng_test_grayscale("ChaCha12", gens::ChaCha12{});

    // --- 3D points test ---
    // ----------------------

    std::filesystem::remove_all(prng_3d_points_directory);
    std::filesystem::create_directories(prng_3d_points_directory);

    generate_random_3d_points("RANDU", RANDU{});
    generate_random_3d_points("Xoshiro256PP", gens::Xoshiro256PP{});

    // --- Normal distr test ---
    // -------------------------

    std::filesystem::remove_all(prng_normal_dist_directory);
    std::filesystem::create_directories(prng_normal_dist_directory);

    generate_normal_distribution("minstd_rand", std::minstd_rand{});
    generate_normal_distribution("Xoshiro256PP", gens::Xoshiro256PP{});
}