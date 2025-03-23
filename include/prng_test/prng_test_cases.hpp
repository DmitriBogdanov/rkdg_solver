#pragma once

#include <random>

#include "firstparty/UTL.hpp"

#include "bad_prngs.hpp"



// =============================
// --- Grayscale matrix test ---
// =============================

const std::string prng_grayscale_directory = "temp/prng_test_grayscale/";

template <class Generator>
void prng_test_grayscale(const std::string prng_name, Generator gen) {
    using namespace utl;

    constexpr std::size_t image_rows = 160;
    constexpr std::size_t image_cols = 200;
    constexpr std::size_t block_size = 1;

    constexpr std::size_t                       color_range = 1;
    std::uniform_int_distribution<unsigned int> dist(0, color_range);

    mvl::Matrix<uint8_t> grayscale_image(image_rows, image_cols, [&](std::size_t i, std::size_t j) {
        if (i % block_size == 0 && j == 0) gen.seed(1 + i);
        return dist(gen);
    });

    const std::string path = prng_grayscale_directory + prng_name + ".txt";
    UTL_LOG_TRACE("PRNG test: Grayscale matrix. PRNG -> ", prng_name, ". path -> ", path);
    std::ofstream(path) << mvl::format::as_raw(grayscale_image);
}

// ======================
// --- 3D points test ---
// ======================

const std::string prng_3d_points_directory = "temp/prng_test_3d_points/";

template <class Generator>
void generate_random_3d_points(const std::string prng_name, Generator gen) {
    constexpr std::uint32_t        seed = 15;
    gen.seed(seed);

    constexpr std::size_t point_count = 10'000; // we generate N^3 points
    
    const std::string path = prng_3d_points_directory + prng_name + ".txt";
    UTL_LOG_TRACE("PRNG test: 3D points. PRNG -> ", prng_name, ". path -> ", path);
    std::ofstream file(path);
    // Plot random points in 3D space,
    // Notice that we specifically use "divide by max" formula instead of 'std::unifor_real_distribution'
    // because std distribution does so smart stuff that makes RANDU badness less apparent.
    const double factor = 1. / (decltype(gen)::max() - decltype(gen)::min());
    for (std::size_t i = 0; i < point_count; ++i) file << gen() * factor << ' ' << gen() * factor << ' ' << gen() * factor << '\n';
}

// =========================
// --- Normal distr test ---
// =========================

const std::string prng_normal_dist_directory = "temp/prng_test_normal_dist/";

template <class Generator>
void generate_normal_distribution(const std::string prng_name, Generator gen) {
    constexpr std::uint32_t        seed = 14;
    gen.seed(seed);
    
    constexpr std::size_t point_count = 1'000; // we generate N^3 points
    
    std::normal_distribution<double> dist;
    
    const std::string path = prng_normal_dist_directory + prng_name + ".txt";
    UTL_LOG_TRACE("PRNG test: Normal distribution. PRNG -> ", prng_name, ". path -> ", path);
    std::ofstream file(path);
    for (std::size_t i = 0; i < point_count; ++i) file << dist(gen) << '\n';
}