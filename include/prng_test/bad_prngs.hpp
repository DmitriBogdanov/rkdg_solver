#pragma once

#include <cstdlib>
#include <random>



// An example of how NOT to do PRNGs (1)
class RANDU {
public:
    using result_type = long int;

private:
    result_type s{};

public:
    constexpr static result_type a = (result_type(1) << 16) + result_type(3);
    constexpr static result_type m = (result_type(1) << 31);

public:
    constexpr explicit RANDU(result_type seed = 1) noexcept { this->seed(seed); }

    [[nodiscard]] static constexpr result_type min() noexcept { return 1; }
    [[nodiscard]] static constexpr result_type max() noexcept { return m - 1; }

    constexpr void seed(result_type seed) noexcept { this->s = seed; }

    constexpr result_type operator()() noexcept {
        this->s = (a * this->s) % m;
        return this->s;
    }
};

// An example of how NOT to do PRNGs (2)
class CRandWrapper {
public:
    using result_type = int;

public:
    explicit CRandWrapper(result_type seed = 1) noexcept { this->seed(seed); }

    [[nodiscard]] static constexpr result_type min() noexcept { return 0; }
    [[nodiscard]] static constexpr result_type max() noexcept { return RAND_MAX; }

    void seed(result_type seed) noexcept { std::srand(seed); }

    result_type operator()() noexcept { return std::rand(); }
};