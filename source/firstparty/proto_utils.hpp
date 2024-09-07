// ~~~~~~~~~~~~~~~~~~~~~~ DmitriBogdanov/prototyping_utils ~~~~~~~~~~~~~~~~~~~~~~
//
// MIT License
//
// Copyright (c) 2023 Dmitri Bogdanov
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~






// -------------------------------
// --------- utl::config ---------
// -------------------------------
#include <cstddef>
#include <random>
#if !defined(UTL_PICK_MODULES) || defined(UTLMODULE_CONFIG)
#ifndef UTLHEADERGUARD_CONFIG
#define UTLHEADERGUARD_CONFIG

// ======== header guard start ========

#include <fstream>
#include <initializer_list>
#include <ios>
#include <iterator>
#include <ostream>
#include <string>
#include <string_view>
#include <tuple>

namespace utl::config {

// _________ DEVELOPER DOCS _________

// Utils for exporting/importing configs.
// Currently suppors JSON.
//
// # ::export_json() #
// Exports JSON with given key-value entries. Entries support homogenous arrays with any level of nesting,
// (after 4 levels some additional boilerplate is needed but who need 4D+ json matrices). Object entries and
// non-homogenous (type-wise) arrays are not supported.
//
// # ::entry() #
// A clutch that helps compiler resolve template types in brace-initialized arguments.
// > export_json(..., { "key", { values } })             // 'value' type doesn't get resolved
// > export_json(..., config::entry("key", { values }))  // 'value' type gets resolved

// _________ IMPLEMENTATION _________

// --- Export config ---
template <typename T>
inline void _write_formatted_value_to_ostream(std::ostream& os, const T& value) {
    constexpr bool is_string = std::is_convertible<T, std::string_view>::value;
    constexpr bool is_array  = !std::is_scalar<T>::value;

    // Array => Expand recursively, add brackets and comma separators
    if constexpr (is_array) {
        os << "[ ";
        for (auto it = value.begin(); it != value.end(); ++it) {
            _write_formatted_value_to_ostream(os, *it);
            if (std::next(it) != value.end()) os << ", "; // prevents trailing comma
        }
        os << " ]";
    }
    // String => Add quotes
    else if constexpr (is_string) {
        os << "\"" << value << "\"";
    }
    // Numeric/Bool => Nothing
    else {
        os << value;
    }
}

inline void _write_entries(std::ostream& os) { os << std::flush; };

template <typename T, typename... TupleTypes>
void _write_entries(std::ostream& os, std::tuple<std::string_view, T> entry, const TupleTypes... other_entries) {
    // Write one entry
    const auto& key   = std::get<0>(entry);
    const auto& value = std::get<1>(entry);

    os << "    \"" << key << "\": ";
    _write_formatted_value_to_ostream(os, value);
    if constexpr (sizeof...(other_entries) != 0) os << ","; // prevents trailing comma
    os << "\n";

    _write_entries(os, other_entries...);
}

template <typename... TupleTypes>
void export_json(std::string_view path, const TupleTypes... entries) {
    std::string   str(path);
    std::ofstream file(str);
    file << std::boolalpha; // in JSON bools are written as true/false

    file << "{\n";
    _write_entries(file, entries...);
    file << "}";
}

template <typename T>
std::tuple<std::string_view, T> entry(std::string_view key, T value) {
    return {key, value};
    // may feel like a redundancy but without it there is no way for compiler to deduce std::tuple<>
}

// - std::initializer_list resolution overloads -
// Contains overloads that specify that brace-initialization should be interpreted as std::initializer_list<>.
// Up to 4 levels of nesting are specified, after that you do that manually on the call site.
// Also I can't imagine an adequate person saving 5D arrays to JSON.
template <typename T>
using _il = const std::initializer_list<T>; // shortcut name to make nested 'std::initializer_list' less verbose
                                            // for some reason without 'const' GCC & clang have issues deducing template
                                            // type in entry(), while MSVC performs type deduction as one would expect

template <typename T>
using _ret = std::tuple<std::string_view, _il<T>>;

// 1D
template <typename T>
std::tuple<std::string_view, _il<T>> entry(std::string_view key, _il<T> value) {
    return {key, value};
}

// 2D
template <typename T>
std::tuple<std::string_view, _il<_il<T>>> entry(std::string_view key, _il<_il<T>> value) {
    return {key, value};
}

// 3D
template <typename T>
std::tuple<std::string_view, _il<_il<_il<T>>>> entry(std::string_view key, _il<_il<_il<T>>> value) {
    return {key, value};
}

// 4D
template <typename T>
std::tuple<std::string_view, _il<_il<_il<_il<T>>>>> entry(std::string_view key, _il<_il<_il<_il<T>>>> value) {
    return {key, value};
}

// ... we can continue so on and so forth but it isn't really necessary
// if someone needs more layers of nesting they can specify type resolutions with std::initializer_list{ ... }

// TODO:
// --- Import config ---
// template<typename T>
// std::tuple<std::string_view, T&> parse(std::string_view key, T &value) {
//	return { key, value };
//	// may feel like a redundancy but without it there is no way for compiler to deduce std::tuple<>
//}

// inline void _read_entries(std::istream &is) {
//
// };

// template<typename T, typename... TupleTypes>
// void _read_entries(std::istream &is, std::tuple<std::string_view, T&> entry, const TupleTypes... other_entries) {
//	// Read one entry
//	std::string key = std::get<0>(entry);
//	auto &value = std::get<1>(entry);

//	std::cout << typeid(value).name() << std::endl;

//	is.ignore('\"');
//	is >> key;
//	is.ignore(':');
//	is >> value;
//	is.ignore(',');

//	_read_entries(is, other_entries...);
//}

// template<typename... TupleTypes>
// void import_json(std::string_view path, const TupleTypes... entries) {
//	std::ifstream file(path);
//
//	_read_entries(file, entries...);
// }

} // namespace utl::config

// ========= header guard end =========

#endif
#endif






// -----------------------------
// --------- utl::math ---------
// -----------------------------
#if !defined(UTL_PICK_MODULES) || defined(UTLMODULE_MATH)
#ifndef UTLHEADERGUARD_MATH
#define UTLHEADERGUARD_MATH

// ======== header guard start ========

#include <cassert>
#include <cstdint>
#include <functional>
#include <type_traits>
#include <utility>
#include <vector>

namespace utl::math {

// _________ DEVELOPER DOCS _________

// Coordinate transformations, mathematical constants and technical helper functions.
//
// # ::PI, ::PI_TWO, ::PI_HALF, ::E, ::GOLDEN_RATION #
// Constants.
//
// # ::is_addable_with_itself<Type> #
// Integral constant, returns in "::value" whether Type supports 'operator()+' with itself.
//
// # ::is_multipliable_by_scalar<Type> #
// Integral constant, returns in "::value" whether Type supports 'operator()*' with double.
//
// # ::is_sized<Type> #
// Integral constant, returns in "::value" whether Type supports '.size()' method.
//
// # ::abs(), ::sign(), ::sqr(), ::cube(), ::midpoint(), deg_to_rad(), rad_to_deg() #
// Constexpr templated math functions, useful when writing expressions with a "textbook form" math.
//
// # ::uint_difference() #
// Returns abs(uint - uint) with respect to uint size and possible overflow.
//
// # ::linspace() #
// Tabulates [min, max] range with N evenly spaced points and returns it as a vector.
//
// # ::ssize() #
// Returns '.size()' of the argument casted to 'int'.
// Essentially a shortcut for verbose 'static_cast<int>(container.size())'.
//
// # ::ternary_branchless() #
// Branchless ternary operator. Slightly slower that regular ternary on most CPUs.
// Should not be used unless branchess qualifier is necessary (like in GPU computation).
//
// # ::ternary_bitselect() #
// Faster branchless ternary for integer types.
// If 2nd return is ommited, 0 is assumed, which allows for significant optimization.

// _________ IMPLEMENTATION _________

// --- Constants ---
constexpr double PI           = 3.14159265358979323846;
constexpr double PI_TWO       = 2. * PI;
constexpr double PI_HALF      = 0.5 * PI;
constexpr double E            = 2.71828182845904523536;
constexpr double GOLDEN_RATIO = 1.6180339887498948482;

// --- type trait: is_addable_with_itself ---
template <typename Type, typename = void>
struct is_addable_with_itself : std::false_type {};

template <typename Type>
struct is_addable_with_itself<Type, std::void_t<decltype(std::declval<Type>() + std::declval<Type>())>
                              // perhaps check that resulting type is same as 'Type', but that can cause issues
                              // with classes like Eigen::MatrixXd that return "foldables" that convert back to
                              // objects in the end
                              > : std::true_type {};

// --- type trait: is_multipliable_by_scalar ---
template <typename Type, typename = void>
struct is_multipliable_by_scalar : std::false_type {};

template <typename Type>
struct is_multipliable_by_scalar<Type, std::void_t<decltype(std::declval<Type>() * std::declval<double>())>>
    : std::true_type {};

// --- type trait: is_sized ---
template <typename Type, typename = void>
struct is_sized : std::false_type {};

template <typename Type>
struct is_sized<Type, std::void_t<decltype(std::declval<Type>().size())>> : std::true_type {};

// --- type trait: is_function_with_signature ---
template <typename FuncType, typename Signature>
using is_function_with_signature = std::is_convertible<FuncType, std::function<Signature>>;

// --- Standard math functions ---
template <typename Type, std::enable_if_t<std::is_scalar<Type>::value, bool> = true>
constexpr Type abs(Type x) {
    return (x > Type(0)) ? x : -x;
}

template <typename Type, std::enable_if_t<std::is_scalar<Type>::value, bool> = true>
constexpr Type sign(Type x) {
    return (x > Type(0)) ? Type(1) : Type(-1);
}

template <typename Type, std::enable_if_t<std::is_arithmetic<Type>::value, bool> = true>
constexpr Type sqr(Type x) {
    return x * x;
}

template <typename Type, std::enable_if_t<std::is_arithmetic<Type>::value, bool> = true>
constexpr Type cube(Type x) {
    return x * x * x;
}

template <typename Type, std::enable_if_t<utl::math::is_addable_with_itself<Type>::value, bool> = true,
          std::enable_if_t<utl::math::is_multipliable_by_scalar<Type>::value, bool> = true>
constexpr Type midpoint(Type a, Type b) {
    return (a + b) * 0.5;
}

template <typename IntegerType, std::enable_if_t<std::is_integral<IntegerType>::value, bool> = true>
constexpr int kronecker_delta(IntegerType i, IntegerType j) {
    // 'IntegerType' here is necessary to prevent enforcing static_cast<int>(...) on the callsite
    return (i == j) ? 1 : 0;
}

template <typename IntegerType, std::enable_if_t<std::is_integral<IntegerType>::value, bool> = true>
constexpr int power_of_minus_one(IntegerType power) {
    return (power % IntegerType(2)) ? -1 : 1; // is there a faster way of doing it?
}


// --- deg/rad conversion ---
template <typename FloatType, std::enable_if_t<std::is_floating_point<FloatType>::value, bool> = true>
constexpr FloatType deg_to_rad(FloatType degrees) {
    constexpr FloatType FACTOR = FloatType(PI / 180.);
    return degrees * FACTOR;
}

template <typename FloatType, std::enable_if_t<std::is_floating_point<FloatType>::value, bool> = true>
constexpr FloatType rad_to_deg(FloatType radians) {
    constexpr FloatType FACTOR = FloatType(180. / PI);
    return radians * FACTOR;
}


// --- Meshing ---

// Semantic helpers that allow user to directly pass both interval/point counts for grid subdivision,
// without thinking about whether function need +1 or -1 to its argument
struct Points {
    std::size_t count;

    Points() = delete;
    explicit Points(std::size_t count) : count(count) {}
};

struct Intervals {
    std::size_t count;

    Intervals() = delete;
    explicit Intervals(std::size_t count) : count(count) {}
    Intervals(Points points) : count(points.count - 1) {}
};

template <typename FloatType, std::enable_if_t<std::is_floating_point<FloatType>::value, bool> = true>
std::vector<FloatType> linspace(FloatType L1, FloatType L2, Intervals N) {
    assert(L1 < L2);
    assert(N.count >= 1);

    const FloatType step = (L2 - L1) / N.count;

    std::vector<FloatType> res(N.count + 1);

    res[0] = L1;
    for (std::size_t i = 1; i < res.size(); ++i) res[i] = res[i - 1] + step;

    return res;
}

template <typename FloatType, typename FuncType,
          std::enable_if_t<std::is_floating_point<FloatType>::value, bool>                          = true,
          std::enable_if_t<is_function_with_signature<FuncType, FloatType(FloatType)>::value, bool> = true>
FloatType integrate_trapezoidal(FuncType f, FloatType L1, FloatType L2, Intervals N) {
    assert(L1 < L2);
    assert(N.count >= 1);

    const FloatType step = (L2 - L1) / N.count;

    FloatType sum = 0;
    FloatType x   = L1;

    for (std::size_t i = 0; i < N.count; ++i, x += step) sum += f(x) + f(x + step);

    return FloatType(0.5) * sum * step;
}

// --- Misc helpers ---
template <typename UintType, std::enable_if_t<std::is_integral<UintType>::value, bool> = true>
constexpr UintType uint_difference(UintType a, UintType b) {
    // Cast to widest type if there is a change values don't fit into a regular 'int'
    using WiderIntType = std::conditional_t<(sizeof(UintType) >= sizeof(int)), int64_t, int>;

    return static_cast<UintType>(utl::math::abs(static_cast<WiderIntType>(a) - static_cast<WiderIntType>(b)));
}

template <typename SizedContainer, std::enable_if_t<utl::math::is_sized<SizedContainer>::value, bool> = true>
int ssize(const SizedContainer& container) {
    return static_cast<int>(container.size());
}

// Branchless ternary
template <typename Type, std::enable_if_t<std::is_arithmetic<Type>::value, bool> = true>
constexpr Type ternary_branchless(bool condition, Type return_if_true, Type return_if_false) {
    return (condition * return_if_true) + (!condition * return_if_false);
}

template <typename IntType, std::enable_if_t<std::is_integral<IntType>::value, bool> = true>
constexpr IntType ternary_bitselect(bool condition, IntType return_if_true, IntType return_if_false) {
    return (return_if_true & -IntType(condition)) | (return_if_false & ~(-IntType(condition)));
}

template <typename IntType, std::enable_if_t<std::is_integral<IntType>::value, bool> = true>
constexpr IntType ternary_bitselect(bool condition, IntType return_if_true) {
    return return_if_true & -IntType(condition);
}

} // namespace utl::math

// ========= header guard end =========

#endif
#endif






// ------------------------------------
// --------- utl::progressbar ---------
// ------------------------------------
#if !defined(UTL_PICK_MODULES) || defined(UTLMODULE_PROGRESSBAR)
#ifndef UTLHEADERGUARD_PROGRESSBAR
#define UTLHEADERGUARD_PROGRESSBAR

// ======== header guard start ========

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <ostream>
#include <sstream>
#include <string>

namespace utl::progressbar {

// _________ DEVELOPER DOCS _________

// Simple progress bars for terminal applications.
//
// # ::set_ostream() #
// Sets ostream used for progress bars.
//
// # ::Percentage #
// Proper progress bar, uses carriage return escape sequence (\r) to render new states in the same spot.
// Shows an estimate of remaining time.
//
// # ::Ruler #
// Primitive & lightweight progress bar, useful when terminal has no proper support for escape sequences.

// _________ IMPLEMENTATION _________

inline std::ostream* _output_stream = &std::cout;

inline void set_ostream(std::ostream& new_ostream) { _output_stream = &new_ostream; }

class Percentage {
private:
    char done_char;
    char not_done_char;
    bool show_time_estimate;

    std::size_t length_total;   // full   bar length
    std::size_t length_current; // filled bar length

    double last_update_percentage;
    double update_rate;

    using Clock     = std::chrono::steady_clock;
    using TimePoint = std::chrono::time_point<Clock>;

    TimePoint timepoint_start; // used to estimate remaining time
    TimePoint timepoint_current;

    int previous_string_length; // used to properly return the carriage when dealing with changed string size

    void draw_progressbar(double percentage) {
        const auto displayed_percentage = this->update_rate * std::floor(percentage / this->update_rate);
        // floor percentage to a closest multiple of 'update_rate' for nicer display
        // since actual updates are only required to happen no more ofter than 'update_rate'
        // and don't have to correspond to exact multiples of it

        // Estimate remaining time (linearly) and format it min + sec
        const auto time_elapsed       = this->timepoint_current - this->timepoint_start;
        const auto estimate_full      = time_elapsed / percentage;
        const auto estimate_remaining = estimate_full - time_elapsed;

        const auto estimate_remaining_sec = std::chrono::duration_cast<std::chrono::seconds>(estimate_remaining);

        const auto displayed_min = (estimate_remaining_sec / 60ll).count();
        const auto displayed_sec = (estimate_remaining_sec % 60ll).count();

        const bool show_min  = (displayed_min != 0);
        const bool show_sec  = (displayed_sec != 0) && !show_min;
        const bool show_time = (estimate_remaining_sec.count() > 0);

        std::stringstream ss;

        // Print bar
        ss << "[";
        std::fill_n(std::ostreambuf_iterator<char>(ss), this->length_current, this->done_char);
        std::fill_n(std::ostreambuf_iterator<char>(ss), this->length_total - this->length_current, this->not_done_char);
        ss << "]";

        // Print percentage
        ss << " " << std::fixed << std::setprecision(2) << 100. * displayed_percentage << "%";

        // Print time estimate
        if (this->show_time_estimate && show_time) {
            ss << " (remaining:";
            if (show_min) ss << " " << displayed_min << " min";
            if (show_sec) ss << " " << displayed_sec << " sec";
            ss << ")";
        }

        const std::string bar_string = ss.str();

        // Add spaces at the end to overwrite the previous string if it was longer that current
        const int current_string_length = static_cast<int>(bar_string.length());
        const int string_length_diff    = this->previous_string_length - current_string_length;

        if (string_length_diff > 0) { std::fill_n(std::ostreambuf_iterator<char>(ss), string_length_diff, ' '); }

        this->previous_string_length = current_string_length;

        // Return the carriage
        (*_output_stream) << ss.str(); // don't reuse 'bar_string', now 'ss' can also contain spaces at the end
        (*_output_stream) << '\r';
        (*_output_stream).flush();
        // '\r' returns cursor to the beginning of the line => most sensible consoles will
        // render render new lines over the last one. Otherwise every update produces a
        // bar on a new line, which looks worse but isn't critical for the purpose.
    }

public:
    Percentage(char done_char = '#', char not_done_char = '.', std::size_t bar_length = 30, double update_rate = 1e-2,
               bool show_time_estimate = true)
        : done_char(done_char), not_done_char(not_done_char), show_time_estimate(show_time_estimate),
          length_total(bar_length), length_current(0), last_update_percentage(0), update_rate(update_rate),
          timepoint_start(Clock::now()), previous_string_length(static_cast<int>(bar_length) + sizeof("[] 100.00%")) {}

    void start() {
        this->last_update_percentage = 0.;
        this->length_current         = 0;
        this->timepoint_start        = Clock::now();
        (*_output_stream) << "\n";
    }

    void set_progress(double percentage) {
        if (percentage - this->last_update_percentage <= this->update_rate) return;

        this->last_update_percentage = percentage;
        this->length_current         = static_cast<std::size_t>(percentage * static_cast<double>(this->length_total));
        this->timepoint_current      = Clock::now();
        this->draw_progressbar(percentage);
    }

    void finish() {
        this->last_update_percentage = 1.;
        this->length_current         = this->length_total;
        this->draw_progressbar(1.);
        (*_output_stream) << "\n";
    }
};

class Ruler {
private:
    char done_char;

    std::size_t length_total;
    std::size_t length_current;

public:
    Ruler(char done_char = '#') : done_char(done_char), length_total(51), length_current(0) {}

    void start() {
        this->length_current = 0;

        (*_output_stream) << "\n"
                          << " 0    10   20   30   40   50   60   70   80   90   100%\n"
                          << " |----|----|----|----|----|----|----|----|----|----|\n"
                          << " ";
    }

    void set_progress(double percentage) {
        const std::size_t length_new = static_cast<std::size_t>(percentage * static_cast<double>(this->length_total));

        if (length_new > length_current) {
            const auto chars_to_add = length_new - this->length_current;
            std::fill_n(std::ostreambuf_iterator<char>(*_output_stream), chars_to_add, this->done_char);
        }

        this->length_current = length_new;
    }

    void finish() {
        if (this->length_total > this->length_current) {
            const auto chars_to_add = this->length_total - this->length_current;
            std::fill_n(std::ostreambuf_iterator<char>(*_output_stream), chars_to_add, this->done_char);
        }

        this->length_current = this->length_total;

        (*_output_stream) << "\n";
    }
};

} // namespace utl::progressbar

// ========= header guard end =========

#endif
#endif






// -------------------------------
// --------- utl::random ---------
// -------------------------------
#if !defined(UTL_PICK_MODULES) || defined(UTLMODULE_RANDOM)
#ifndef UTLHEADERGUARD_RANDOM
#define UTLHEADERGUARD_RANDOM

// ======== header guard start ========

#include <ctime>
#include <initializer_list>
#include <limits>
#include <random>

namespace utl::random {

// _________ DEVELOPER DOCS _________

// Various convenient random functions, utilizes rand() internally.
// If good quality random is necessary use std::random generators instead.
//
// # XorShift64StarGenerator #
// Random 'std::uint64_t' generator that satisfies uniform random number generator requirements
// from '<random>' (see https://en.cppreference.com/w/cpp/named_req/UniformRandomBitGenerator).
// Implementation "XorShift64* suggested by Marsaglia G. in 2003 "Journal of Statistical Software"
// (see https://www.jstatsoft.org/article/view/v008i14).
// Slightly faster than 'std::rand()', while providing higher quality random.
// Significantly faster than 'std::mt19937'.
// State consists of a single 'std::uint64_t', requires seed >= 1.
//
// # xorshift64star #
// Global instance if XorShift64StarGenerator.
//
// # ::seed(), ::seed_with_time(), ::seed_with_random_device() #
// Seeds random with value/current_time/random_device.
// Random device is a better source of entropy, however it is more expensive to initialize
// than just taking current time with <ctime>.
//
// # ::rand_int(), ::rand_uint(), ::rand_float(), ::rand_double() #
// Random value in [min, max] range.
// Floats with no parameters assume range [0, 1].
//
// # ::rand_bool() #
// Randomly chooses 0 or 1.
//
// # ::rand_choise() #
// Randomly chooses a value from initializer list.
//
// # ::rand_linear_combination() #
// Produces "c A + (1-c) B" with random "0 < c < 1" assuming objects "A", "B" support arithmetic operations.
// Useful for vector and color operations.

// _________ IMPLEMENTATION _________

class XorShift64StarGenerator {
    // meets uniform random number generator requirements
    // (see https://en.cppreference.com/w/cpp/named_req/UniformRandomBitGenerator)
    //
    // implementation XorShift64* suggested by Marsaglia G. in 2003 "Journal of Statistical Software"
    // (see https://www.jstatsoft.org/article/view/v008i14)
public:
    using result_type = std::uint64_t;

private:
    result_type state;

public:
    XorShift64StarGenerator(result_type seed = 0) : state(seed + 1) {}

    static constexpr result_type min() { return std::numeric_limits<result_type>::min(); }
    static constexpr result_type max() { return std::numeric_limits<result_type>::max(); }

    void seed(result_type seed) { this->state = seed + 1; } // enforce seed >= 1 by always adding +1 to uint

    result_type operator()() {
        this->state ^= this->state >> 12;
        this->state ^= this->state << 25;
        this->state ^= this->state >> 27;
        return this->state * 0x2545F4914F6CDD1DULL;
    }
};

inline XorShift64StarGenerator xorshift64star;

inline void seed(XorShift64StarGenerator::result_type random_seed) { xorshift64star.seed(random_seed); }

inline void seed_with_time() { utl::random::seed(static_cast<XorShift64StarGenerator::result_type>(std::time(NULL))); }

inline void seed_with_random_device() {
    std::random_device rd;
    utl::random::seed(static_cast<XorShift64StarGenerator::result_type>(rd()));
}

// --- Uniform distribution shortcuts ---
// Note: despite common sense benchmarks show that creating 'std::uniform_..._distribution<>'
// repeatedly doesn't introduce any noticeble overhead
// sizeof(std::uniform_int_distribution<int>) == 8
// sizeof(std::uniform_real_distribution<double>) == 16

inline int rand_int(int min, int max) {
    std::uniform_int_distribution<int> distr{min, max};
    return distr(xorshift64star);
}

inline int rand_uint(unsigned int min, unsigned int max) {
    std::uniform_int_distribution<unsigned int> distr{min, max};
    return distr(xorshift64star);
}

inline float rand_float() { return std::generate_canonical<float, std::numeric_limits<float>::digits>(xorshift64star); }

inline float rand_float(float min, float max) {
    std::uniform_real_distribution<float> distr{min, max};
    return distr(xorshift64star);
}

inline double rand_double() {
    return std::generate_canonical<double, std::numeric_limits<double>::digits>(xorshift64star);
}

inline double rand_double(double min, double max) {
    std::uniform_real_distribution<double> distr{min, max};
    return distr(xorshift64star);
}

inline bool rand_bool() { return static_cast<bool>(xorshift64star() % 2); }

template <class T>
const T& rand_choise(std::initializer_list<T> objects) {
    const int random_index = rand_int(0, static_cast<int>(objects.size()) - 1);
    return objects.begin()[random_index];
}

template <class T>
T rand_linear_combination(const T& A, const T& B) { // random linear combination of 2 colors/vectors/etc
    const auto coef = rand_double();
    return A * coef + B * (1. - coef);
}

} // namespace utl::random

// ========= header guard end =========

#endif
#endif






// ------------------------------
// --------- utl::shell ---------
// ------------------------------
#if !defined(UTL_PICK_MODULES) || defined(UTLMODULE_SHELL)
#ifndef UTLHEADERGUARD_SHELL
#define UTLHEADERGUARD_SHELL

// ======== header guard start ========

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <ostream>
#include <sstream>
#include <string>
#include <string_view>
#include <unordered_set>
#include <vector>

namespace utl::shell {

// _________ DEVELOPER DOCS _________

// Command line utils that allow simple creation of temporary files and command line
// calls with stdout and stderr piping (a task surprisingly untrivial in standard C++).
//
// # ::random_ascii_string() #
// Creates random ASCII string of given length.
// Uses chars in ['a', 'z'] range.
//
// # ::generate_temp_file() #
// Generates temporary .txt file with a random unique name, and returns it's filepath.
// Files generated during current runtime can be deleted with '::clear_temp_files()'.
// If '::clear_temp_files()' wasn't called manually, it gets called automatically upon exiting 'main()'.
// Uses relative path internally.
//
// # ::clear_temp_files() #
// Clears temporary files generated during current runtime.
//
// # ::erase_temp_file() #
// Clears a single temporary file with given filepath.
//
// # ::run_command() #
// Runs a command using the default system shell.
// Returns piped status (error code), stdout and stderr.
//
// # ::exe_path() #
// Parses executable path from argcv as std::string_view.
//
// # ::command_line_args() #
// Parses command line arguments from argcv as std::string_view.

// _________ IMPLEMENTATION _________

// --- Random ASCII strings ---
inline std::string random_ascii_string(std::size_t length) {
    constexpr char min_char = 'a';
    constexpr char max_char = 'z';

    std::string result(length, '0');
    for (std::size_t i = 0; i < length; ++i)
        result[i] = static_cast<char>(min_char + rand() % (max_char - min_char + 1));
    return result;
}


// --- Temp files ---
inline std::unordered_set<std::string> _temp_files; // currently existing temp files
inline bool                            _temp_files_cleanup_registered = false;

inline void clear_temp_files() {
    for (const auto& file : _temp_files) std::filesystem::remove(file);
    _temp_files.clear();
}

inline void erase_temp_file(const std::string& file) {
    // we take 'file' as 'std::string&' instead of 'std::string_view' because it is
    // used to call '.erase()' on the map of 'std::string', which does not take string_view
    std::filesystem::remove(file);
    _temp_files.erase(file);
}

inline std::string generate_temp_file() {
    constexpr std::size_t MAX_ATTEMPTS = 500; // shouldn't realistically be encountered but still
    constexpr std::size_t NAME_LENGTH  = 30;

    // Register std::atexit() if not already registered
    if (!_temp_files_cleanup_registered) {
        const bool success             = (std::atexit(clear_temp_files) == 0);
        _temp_files_cleanup_registered = success;
    }

    // Try creating files until unique name is found
    for (std::size_t i = 0; i < MAX_ATTEMPTS; ++i) {
        const std::filesystem::path temp_path(random_ascii_string(NAME_LENGTH) + ".txt");

        if (std::filesystem::exists(temp_path)) continue;

        const std::ofstream temp_file(temp_path);

        if (temp_file.good()) {
            _temp_files.insert(temp_path.string());
            return temp_path.string();
        } else {
            return std::string();
        }
    }

    return std::string();
}


// --- Command line operations ---
struct CommandResult {
    int         status; // error code
    std::string stdout_output;
    std::string stderr_output;
};

inline CommandResult run_command(const std::string& command) {
    // we take 'std::string&' instead of 'std::string_view' because there
    // has to be a guarantee that contained string is null-terminated

    const auto stdout_file = utl::shell::generate_temp_file();
    const auto stderr_file = utl::shell::generate_temp_file();

    // Redirect stdout and stderr of the command to temporary files
    std::ostringstream ss;
    ss << command.c_str() << " >" << stdout_file << " 2>" << stderr_file;
    const std::string modified_command = ss.str();

    // Call command
    const auto status = std::system(modified_command.c_str());

    // Read stdout and stderr from temp files and remove them
    std::ostringstream stdout_stream;
    std::ostringstream stderr_stream;
    stdout_stream << std::ifstream(stdout_file).rdbuf();
    stderr_stream << std::ifstream(stderr_file).rdbuf();
    utl::shell::erase_temp_file(stdout_file);
    utl::shell::erase_temp_file(stderr_file);

    // Return
    CommandResult result = {status, stdout_stream.str(), stderr_stream.str()};

    return result;
}


// --- Argc/Argv parsing ---
inline std::string_view get_exe_path(char** argv) {
    // argc == 1 is a reasonable assumption since the only way to achieve such launch
    // is to run executable through a null-execv, most command-line programs assume
    // such scenario to be either impossible or an error on user side
    return std::string_view(argv[0]);
}

inline std::vector<std::string_view> get_command_line_args(int argc, char** argv) {
    std::vector<std::string_view> arguments(argc - 1);
    for (size_t i = 0; i < arguments.size(); ++i) arguments.emplace_back(argv[i]);
    return arguments;
}

} // namespace utl::shell

// ========= header guard end =========

#endif
#endif






// ------------------------------
// --------- utl::sleep ---------
// ------------------------------
#if !defined(UTL_PICK_MODULES) || defined(UTLMODULE_SLEEP)
#ifndef UTLHEADERGUARD_SLEEP
#define UTLHEADERGUARD_SLEEP

// ======== header guard start ========

#include <chrono>
#include <cmath>
#include <thread>

namespace utl::sleep {

// _________ DEVELOPER DOCS _________

// Various implementation of sleep(), used for precise delays.
//
// # ::spinlock() #
// Best precision, uses CPU.
//
// # ::hybrid() #
// Recommended option, similar precision to spinlock with minimal CPU usage.
// Loops short system sleep while statistically estimating its error on the fly and once within error
// margin of the end time, finished with spinlock sleep (essentialy negating system sleep error).
//
// # ::system() #
// Worst precision, frees CPU.

// _________ IMPLEMENTATION _________

using _clock     = std::chrono::steady_clock;
using _chrono_ns = std::chrono::nanoseconds;

inline void spinlock(double ms) {
    const long long ns              = static_cast<int64_t>(ms * 1e6);
    const auto      start_timepoint = _clock::now();

    volatile int i = 0; // volatile prevents standard-compliant compilers from optimizing away the loop
    while (std::chrono::duration_cast<_chrono_ns>(_clock::now() - start_timepoint).count() < ns) { ++i; }
}

inline void hybrid(double ms) {
    static double  estimate = 5e-3; // initial sleep_for() error estimate
    static double  mean     = estimate;
    static double  m2       = 0;
    static int64_t count    = 1;

    // We treat sleep_for(1 ms) as a random variate "1 ms + random_value()"
    while (ms > estimate) {
        const auto start = _clock::now();
        std::this_thread::sleep_for(_chrono_ns(static_cast<int64_t>(1e6)));
        const auto end = _clock::now();

        const double observed = std::chrono::duration_cast<_chrono_ns>(end - start).count() / 1e6;
        ms -= observed;

        ++count;

        // Welford's algorithm for mean and unbiased variance estimation
        const double delta = observed - mean;
        mean += delta / static_cast<double>(count);
        m2 += delta * (observed - mean); // intermediate values 'm2' reduce numerical instability
        const double variance = std::sqrt(m2 / static_cast<double>(count - 1));

        estimate = mean + variance; // set estimate 1 standard deviation above the mean
        // can be adjusted to make estimate more or less pessimistic
    }

    utl::sleep::spinlock(ms);
}

inline void system(double ms) { std::this_thread::sleep_for(_chrono_ns(static_cast<int64_t>(ms * 1e6))); }

} // namespace utl::sleep

// ========= header guard end =========

#endif
#endif






// ----------------------------
// --------- utl::mvl ---------
// ----------------------------
#if !defined(UTL_PICK_MODULES) || defined(UTLMODULE_STORAGE)
#ifndef UTLHEADERGUARD_STORAGE
#define UTLHEADERGUARD_STORAGE

// ======== header guard start ========

#include <cstddef>
#include <functional>
#include <iterator>
#include <memory>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace utl::mvl {

// _________ DEVELOPER DOCS _________

// NOTE: Docs here (I dread the day I'll have to explain this)
//
// Here matrix/vector/view/etc code reuse is implemented though macro 'mixins', while preprocessor appoarch is
// obviously not ideal, the alternatives have proven themselves either unsuitable or unrealistically cumbersome.
// Below is a list of possible approaches:
//
// 1) No code reuse, all classes implemented manually - unrealistically cumbersome, huge code duplication
//
// 2) Regular OOP with virtual classes - having vtables in lightweigh containers is highly undesirable, we can
//    avoid this using CRTP however we run into issues with multiple inheritance (see below)
//
// 3) CRTP - we unavoidably run into multiple inheritance issues (diamond problem) with concepts like:
//       [mutable-matrix-like] : [const-matrix-like] & [mutable-iterable] where both parents derive from
//       [const-iterable]
//
//    While those issues are resolvable, their handling becomes more and more difficult as parent glasses grow larger.
//    A usual solution for this would be virtual inheritance, however it prevents us from calling derived methods inside
//    the base class due to "static downcast via virtual inheritance" rule, which makes 'static_cast<Derived*>(this)'
//    illegal. We can resolve all name collisions manually, however it leads to a code that is absolutely horrific and
//    easy to mess up.
//
//    There is also an issue of extreme boilerplate, since we want to propagate '_types<>' wrapper and the 'Final'
//    return type to all base classes we end up with extremely unwieldy template arguments that make compile errors and
//    refactoring a challenge.
//
// 4) Template-based mixins without CRTP - these allow to add functionality on top of a class, essentially it works like
// "inverted" CRTP
//    from the code logic POV. While being almost a perfect solution, it doesn't allow us to return objects of derived
//    type from base class methods, which is a requirement for the intended matrix API.
//
// 5) Macro-based mixins - insert functionality directly in a class, a simple way of doing exactly what we want.
//    "Mixed-in" methods are simply functions implemented in terms of certain "required" methods.
//    Doesn't lead to huge template chains and doesn't have any inheritance at all. The main downside is being
//    preprocessor based, which exposes some identifiers to a global namespace (even if we #undef them later).

// _________ IMPLEMENTATION _________

// - Type pack wrapper -
// Wrapper that allows passing standard set of member types down the CRTP chain
// (types have to be passed through template args one way or another. Trying to access
// 'derived_type::value_type' from the CRTP base class will fail, since in order for it to be available
// the derived class has to be "implemented" which can only happen after the base class is "implemented",
// which requires the type => a logical loop. So we wrap all the types in a class a pass them down
// the chain a "pack with everything" type of deal. This pack then gets "unwrapped" with
// '_utl_storage_define_types' macro in every class => we have all the member types defined).
//
// Since CRTP approach has been deprecated in favor of mixins, original reasoning is no longer accurate,
// however it is this a useful class whenever a pack of types have to be passed through a template
// (like in iterator implementation)
template <class T>
struct _types {
    using value_type      = T;
    using size_type       = std::size_t;
    using difference_type = std::ptrdiff_t;
    using reference       = T&;
    using const_reference = const T&;
    using pointer         = T*;
    using const_pointer   = const T*;
};


// - Observer Ptr -

// Utility used in views to allow observer pointers with the same interface as std::unique_ptr,
// which means a generic '.data()` can be implemented without having to make a separate version for views and containers
// NOTE: A clutch really, but it works.
template <class T>
class _observer_ptr {
public:
    using element_type = T;

private:
    element_type* _data = nullptr;

public:
    // - Constructors -
    constexpr _observer_ptr() noexcept = default;
    constexpr _observer_ptr(std::nullptr_t) noexcept {}
    explicit _observer_ptr(element_type* ptr) : _data(ptr) {}

    // - Interface -
    constexpr element_type* get() const noexcept { return _data; }
};


// - Enums -

// Parameter enums
// Their combination + value type fully defines the type of matrix/matri-view
enum class Dimension { VECTOR, MATRIX };
enum class Type { DENSE, STRIDED, SPARSE };
enum class Ownership { CONTAINER, VIEW, CONST_VIEW };
enum class Checking { NONE, BOUNDS };
enum class Layout { /* 1D */ FLAT, /* 2D */ RC, CR, /* Other */ SPARSE };

// Overload tags (dummy types used to create "overloads" of .for_each(), which changes the way compiler handles
// name based lookup, preventing shadowing of base class methods)
// NOTE: Might be unnecessary after a switch to mixins
enum class _for_each_tag { DUMMY };
enum class _for_each_idx_tag { DUMMY };
enum class _for_each_ij_tag { DUMMY };


// - Iterator template -
// (reduces code duplication for const & non-const variants)
// Uses 'operator[]' of 'ParentPointer' for all logic, thus allowing arbitrary containers that
// don't have to be contiguous or even ordered in memory
template <class ParentPointer, class Types, bool is_const_iter = false>
class _flat_iterator {
    using parent_pointer = ParentPointer;

public:
    // Iterator reqs: [General]
    // Contains member types      =>
    using iterator_category = std::random_access_iterator_tag;
    using difference_type   = typename Types::difference_type;
    using value_type        = typename Types::value_type;
    using pointer = typename std::conditional_t<is_const_iter, typename Types::const_pointer, typename Types::pointer>;
    using reference =
        typename std::conditional_t<is_const_iter, typename Types::const_reference, typename Types::reference>;

    // Iterator reqs: [General]
    // Constructible
    _flat_iterator(parent_pointer parent, difference_type idx) : _parent(parent), _idx(idx) {}
    // Copy-constructible => by default
    // Copy-assignable    => by default
    // Destructible       => by default
    // Swappable
    friend void swap(_flat_iterator& lhs, _flat_iterator& rhs) { std::swap(lhs._idx, rhs._idx); }

    // Iterator reqs: [General]
    // Can be incremented (prefix & postfix)
    _flat_iterator& operator++() {
        ++_idx;
        return *this;
    } // prefix
    _flat_iterator operator++(int) {
        _flat_iterator temp = *this;
        ++_idx;
        return temp;
    } // postfix

    // Iterator reqs: [Input iterator]
    // Supports equality/inequality comparisons
    friend bool operator==(const _flat_iterator& it1, const _flat_iterator& it2) { return it1._idx == it2._idx; };
    friend bool operator!=(const _flat_iterator& it1, const _flat_iterator& it2) { return it1._idx != it2._idx; };
    // Can be dereferenced as an rvalue
    reference   operator*() const { return _parent->operator[](_idx); }
    pointer     operator->() const { return &_parent->operator[](_idx); }

    // Iterator reqs: [Output iterator]
    // Can be dereferenced as an lvalue (only for mutable iterator types) => not needed

    // Iterator reqs: [Forward iterator]
    // Default-constructible
    _flat_iterator() : _parent(nullptr), _idx(0) {}
    // "Multi-pass" - dereferencing & incrementing does not affects dereferenceability => satisfied

    // Iterator reqs: [Bidirectional iterator]
    // Can be decremented
    _flat_iterator& operator--() {
        --_idx;
        return *this;
    } // prefix
    _flat_iterator operator--(int) {
        _flat_iterator temp = *this;
        --_idx;
        return temp;
    } // postfix

    // Iterator reqs: [Random Access iterator]
    // See: https://en.cppreference.com/w/cpp/named_req/RandomAccessIterator
    // Supports arithmetic operators: it + n, n + it, it - n, it1 - it2
    friend _flat_iterator operator+(const _flat_iterator& it, difference_type diff) {
        return _flat_iterator(it._parent, it._idx + diff);
    }
    friend _flat_iterator operator+(difference_type diff, const _flat_iterator& it) {
        return _flat_iterator(it._parent, it._idx + diff);
    }
    friend _flat_iterator operator-(const _flat_iterator& it, difference_type diff) {
        return _flat_iterator(it._parent, it._idx - diff);
    }
    friend difference_type operator-(const _flat_iterator& it1, const _flat_iterator& it2) {
        return it1._idx - it2._idx;
    }
    // Supports inequality comparisons (<, >, <= and >=) between iterators
    friend bool     operator<(const _flat_iterator& it1, const _flat_iterator& it2) { return it1._idx < it2._idx; }
    friend bool     operator>(const _flat_iterator& it1, const _flat_iterator& it2) { return it1._idx > it2._idx; }
    friend bool     operator<=(const _flat_iterator& it1, const _flat_iterator& it2) { return it1._idx <= it2._idx; }
    friend bool     operator>=(const _flat_iterator& it1, const _flat_iterator& it2) { return it1._idx >= it2._idx; }
    // Standard assumption: Both iterators are from the same container => no need to compare '_parent'
    // Supports compound assignment operations += and -=
    _flat_iterator& operator+=(difference_type diff) {
        _idx += diff;
        return *this;
    }
    _flat_iterator& operator-=(difference_type diff) {
        _idx -= diff;
        return *this;
    }
    // Supports offset dereference operator ([])
    reference operator[](difference_type diff) const { return _parent->operator[](_idx + diff); }

private:
    parent_pointer  _parent;
    difference_type _idx; // not size_type because we have to use it in 'difference_type' operations most of the time
};


// - Helpers -

// Utility for creating unique_ptr arrays
template <class T>
std::unique_ptr<T[]> make_unique_ptr_array(size_t size) {
    return std::unique_ptr<T[]>(new T[size]);
}

// Shortuct for labda-type-based SFINAE
template <class FuncType, class Signature>
using _enable_if_signature = std::enable_if_t<std::is_convertible_v<FuncType, std::function<Signature>>, bool>;

// Variadic stringification
template <typename... Args>
std::string _stringify(const Args&... args) {
    std::stringstream ss;
    (ss << ... << args);
    return ss.str();
}

// Marker for uncreachable code
[[noreturn]] inline void _unreachable() {
// (Implementation from https://en.cppreference.com/w/cpp/utility/unreachable)
// Use compiler specific extensions if possible.
// Even if no extension is used, undefined behavior is still raised by
// an empty function body and the noreturn attribute.
#if defined(_MSC_VER) && !defined(__clang__) // MSVC
    __assume(false);
#else // GCC, Clang
    __builtin_unreachable();
#endif
}

// Pairs & Triplets for sparse matrices
template <typename T>
struct SparseEntry1D {
    size_t i;
    T      value;
};

template <typename T>
struct SparseEntry2D {
    size_t i;
    size_t j;
    T      value;
};

struct Index2D {
    size_t i;
    size_t j;
};

// Shortcut template used to deduce type of '_data' based on ownership inside mixins
template <typename MatrixType, typename ContainerResult, typename ViewResult, typename ConstViewResult>
using _choose_based_on_ownership = std::conditional_t<
    MatrixType::params::ownership == Ownership::CONTAINER, ContainerResult,
    std::conditional_t<MatrixType::params::ownership == Ownership::VIEW, ViewResult, ConstViewResult>>;


// - Type traits -

#define _utl_define_operator_support_type_trait(trait_name_, operator_)                                                \
    template <typename Type, typename = void>                                                                          \
    struct trait_name_ : std::false_type {};                                                                           \
                                                                                                                       \
    template <typename Type>                                                                                           \
    struct trait_name_<Type, std::void_t<decltype(std::declval<Type>() operator_ std::declval<Type>())>>               \
        : std::true_type {};                                                                                           \
                                                                                                                       \
    static_assert(true)

_utl_define_operator_support_type_trait(_supports_addition, +);
_utl_define_operator_support_type_trait(_supports_multiplication, *);
_utl_define_operator_support_type_trait(_supports_comparison, <);
_utl_define_operator_support_type_trait(_supports_eq_comparison, ==);

// Shortcut for enable_if'ing based on value_type traits
// Note: 'value_type = value_type' is crucial here, it allows function to be rejected without rejecting the whole class
#define _utl_require(condition_)                                                                                       \
    template <typename value_type = value_type, typename std::enable_if_t<condition_, bool> = true>

// - Template macros -

#define _utl_template_full_args                                                                                        \
    template <class T, Dimension _dimension, Type _type, Ownership _ownership, Checking _checking, Layout _layout>

#define _utl_template_param_args                                                                                       \
    template <Dimension _dimension, Type _type, Ownership _ownership, Checking _checking, Layout _layout>

// NOTE: This trick with inheritance allows us to forward template arguments to member variables before entering
// class definition, without this, 'params' struct would have no automatic way of getting 'dimension_', 'type_'
// and 'ownership_' of template specialization
#define _utl_template_specialization(dimension_, type_, ownership_)                                                    \
    template <class T, Checking _checking, Layout _layout>                                                             \
    class GenericTensor<T, dimension_, type_, ownership_, _checking, _layout>                                          \
        : public _params_wrapper<dimension_, type_, ownership_, _checking, _layout>


// - Mixins -

// - Types -

#define _utl_mixin_types                                                                                               \
private:                                                                                                               \
    using type_wrapper = _types<T>;                                                                                    \
                                                                                                                       \
public:                                                                                                                \
    using self = GenericTensor;                                                                                        \
                                                                                                                       \
    using value_type      = typename type_wrapper::value_type;                                                         \
    using size_type       = typename type_wrapper::size_type;                                                          \
    using difference_type = typename type_wrapper::difference_type;                                                    \
    using reference       = typename type_wrapper::reference;                                                          \
    using const_reference = typename type_wrapper::const_reference;                                                    \
    using pointer         = typename type_wrapper::pointer;                                                            \
    using const_pointer   = typename type_wrapper::const_pointer;                                                      \
                                                                                                                       \
private:                                                                                                               \
    static_assert(true)


// - Members -

// Requirements:
//    > NONE
#define _utl_mixin_extents_2d                                                                                          \
private:                                                                                                               \
    size_type _rows = 0;                                                                                               \
    size_type _cols = 0;                                                                                               \
                                                                                                                       \
public:                                                                                                                \
    size_type rows() const noexcept { return this->_rows; }                                                            \
    size_type cols() const noexcept { return this->_cols; }                                                            \
                                                                                                                       \
private:                                                                                                               \
    static_assert(true)


// Requirements:
//    > NONE
#define _utl_mixin_strides_2d                                                                                          \
private:                                                                                                               \
    size_type _row_stride;                                                                                             \
    size_type _col_stride;                                                                                             \
                                                                                                                       \
public:                                                                                                                \
    size_type row_stride() const noexcept { return this->_row_stride; }                                                \
    size_type col_stride() const noexcept { return this->_col_stride; }                                                \
                                                                                                                       \
private:                                                                                                               \
    static_assert(true)


// Requirements:
//    > NONE
#define _utl_mixin_data_2d_dense                                                                                       \
private:                                                                                                               \
    using _data_t = _choose_based_on_ownership<self, std::unique_ptr<value_type[]>, _observer_ptr<value_type>,         \
                                               _observer_ptr<const value_type>>;                                       \
                                                                                                                       \
public:                                                                                                                \
    _data_t _data;                                                                                                     \
                                                                                                                       \
private:                                                                                                               \
    static_assert(true)


// Requirements:
//    > NONE
#define _utl_mixin_data_2d_sparse                                                                                      \
private:                                                                                                               \
    using _data_t = _choose_based_on_ownership<self, std::vector<SparseEntry2D<value_type>>,                           \
                                               std::vector<SparseEntry2D<std::reference_wrapper<value_type>>>,         \
                                               std::vector<SparseEntry2D<std::reference_wrapper<const value_type>>>>;  \
                                                                                                                       \
public:                                                                                                                \
    _data_t _data;                                                                                                     \
                                                                                                                       \
private:                                                                                                               \
    static_assert(true)


/// - Methods -

// Requirements:
//    > const_reference operator()(size_type i, size_type j) const;
//    > size_type rows() const;
//    > size_type cols() const;
//    > Checking checking; // member or template arg
//    > Layout layout;     // member or template arg
#define _utl_mixin_conversions_2d_dense                                                                                \
private:                                                                                                               \
    size_type _unchecked_get_idx_of_ij(size_type i, size_type j) const {                                               \
        if constexpr (self::params::layout == Layout::RC) return i * this->cols() + j;                                 \
        if constexpr (self::params::layout == Layout::CR) return j * this->rows() + i;                                 \
        _unreachable();                                                                                                \
    }                                                                                                                  \
                                                                                                                       \
    Index2D _unchecked_get_ij_of_idx(size_type idx) const {                                                            \
        if constexpr (self::params::layout == Layout::RC) return {idx / this->cols(), idx % this->cols()};             \
        if constexpr (self::params::layout == Layout::CR) return {idx % this->rows(), idx / this->rows()};             \
        _unreachable();                                                                                                \
    }                                                                                                                  \
                                                                                                                       \
public:                                                                                                                \
    size_type get_idx_of_ij(size_type i, size_type j) const {                                                          \
        if constexpr (self::params::checking == Checking::BOUNDS) this->_bound_check_ij(i, j);                         \
        return _unchecked_get_idx_of_ij(i, j);                                                                         \
    }                                                                                                                  \
                                                                                                                       \
    Index2D get_ij_of_idx(size_type idx) const {                                                                       \
        if constexpr (self::params::checking == Checking::BOUNDS) this->_bound_check_idx(idx);                         \
        return _unchecked_get_ij_of_idx(idx);                                                                          \
    }                                                                                                                  \
                                                                                                                       \
    size_type extent_major() const noexcept {                                                                          \
        if constexpr (self::params::layout == Layout::RC) return this->rows();                                         \
        if constexpr (self::params::layout == Layout::CR) return this->cols();                                         \
        _unreachable();                                                                                                \
    }                                                                                                                  \
                                                                                                                       \
    size_type extent_minor() const noexcept {                                                                          \
        if constexpr (self::params::layout == Layout::RC) return this->cols();                                         \
        if constexpr (self::params::layout == Layout::CR) return this->rows();                                         \
        _unreachable();                                                                                                \
    }                                                                                                                  \
                                                                                                                       \
private:                                                                                                               \
    static_assert(true)


// Requirements:
//    > _data.get()                                         ~ provided by '_utl_mixin_data_2d_dense'
//    > rows()                                              ~ provided by '_utl_mixin_extents_2d'
//    > cols()                                              ~ provided by '_utl_mixin_extents_2d'
//    > _bound_check_idx(size_type i, size_type j)          ~ provided by '_utl_mixin_indexable_2d_const'
//    > get_idx_of_ij(i, j)                                 ~ provided by '_utl_mixin_conversions_2d_dense'
//    > self::params::checking                              ~ provided by base class
// Note: Bound checking in operator[] & operator() happens during index conversion when computing memory offset
#define _utl_mixin_indexation_2d_dense_const                                                                           \
public:                                                                                                                \
    size_type get_memory_offset_of_idx(size_type idx) const {                                                          \
        if constexpr (self::params::checking == Checking::BOUNDS) this->_bound_check_idx(idx);                         \
        return idx;                                                                                                    \
    }                                                                                                                  \
                                                                                                                       \
    size_type get_memory_offset_of_ij(size_type i, size_type j) const { return this->get_idx_of_ij(i, j); }            \
                                                                                                                       \
public:                                                                                                                \
    const_pointer data() const noexcept { return this->_data.get(); }                                                  \
                                                                                                                       \
    size_type size() const noexcept { return this->rows() * this->cols(); }                                            \
                                                                                                                       \
    const_reference operator[](size_type idx) const { return this->data()[this->get_memory_offset_of_idx(idx)]; }      \
                                                                                                                       \
    const_reference operator()(size_type i, size_type j) const {                                                       \
        return this->data()[this->get_memory_offset_of_ij(i, j)];                                                      \
    }                                                                                                                  \
                                                                                                                       \
private:                                                                                                               \
    static_assert(true)


// Requirements:
//    > _data.get()                                         ~ provided by '_utl_mixin_data_2d_dense'
//    > get_memory_offset_of_idx(idx)                       ~ provided by '_utl_mixin_indexation_2d_dense_const'
//    > this->get_memory_offset_of_ij(i, j)                 ~ provided by '_utl_mixin_indexation_2d_dense_const'
#define _utl_mixin_indexation_2d_dense_mutable                                                                         \
public:                                                                                                                \
    pointer data() noexcept { return this->_data.get(); }                                                              \
                                                                                                                       \
    reference operator[](size_type idx) { return this->data()[this->get_memory_offset_of_idx(idx)]; }                  \
                                                                                                                       \
    reference operator()(size_type i, size_type j) { return this->data()[this->get_memory_offset_of_ij(i, j)]; }       \
                                                                                                                       \
private:                                                                                                               \
    static_assert(true)


// Requirements:
//    > _data.get()                                         ~ provided by '_utl_mixin_data_2d_dense'
//    > rows()                                              ~ provided by '_utl_mixin_extents_2d'
//    > cols()                                              ~ provided by '_utl_mixin_extents_2d'
//    > row_stride()                                        ~ provided by '_utl_mixin_extents_2d'
//    > col_stride()                                        ~ provided by '_utl_mixin_extents_2d'
//    > get_idx_of_ij(i, j)                                 ~ provided by '_utl_mixin_conversions_2d_dense'
//    > get_ij_of_idx(idx)                                  ~ provided by '_utl_mixin_conversions_2d_dense'
//    > self::params::layout                                ~ provided by base class
// Note: Bound checking in operator[] & operator() happens during index conversion when computing memory offset
#define _utl_mixin_indexation_2d_strided_const                                                                         \
private:                                                                                                               \
    size_type _get_memory_offset_impl(size_type idx, size_type i, size_type j) const {                                 \
        if constexpr (self::params::layout == Layout::RC) return idx * this->col_stride() + this->row_stride() * i;    \
        if constexpr (self::params::layout == Layout::CR) return idx * this->row_stride() + this->col_stride() * j;    \
        _unreachable();                                                                                                \
    }                                                                                                                  \
                                                                                                                       \
public:                                                                                                                \
    size_type get_memory_offset_of_idx(size_type idx) const {                                                          \
        const auto ij = this->get_ij_of_idx(idx);                                                                      \
        return _get_memory_offset_impl(idx, ij.i, ij.j);                                                               \
    }                                                                                                                  \
                                                                                                                       \
    size_type get_memory_offset_of_ij(size_type i, size_type j) const {                                                \
        const auto idx = this->get_idx_of_ij(i, j);                                                                    \
        return _get_memory_offset_impl(idx, i, j);                                                                     \
    }                                                                                                                  \
                                                                                                                       \
public:                                                                                                                \
    const_pointer data() const noexcept { return this->_data.get(); }                                                  \
                                                                                                                       \
    size_type size() const noexcept { return this->rows() * this->cols(); }                                            \
                                                                                                                       \
    const_reference operator[](size_type idx) const { return this->data()[this->get_memory_offset_of_idx(idx)]; }      \
                                                                                                                       \
    const_reference operator()(size_type i, size_type j) const {                                                       \
        return this->data()[this->get_memory_offset_of_ij(i, j)];                                                      \
    }                                                                                                                  \
                                                                                                                       \
private:                                                                                                               \
    static_assert(true)


// Requirements:
//    > _data.get()                                         ~ provided by '_utl_mixin_data_2d_dense'
//    > get_memory_offset_of_idx(idx)                       ~ provided by '_utl_mixin_indexation_2d_dense_const'
//    > this->get_memory_offset_of_ij(i, j)                 ~ provided by '_utl_mixin_indexation_2d_dense_const'
#define _utl_mixin_indexation_2d_strided_mutable                                                                       \
public:                                                                                                                \
    pointer data() noexcept { return this->_data.get(); }                                                              \
                                                                                                                       \
    reference operator[](size_type idx) { return this->data()[this->get_memory_offset_of_idx(idx)]; }                  \
                                                                                                                       \
    reference operator()(size_type i, size_type j) { return this->data()[this->get_memory_offset_of_ij(i, j)]; }       \
                                                                                                                       \
private:                                                                                                               \
    static_assert(true)


// Requirements:
//    > const_reference operator[](size_type idx) const;
//    > size_type size() const;
#define _utl_mixin_indexable_1d_const                                                                                  \
private:                                                                                                               \
    void _bound_check_idx(size_type idx) const {                                                                       \
        if (idx >= this->size())                                                                                       \
            throw std::out_of_range(                                                                                   \
                _stringify("idx (which is ", idx, ") >= this->size() (which is ", this->size(), ")"));                 \
    }                                                                                                                  \
                                                                                                                       \
public:                                                                                                                \
    bool empty() const noexcept { return (this->size() == 0); }                                                        \
                                                                                                                       \
    const_reference front() const { return this->operator[](0); }                                                      \
    const_reference back() const { return this->operator[](this->size() - 1); }                                        \
                                                                                                                       \
    std::vector<value_type> to_std_vector() const { return std::vector(this->cbegin(), this->cend()); }                \
                                                                                                                       \
    template <class FuncType, _enable_if_signature<FuncType, void(const_reference)> = true>                            \
    const self& for_each(FuncType func, _for_each_tag = _for_each_tag::DUMMY) const {                                  \
        for (size_type idx = 0; idx < this->size(); ++idx) func(this->operator[](idx));                                \
        return *this;                                                                                                  \
    }                                                                                                                  \
                                                                                                                       \
    template <class FuncType, _enable_if_signature<FuncType, void(const_reference, size_type)> = true>                 \
    const self& for_each(FuncType func, _for_each_idx_tag = _for_each_idx_tag::DUMMY) const {                          \
        for (size_type idx = 0; idx < this->size(); ++idx) func(this->operator[](idx), idx);                           \
        return *this;                                                                                                  \
    }                                                                                                                  \
                                                                                                                       \
    template <class PredType, _enable_if_signature<PredType, bool(const_reference)> = true>                            \
    bool true_for_any(PredType predicate) const {                                                                      \
        for (size_type idx = 0; idx < this->size(); ++idx)                                                             \
            if (predicate(this->operator[](idx), idx)) return true;                                                    \
        return false;                                                                                                  \
    }                                                                                                                  \
                                                                                                                       \
    template <class PredType, _enable_if_signature<PredType, bool(const_reference, size_type)> = true>                 \
    bool true_for_any(PredType predicate) const {                                                                      \
        for (size_type idx = 0; idx < this->size(); ++idx)                                                             \
            if (predicate(this->operator[](idx), idx)) return true;                                                    \
        return false;                                                                                                  \
    }                                                                                                                  \
                                                                                                                       \
    template <class PredType, _enable_if_signature<PredType, bool(const_reference)> = true>                            \
    bool true_for_all(PredType predicate) const {                                                                      \
        auto inversed_predicate = [&](const_reference e) -> bool { return !predicate(e); };                            \
        return !this->true_for_any(inversed_predicate);                                                                \
    }                                                                                                                  \
                                                                                                                       \
    template <class PredType, _enable_if_signature<PredType, bool(const_reference, size_type)> = true>                 \
    bool true_for_all(PredType predicate) const {                                                                      \
        auto inversed_predicate = [&](const_reference e, size_type idx) -> bool { return !predicate(e, idx); };        \
        return !this->true_for_any(inversed_predicate);                                                                \
    }                                                                                                                  \
                                                                                                                       \
    _utl_require(_supports_addition<value_type>::value) value_type sum() const {                                       \
        return std::accumulate(this->cbegin(), this->cend(), value_type());                                            \
    }                                                                                                                  \
                                                                                                                       \
    _utl_require(_supports_multiplication<value_type>::value) value_type product() const {                             \
        return std::accumulate(this->cbegin(), this->cend(), value_type(), std::multiplies<value_type>());             \
    }                                                                                                                  \
                                                                                                                       \
    _utl_require(_supports_comparison<value_type>::value) value_type min() const {                                     \
        return *std::min_element(this->cbegin(), this->cend());                                                        \
    }                                                                                                                  \
                                                                                                                       \
    _utl_require(_supports_comparison<value_type>::value) value_type max() const {                                     \
        return *std::max_element(this->cbegin(), this->cend());                                                        \
    }                                                                                                                  \
                                                                                                                       \
    _utl_require(_supports_comparison<value_type>::value) bool is_sorted() const {                                     \
        return std::is_sorted(this->cbegin(), this->cend());                                                           \
    }                                                                                                                  \
                                                                                                                       \
    template <typename Compare>                                                                                        \
    bool is_sorted(Compare cmp) const {                                                                                \
        return std::is_sorted(this->cbegin(), this->cend(), cmp);                                                      \
    }                                                                                                                  \
                                                                                                                       \
    _utl_require(_supports_eq_comparison<value_type>::value) bool contains(const_reference value) const {              \
        return std::find(this->cbegin(), this->cend(), value) != this->cend();                                         \
    }                                                                                                                  \
                                                                                                                       \
    _utl_require(_supports_eq_comparison<value_type>::value) size_type count(const_reference value) const {            \
        return std::count(this->cbegin(), this->cend(), value);                                                        \
    }                                                                                                                  \
                                                                                                                       \
    using const_iterator         = _flat_iterator<const self*, type_wrapper, true>;                                    \
    using const_reverse_iterator = std::reverse_iterator<const_iterator>;                                              \
                                                                                                                       \
    const_iterator cbegin() const { return const_iterator(this, 0); }                                                  \
    const_iterator cend() const { return const_iterator(this, this->size()); }                                         \
                                                                                                                       \
    const_iterator begin() const { return this->cbegin(); }                                                            \
    const_iterator end() const { return this->cend(); }                                                                \
                                                                                                                       \
    const_reverse_iterator crbegin() const { return const_reverse_iterator(this->cend()); }                            \
    const_reverse_iterator crend() const { return const_reverse_iterator(this->cbegin()); }                            \
                                                                                                                       \
    const_reverse_iterator rbegin() const { return this->crbegin(); }                                                  \
    const_reverse_iterator rend() const { return this->crend(); }                                                      \
                                                                                                                       \
private:                                                                                                               \
    static_assert(true)


// Requirements:
//    > reference operator[](size_type idx);
//    > size_type size() const;
#define _utl_mixin_indexable_1d_mutable                                                                                \
public:                                                                                                                \
    reference front() { return this->operator[](0); }                                                                  \
    reference back() { return this->operator[](this->size() - 1); }                                                    \
                                                                                                                       \
    template <class FuncType, _enable_if_signature<FuncType, void(reference)> = true>                                  \
    self& for_each(FuncType func) {                                                                                    \
        for (size_type idx = 0; idx < this->size(); ++idx) func(this->operator[](idx));                                \
        return *this;                                                                                                  \
    }                                                                                                                  \
                                                                                                                       \
    template <class FuncType, _enable_if_signature<FuncType, void(reference, size_type)> = true>                       \
    self& for_each(FuncType func) {                                                                                    \
        for (size_type idx = 0; idx < this->size(); ++idx) func(this->operator[](idx), idx);                           \
        return *this;                                                                                                  \
    }                                                                                                                  \
                                                                                                                       \
    self& fill(const_reference value) {                                                                                \
        for (size_type idx = 0; idx < this->size(); ++idx) this->operator[](idx) = value;                              \
        return *this;                                                                                                  \
    }                                                                                                                  \
                                                                                                                       \
    _utl_require(_supports_comparison<value_type>::value) self& sort() {                                               \
        std::sort(this->begin(), this->end());                                                                         \
        return *this;                                                                                                  \
    }                                                                                                                  \
                                                                                                                       \
    template <typename Compare>                                                                                        \
    self& sort(Compare cmp) {                                                                                          \
        std::sort(this->begin(), this->end(), cmp);                                                                    \
        return *this;                                                                                                  \
    }                                                                                                                  \
    _utl_require(_supports_comparison<value_type>::value) self& stable_sort() {                                        \
        std::stable_sort(this->begin(), this->end());                                                                  \
        return *this;                                                                                                  \
    }                                                                                                                  \
                                                                                                                       \
    template <typename Compare>                                                                                        \
    self& stable_sort(Compare cmp) {                                                                                   \
        std::stable_sort(this->begin(), this->end(), cmp);                                                             \
        return *this;                                                                                                  \
    }                                                                                                                  \
                                                                                                                       \
    using iterator         = _flat_iterator<self*, type_wrapper>;                                                      \
    using reverse_iterator = std::reverse_iterator<iterator>;                                                          \
                                                                                                                       \
    iterator begin() { return iterator(this, 0); }                                                                     \
    iterator end() { return iterator(this, this->size()); }                                                            \
                                                                                                                       \
    reverse_iterator rbegin() { return reverse_iterator(this->end()); }                                                \
    reverse_iterator rend() { return reverse_iterator(this->begin()); }                                                \
                                                                                                                       \
private:                                                                                                               \
    static_assert(true)
// 1D iterators
// std::vector-like API and conversions
// 1D bound checking
// 1D .for_each() and algos


// Requirements:
//    > const_reference operator()(size_type i, size_type j) const;
//    > size_type rows() const;
//    > size_type cols() const;
#define _utl_mixin_indexable_2d_const                                                                                  \
private:                                                                                                               \
    void _bound_check_ij(size_type i, size_type j) const {                                                             \
        if (i >= this->rows())                                                                                         \
            throw std::out_of_range(_stringify("i (which is ", i, ") >= this->rows() (which is ", this->rows(), ")")); \
        else if (j >= this->cols())                                                                                    \
            throw std::out_of_range(_stringify("j (which is ", j, ") >= this->cols() (which is ", this->cols(), ")")); \
    }                                                                                                                  \
                                                                                                                       \
public:                                                                                                                \
    template <class FuncType, _enable_if_signature<FuncType, void(const_reference, size_type, size_type)> = true>      \
    const self& for_each(FuncType func) const {                                                                        \
        /* Loop over all 2D indices using 1D loop with idx->ij conversion */                                           \
        /* This is just as fast and ensures looping only over existing elements in non-dense matrices */               \
        for (size_type idx = 0; idx < this->size(); ++idx) {                                                           \
            const auto ij = this->get_ij_of_idx(idx);                                                                  \
            func(this->operator[](idx), ij.i, ij.j);                                                                   \
        }                                                                                                              \
        return *this;                                                                                                  \
    }                                                                                                                  \
                                                                                                                       \
    template <class PredType, _enable_if_signature<PredType, bool(const_reference, size_type, size_type)> = true>      \
    bool true_for_any(PredType predicate) const {                                                                      \
        /* Loop over all 2D indices using 1D loop with idx->ij conversion */                                           \
        /* This is just as fast and ensures looping only over existing elements in non-dense matrices */               \
        for (size_type idx = 0; idx < this->size(); ++idx) {                                                           \
            const auto ij = this->get_ij_of_idx(idx);                                                                  \
            if (predicate(this->operator[](idx), ij.i, ij.j)) return true;                                             \
        }                                                                                                              \
        return false;                                                                                                  \
    }                                                                                                                  \
                                                                                                                       \
    template <class PredType, _enable_if_signature<PredType, bool(const_reference, size_type, size_type)> = true>      \
    bool true_for_all(PredType predicate) const {                                                                      \
        /* We can reuse .true_for_any() with inverted predicate due to following conjecture: */                        \
        /* FOR_ALL (predicate)  ~  ! FOR_ANY (!predicate) */                                                           \
        auto inversed_predicate = [&](const_reference e, size_type i, size_type j) -> bool {                           \
            return !predicate(e, i, j);                                                                                \
        };                                                                                                             \
        return !this->true_for_any(inversed_predicate);                                                                \
    }                                                                                                                  \
                                                                                                                       \
    template <Type other_type, Ownership other_ownership, Checking other_checking, Layout other_layout>                \
    bool compare_contents(const GenericTensor<value_type, self::params::dimension, other_type, other_ownership,        \
                                              other_checking, other_layout>& other) const {                            \
        /* Surface-level checks */                                                                                     \
        if ((this->rows() != other.rows()) || (this->cols() != other.cols())) return false;                            \
        /* Compare while respecting sparsity */                                                                        \
        constexpr bool is_sparse_l = (self::params::type == Type::SPARSE);                                             \
        constexpr bool is_sparse_r = (std::remove_reference_t<decltype(other)>::params::type == Type::SPARSE);         \
        /* Same sparsity comparison */                                                                                 \
        if constexpr (is_sparse_l == is_sparse_r) {                                                                    \
            return this->size() == other.size() &&                                                                     \
                   this->true_for_all([&](const_reference e, size_type i, size_type j) { return e == other(i, j); });  \
        }                                                                                                              \
        /* Different sparsity comparison */                                                                            \
        /* TODO: Impl here and use .all_of() OR .any_of() */                                                           \
        return true;                                                                                                   \
    }                                                                                                                  \
                                                                                                                       \
private:                                                                                                               \
    static_assert(true)


// Requirements:
//    > const_reference operator()(size_type i, size_type j) const;
//    > size_type rows() const;
//    > size_type cols() const;
#define _utl_mixin_indexable_2d_mutable                                                                                \
public:                                                                                                                \
    template <class FuncType, _enable_if_signature<FuncType, void(reference, size_type, size_type)> = true>            \
    self& for_each(FuncType func, _for_each_ij_tag = _for_each_ij_tag::DUMMY) {                                        \
        for (size_type idx = 0; idx < this->size(); ++idx) {                                                           \
            const auto ij = this->get_ij_of_idx(idx);                                                                  \
            func(this->operator[](idx), ij.i, ij.j);                                                                   \
        }                                                                                                              \
        return *this;                                                                                                  \
    }                                                                                                                  \
                                                                                                                       \
private:                                                                                                               \
    static_assert(true)
// 1D <-> 2D index conversions
// 2D bound checking
// 2D .for_each() and algos


// Requirements:
//    > self& operator=(const self& other)
#define _utl_mixin_default_constructors_dense                                                                          \
public:                                                                                                                \
    GenericTensor(const self& other) { *this = other; }                                                                \
                                                                                                                       \
    self& operator=(self&& other) noexcept = default;                                                                  \
    GenericTensor(self&& other) noexcept   = default;                                                                  \
                                                                                                                       \
    template <Type other_type, Ownership other_ownership, Checking other_checking, Layout other_layout>                \
    GenericTensor(const GenericTensor<value_type, self::params::dimension, other_type, other_ownership,                \
                                      other_checking, other_layout>& other)                                            \
        : GenericTensor(other.rows(), other.cols()) {                                                                  \
        other.for_each(                                                                                                \
            [&](const value_type& element, size_type i, size_type j) { this->operator()(i, j) = element; });           \
    }                                                                                                                  \
                                                                                                                       \
private:                                                                                                               \
    static_assert(true)
//
// copy ctor induced from copy assignment.
//
// move ctor and assignment are trivial.
//
// copy-from-tensor ctor: Any matrix-like object with the same 'value_type' and dimension.
// This is not a copy constructor even if it logically implements it. Copy ctor needs to take 'const self&'.


// - Tensor Definition -

// Needed in order to forward template args to properties in specialization macro
_utl_template_param_args struct _params_wrapper {
    struct params {
        constexpr static auto dimension = _dimension;
        constexpr static auto type      = _type;
        constexpr static auto ownership = _ownership;
        constexpr static auto checking  = _checking;
        constexpr static auto layout    = _layout;

        // Prevent impossible layouts
        static_assert((dimension == Dimension::VECTOR) == (layout == Layout::FLAT), "Flat layout <=> matrix is 1D.");

        static_assert((type == Type::SPARSE) == (layout == Layout::SPARSE), "Sparse layout <=> matrix is sparse.");
    };
};

_utl_template_full_args class GenericTensor;


// - Partial Specializations -

// Matrix
_utl_template_specialization(Dimension::MATRIX, Type::DENSE, Ownership::CONTAINER) {
    _utl_mixin_types;

    _utl_mixin_extents_2d;
    _utl_mixin_data_2d_dense;

    _utl_mixin_indexation_2d_dense_const;
    _utl_mixin_indexation_2d_dense_mutable;
    _utl_mixin_indexable_1d_const;
    _utl_mixin_indexable_1d_mutable;
    _utl_mixin_indexable_2d_const;
    _utl_mixin_indexable_2d_mutable;
    _utl_mixin_conversions_2d_dense;
    _utl_mixin_default_constructors_dense;

public:
    // Default ctor
    GenericTensor() noexcept = default;

    // Copy-ctor
    self& operator=(const self& other) {
        this->_rows = other.rows();
        this->_cols = other.cols();
        this->_data = make_unique_ptr_array<value_type>(this->size());
        std::copy(other.begin(), other.end(), this->begin());
        return *this;
    }

    // Init-with-value
    explicit GenericTensor(size_type rows, size_type cols, const_reference value = value_type()) noexcept
        : _rows(rows), _cols(cols) {
        this->_data = make_unique_ptr_array<value_type>(this->size());
        this->fill(value);
    }

    // Init-with-data
    explicit GenericTensor(size_type rows, size_type cols, pointer data_ptr) noexcept
        : _data(data_ptr), _rows(rows), _cols(cols) {}

    // Init-with-ilist
    GenericTensor(std::initializer_list<std::initializer_list<value_type>> init) {
        this->_rows = init.size();
        this->_cols = (*init.begin()).size();
        this->_data = make_unique_ptr_array<value_type>(this->size());

        // Check dimensions (throw if cols have different dimensions)
        for (auto row_it = init.begin(); row_it < init.end(); ++row_it)
            if ((*row_it).end() - (*row_it).begin() != this->_cols)
                throw std::invalid_argument("Initializer list dimensions don't match.");

        // Copy elements
        for (size_type i = 0; i < this->rows(); ++i)
            for (size_type j = 0; j < this->cols(); ++j) this->operator()(i, j) = (init.begin()[i]).begin()[j];
    }

    // Move-from-other-matrix
    /* Accepts: Any matrix of the same partial specialization */
    template <Checking other_checking, Layout other_layout>
    GenericTensor(GenericTensor<value_type, self::params::dimension, self::params::type, self::params::ownership,
                                other_checking, other_layout> &&
                  other)
        : _data(std::move(other._data)), _rows(other.rows()), _cols(other.cols()) {}

    // NOTE: Move to mixins
    self transposed() const {
        self res(this->cols(), this->rows());
        this->for_each([&](const value_type& element, size_type i, size_type j) { res(j, i) = element; });
        return res;
    }
};

// MatrixView
_utl_template_specialization(Dimension::MATRIX, Type::DENSE, Ownership::VIEW) {
    _utl_mixin_types;

    _utl_mixin_extents_2d;
    _utl_mixin_data_2d_dense;

    _utl_mixin_indexation_2d_dense_const;
    _utl_mixin_indexation_2d_dense_mutable;
    _utl_mixin_indexable_1d_const;
    _utl_mixin_indexable_1d_mutable;
    _utl_mixin_indexable_2d_const;
    _utl_mixin_indexable_2d_mutable;
    _utl_mixin_conversions_2d_dense;

public:
    // Init-from-tensor
    /* Accepts: Any tensor of the same type */
    template <Ownership other_ownership, Checking other_checking, Layout other_layout>
    GenericTensor(GenericTensor<value_type, self::params::dimension, self::params::type, other_ownership,
                                other_checking, other_layout> &
                  other)
        : _data(other.data()), _rows(other.rows()), _cols(other.cols()) {}

    // Init-from-data
    GenericTensor(size_type rows, size_type cols, pointer data_ptr) : _data(data_ptr), _rows(rows), _cols(cols) {}
};

// ConstMatrixView
_utl_template_specialization(Dimension::MATRIX, Type::DENSE, Ownership::CONST_VIEW) {
    _utl_mixin_types;

    _utl_mixin_extents_2d;
    _utl_mixin_data_2d_dense;

    _utl_mixin_indexation_2d_dense_const;
    _utl_mixin_indexable_1d_const;
    _utl_mixin_indexable_2d_const;
    _utl_mixin_conversions_2d_dense;

public:
    // Init-from-tensor
    /* Accepts: Any tensor of the same type */
    template <Ownership other_ownership, Checking other_checking, Layout other_layout>
    GenericTensor(const GenericTensor<value_type, self::params::dimension, self::params::type, other_ownership,
                                      other_checking, other_layout>& other)
        : _data(other.data()), _rows(other.rows()), _cols(other.cols()) {}

    // Init-from-data
    GenericTensor(size_type rows, size_type cols, const_pointer data_ptr) : _data(data_ptr), _rows(rows), _cols(cols) {}

    // TEMP:
    std::string dump() const {
        std::stringstream ss;
        for (size_type i = 0; i < this->rows(); ++i)
            for (size_type j = 0; j < this->cols(); ++j)
                ss << "(" << i << ", " << j << ") ~ [" << this->get_idx_of_ij(i, j) << "] = " << this->operator()(i, j)
                   << "\n";
        return ss.str();
    }
    // TEMP:
};

// StridedMatrixView
_utl_template_specialization(Dimension::MATRIX, Type::STRIDED, Ownership::VIEW) {
    _utl_mixin_types;

    _utl_mixin_extents_2d;
    _utl_mixin_strides_2d;
    _utl_mixin_data_2d_dense;

    _utl_mixin_indexation_2d_strided_const;
    _utl_mixin_indexation_2d_strided_mutable;
    _utl_mixin_indexable_1d_const;
    _utl_mixin_indexable_1d_mutable;
    _utl_mixin_indexable_2d_const;
    _utl_mixin_indexable_2d_mutable;
    _utl_mixin_conversions_2d_dense;

public:
    // Init-from-tensor
    /* Accepts: Any tensor of the same type */
    template <Ownership other_ownership, Checking other_checking, Layout other_layout>
    GenericTensor(GenericTensor<value_type, self::params::dimension, self::params::type, other_ownership,
                                other_checking, other_layout> &
                  other)
        : _data(other.data()), _rows(other.rows()), _cols(other.cols()), _row_stride(other.row_stride()),
          _col_stride(other.col_stride()) {}

    // Init-from-data
    GenericTensor(size_type rows, size_type cols, size_type row_stride, size_type col_stride, pointer data_ptr)
        : _data(data_ptr), _rows(rows), _cols(cols), _row_stride(row_stride), _col_stride(col_stride) {}
};

// SparseMatrix
_utl_template_specialization(Dimension::MATRIX, Type::SPARSE, Ownership::CONTAINER) {
    _utl_mixin_types;

    _utl_mixin_extents_2d;
    _utl_mixin_data_2d_sparse;

    _utl_mixin_indexable_1d_const;
    _utl_mixin_indexable_1d_mutable;
    _utl_mixin_indexable_2d_const;
    _utl_mixin_indexable_2d_mutable;
    _utl_mixin_conversions_2d_dense;
};

// SparseMatrixView
_utl_template_specialization(Dimension::MATRIX, Type::SPARSE, Ownership::CONST_VIEW) {
    _utl_mixin_types;

    _utl_mixin_extents_2d;
    _utl_mixin_data_2d_sparse;

    _utl_mixin_indexable_1d_const;
    _utl_mixin_indexable_2d_const;
};


// - Implementation wrappers -

constexpr auto _default_checking        = Checking::NONE;
constexpr auto _default_layout_dense_2d = Layout::RC;

// Dense 2D
template <typename T, Checking checking = _default_checking, Layout layout = _default_layout_dense_2d>
using Matrix = GenericTensor<T, Dimension::MATRIX, Type::DENSE, Ownership::CONTAINER, checking, layout>;

template <typename T, Checking checking = _default_checking, Layout layout = _default_layout_dense_2d>
using MatrixView = GenericTensor<T, Dimension::MATRIX, Type::DENSE, Ownership::VIEW, checking, layout>;

template <typename T, Checking checking = _default_checking, Layout layout = _default_layout_dense_2d>
using ConstMatrixView = GenericTensor<T, Dimension::MATRIX, Type::DENSE, Ownership::CONST_VIEW, checking, layout>;

// Strided 2D
template <typename T, Checking checking = _default_checking, Layout layout = _default_layout_dense_2d>
using StridedMatrix = GenericTensor<T, Dimension::MATRIX, Type::STRIDED, Ownership::CONTAINER, checking, layout>;

template <typename T, Checking checking = _default_checking, Layout layout = _default_layout_dense_2d>
using StridedMatrixView = GenericTensor<T, Dimension::MATRIX, Type::STRIDED, Ownership::VIEW, checking, layout>;

template <typename T, Checking checking = _default_checking, Layout layout = _default_layout_dense_2d>
using ConstStridedMatrixView =
    GenericTensor<T, Dimension::MATRIX, Type::STRIDED, Ownership::CONST_VIEW, checking, layout>;

// Sparse 2D
template <typename T, Checking checking = _default_checking>
using SparseMatrix = GenericTensor<T, Dimension::MATRIX, Type::SPARSE, Ownership::CONTAINER, checking, Layout::SPARSE>;

template <typename T, Checking checking = _default_checking>
using SparseMatrixView = GenericTensor<T, Dimension::MATRIX, Type::SPARSE, Ownership::VIEW, checking, Layout::SPARSE>;

template <typename T, Checking checking = _default_checking>
using ConstSparseMatrixView =
    GenericTensor<T, Dimension::MATRIX, Type::SPARSE, Ownership::CONST_VIEW, checking, Layout::SPARSE>;

} // namespace utl::mvl

// ========= header guard end =========

#endif
#endif






// -----------------------------
// --------- utl::stre ---------
// -----------------------------
#if !defined(UTL_PICK_MODULES) || defined(UTLMODULE_STRE)
#ifndef UTLHEADERGUARD_STRE
#define UTLHEADERGUARD_STRE

// ======== header guard start ========

#include <iomanip>
#include <sstream>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>

namespace utl::stre {

// _________ DEVELOPER DOCS _________

// String utility extensions, mainly a template ::to_str() method which works with all STL containers,
// including maps, sets and tuples with any level of mutual nesting. Also includes some
// expansions of <type_traits> header that allow categorizing types at compile-time.
//
// # ::is_printable<Type> #
// Integral constant, returns in "::value" whether Type can be printed through std::cout.
// Criteria: Existance of operator 'ANY_TYPE operator<<(std::ostream&, Type&)'
//
// # ::is_iterable_through<Type> #
// Integral constant, returns in "::value" whether Type can be iterated through.
// Criteria: Existance of .begin() and .end() with applicable operator()++
//
// # ::is_const_iterable_through<Type> #
// Integral constant, returns in "::value" whether Type can be const-iterated through.
// Criteria: Existance of .cbegin() and .cend() with applicable operator()++
//
// # ::is_tuple_like<Type> #
// Integral constant, returns in "::value" whether Type has a tuple-like structure.
// Tuple-like structure include std::tuple, std::pair, std::array, std::ranges::subrange (since C++20)
// Criteria: Existance of applicable std::get<0>() and std::tuple_size()
//
// # ::is_string<Type> #
// Integral constant, returns in "::value" whether Type is a char string.
// Criteria: Type can be decayed to std::string or a char* pointer
//
// # ::is_to_str_convertible<Type> #
// Integral constant, returns in "::value" whether Type can be converted to string through ::to_str().
// Criteria: Existance of a valid utl::stre::to_str() overload
//
// # ::to_str() #
// Converts any standard container or a custom container with necessary member functions to std::string.
// Works with tuples and tuple-like classes.
// Works with nested containers/tuples through recursive template instantiation, which
// resolves as long as types at the end of recursion have a valid operator<<() for ostreams.
//
// # ::InlineSStream #
// Inline 'std::stringstream' construction with implicit conversion to 'std::string'.
// Rather unperformant, but convenient for using stream formating during string construction.
// Example: std::string str = (stre::InlineStream() << "Value " << 3.14 << " is smaller than " << 6.28);
//
// # ::repeat_symbol(), ::repeat_string() #
// Repeats character/string a given number of times.
//
// # ::pad_with_zeroes() #
// Pads given integer with zeroes untill a certain lenght.
// Useful when saving data in files like 'data_0001.txt', 'data_0002.txt', '...' so they get properly sorted.

// _________ IMPLEMENTATION _________

// --- type trait: is_printable ---
template <typename Type, typename = void>
struct is_printable : std::false_type {};

template <typename Type>
struct is_printable<Type, std::void_t<decltype(std::declval<std::ostream&>() << std::declval<Type>())>>
    : std::true_type {};

// --- type trait: is_iterable_through ---
template <typename Type, typename = void, typename = void>
struct is_iterable_through : std::false_type {};

template <typename Type>
struct is_iterable_through<Type, std::void_t<decltype(std::declval<Type>().begin().operator++())>,
                           std::void_t<decltype(std::declval<Type>().end().operator++())>> : std::true_type {};

// --- type trait: is_const_iterable_through ---
template <typename Type, typename = void, typename = void>
struct is_const_iterable_through : std::false_type {};

template <typename Type>
struct is_const_iterable_through<Type, std::void_t<decltype(std::declval<Type>().cbegin().operator++())>,
                                 std::void_t<decltype(std::declval<Type>().cend().operator++())>> : std::true_type {};

// --- type trait: is_tuple_like ---
template <typename Type, typename = void, typename = void>
struct is_tuple_like : std::false_type {};

template <typename Type>
struct is_tuple_like<Type, std::void_t<decltype(std::get<0>(std::declval<Type>()))>,
                     std::void_t<decltype(std::tuple_size<Type>::value)>> : std::true_type {};

// --- type trait: is_string ---
template <typename Type>
struct is_string : std::is_convertible<Type, std::string_view> {};

// --- to_str() ---

// - delimers -
constexpr auto _CONTAINER_DELIMER_L = "[ ";
constexpr auto _CONTAINER_DELIMER_M = ", ";
constexpr auto _CONTAINER_DELIMER_R = " ]";
constexpr auto _TUPLE_DELIMER_L     = "< ";
constexpr auto _TUPLE_DELIMER_M     = ", ";
constexpr auto _TUPLE_DELIMER_R     = " >";

// - predeclarations (is_to_str_convertible) -
template <typename Type, typename = void>
struct is_to_str_convertible : std::false_type {};
// false_type half should be declared before to_str() to resolve circular dependency

// - predeclarations (to_str(tuple)) -
template <template <typename... Params> class TupleLikeType, typename... Args,
          std::enable_if_t<!utl::stre::is_printable<TupleLikeType<Args...>>::value, bool>              = true,
          std::enable_if_t<!utl::stre::is_const_iterable_through<TupleLikeType<Args...>>::value, bool> = true,
          std::enable_if_t<utl::stre::is_tuple_like<TupleLikeType<Args...>>::value, bool>              = true>
std::string to_str(const TupleLikeType<Args...>& tuple);
// predeclare to resolve circular dependency between to_str(container) and to_str(tuple)

// - to_str(printable) -
template <typename Type, std::enable_if_t<utl::stre::is_printable<Type>::value, bool> = true>
std::string to_str(const Type& value) {
    std::stringstream ss;

    ss << value;

    return ss.str();
}

// - to_str(container) -
template <typename ContainerType, std::enable_if_t<!utl::stre::is_printable<ContainerType>::value, bool> = true,
          std::enable_if_t<utl::stre::is_const_iterable_through<ContainerType>::value, bool>                  = true,
          std::enable_if_t<utl::stre::is_to_str_convertible<typename ContainerType::value_type>::value, bool> = true>
std::string to_str(const ContainerType& container) {

    std::stringstream ss;

    // Special case for empty containers
    if (container.cbegin() == container.cend()) {
        ss << _CONTAINER_DELIMER_L << _CONTAINER_DELIMER_R;
        return ss.str();
    }

    // Iterate throught the container 'looking forward' by one step
    // so we can know not to place the last delimer. Using -- or std::prev()
    // is not and option since we only require iterators to be forward-iterable
    ss << _CONTAINER_DELIMER_L;

    auto it_next = container.cbegin(), it = it_next++;
    for (; it_next != container.cend(); ++it_next, ++it) ss << utl::stre::to_str(*it) << _CONTAINER_DELIMER_M;

    ss << utl::stre::to_str(*it) << _CONTAINER_DELIMER_R;

    return ss.str();
}

// - to_str(tuple) helpers -
template <typename TupleElemType>
std::string _deduce_and_perform_string_conversion(const TupleElemType& elem) {
    std::stringstream temp_ss;

    if constexpr (utl::stre::is_to_str_convertible<TupleElemType>::value) temp_ss << utl::stre::to_str(elem);
    else temp_ss << elem;

    return temp_ss.str();
}

template <typename TupleLikeType, std::size_t... Is>
void _print_tuple_fold(std::stringstream& ss, const TupleLikeType& tuple, std::index_sequence<Is...>) {
    ((ss << (Is == 0 ? "" : _TUPLE_DELIMER_M) << _deduce_and_perform_string_conversion(std::get<Is>(tuple))), ...);
} // prints tuple to stream

// - to_str(tuple) -
template <template <typename... Params> class TupleLikeType, typename... Args,
          std::enable_if_t<!utl::stre::is_printable<TupleLikeType<Args...>>::value, bool>,
          std::enable_if_t<!utl::stre::is_const_iterable_through<TupleLikeType<Args...>>::value, bool>,
          std::enable_if_t<utl::stre::is_tuple_like<TupleLikeType<Args...>>::value, bool>>
std::string to_str(const TupleLikeType<Args...>& tuple) {

    std::stringstream ss;

    // Print tuple using C++17 variadic folding with index sequence
    ss << _TUPLE_DELIMER_L;
    _print_tuple_fold(ss, tuple, std::index_sequence_for<Args...>{});
    ss << _TUPLE_DELIMER_R;

    return ss.str();
}

// - is_to_str_convertible -
template <typename Type>
struct is_to_str_convertible<Type, std::void_t<decltype(utl::stre::to_str(std::declval<Type>()))>> : std::true_type {};


// --- Inline stringstream ---
class InlineStream {
public:
    template <typename Type>
    InlineStream& operator<<(const Type& arg) {
        this->ss << arg;
        return *this;
    }

    inline operator std::string() const { return this->ss.str(); }

private:
    std::stringstream ss;
};

// Repeat symbol/string
inline std::string repeat_symbol(char symbol, size_t repeats) { return std::string(repeats, symbol); }
inline std::string repeat_string(std::string_view str, size_t repeats) {
    std::string res;
    res.reserve(str.size() * repeats);
    while (repeats--) res += str;
    return res;
}

// - Small utils -
template <typename IntegerType, std::enable_if_t<std::is_integral<IntegerType>::value, bool> = true>
std::string pad_with_zeroes(IntegerType number, std::streamsize total_size = 10) {
    std::stringstream ss;
    ss << std::setfill('0') << std::setw(total_size) << number;
    return ss.str();
}

} // namespace utl::stre

// ========= header guard end =========

#endif
#endif






// ------------------------------
// --------- utl::table ---------
// ------------------------------
#if !defined(UTL_PICK_MODULES) || defined(UTLMODULE_TABLE)
#ifndef UTLHEADERGUARD_TABLE
#define UTLHEADERGUARD_TABLE

// ======== header guard start ========

#include <initializer_list>
#include <iomanip>
#include <ios>
#include <iostream>
#include <ostream>
#include <string>
#include <vector>

namespace utl::table {

// _________ DEVELOPER DOCS _________

// Functions used to display data in a tabular fashion.
//
// # ::create() #
// Sets up table with given number of columns and their widths.
//
// # ::set_formats() #
// (optional) Sets up column formats for better display
//
// # ::set_ostream() #
// (optional) Select 'std::ostream' to which all output gets forwarded. By default 'std::cout' is selected.
//
// # ::NONE, ::FIXED(), ::DEFAULT(), ::SCIENTIFIC(), ::BOOL() #
// Format flags with following effects:
// > NONE          - Use default C++ formats
// > FIXED(N)      - Display floats in fixed form with N decimals, no argument assumes N = 3
// > DEFAULT(N)    - Display floats in default form with N decimals, no argument assumes N = 6
// > SCIENTIFIC(N) - Display floats in scientific form with N decimals, no argument assumes N = 3
// > BOOL          - Display booleans as text
//
// # ::cell() #
// Draws a single table cell, if multiple arguments are passed, draws each one in a new cell.
// Accepts any type with a defined "<<" ostream operator.
//
// # ::hline() #
// Draws horizontal line in a table. Similar to LaTeX '\hline'.

// _________ IMPLEMENTATION _________

// --- Types ---
using uint       = std::streamsize;
using _ios_flags = std::ios_base::fmtflags;

struct ColumnFormat {
    _ios_flags flags;
    uint       precision;
};

struct _Column {
    uint         width;
    ColumnFormat col_format;
};


// --- Predefined formats ---
constexpr ColumnFormat NONE = {std::ios::showpoint, 6};

constexpr ColumnFormat FIXED(uint decimals = 3) { return {std::ios::fixed, decimals}; }
constexpr ColumnFormat DEFAULT(uint decimals = 6) { return {std::ios::showpoint, decimals}; }
constexpr ColumnFormat SCIENTIFIC(uint decimals = 3) { return {std::ios::scientific, decimals}; }

constexpr ColumnFormat BOOL = {std::ios::boolalpha, 3};


// --- Internal state ---
inline std::vector<_Column> _columns;
inline std::size_t          _current_column = 0;
inline std::ostream*        _output_stream  = &std::cout;


// --- Table setup ---
inline void create(std::initializer_list<uint>&& widths) {
    _columns.resize(widths.size());
    for (size_t i = 0; i < _columns.size(); ++i) {
        _columns[i].width      = widths.begin()[i];
        _columns[i].col_format = DEFAULT();
    }
}

inline void set_formats(std::initializer_list<ColumnFormat>&& formats) {
    for (size_t i = 0; i < _columns.size(); ++i) _columns[i].col_format = formats.begin()[i];
}

inline void set_ostream(std::ostream& new_ostream) { _output_stream = &new_ostream; }

inline void cell(){};

template <typename T, typename... Types>
void cell(T value, const Types... other_values) {
    const std::string left_cline      = (_current_column == 0) ? "|" : "";
    const std::string right_cline     = (_current_column == _columns.size() - 1) ? "|\n" : "|";
    const _ios_flags  format          = _columns[_current_column].col_format.flags;
    const uint        float_precision = _columns[_current_column].col_format.precision;

    // Save old stream state
    std::ios old_state(nullptr);
    old_state.copyfmt(*_output_stream);

    // Set table formatting
    (*_output_stream) << std::resetiosflags((*_output_stream).flags());
    (*_output_stream).flags(format);
    (*_output_stream).precision(float_precision);

    // Print
    (*_output_stream) << left_cline << std::setw(_columns[_current_column].width) << value << right_cline;

    // Return old stream state
    (*_output_stream).copyfmt(old_state);

    // Advance column counter
    _current_column = (_current_column == _columns.size() - 1) ? 0 : _current_column + 1;

    cell(other_values...);
}

inline void hline() {
    (*_output_stream) << "|";
    for (const auto& col : _columns) (*_output_stream) << std::string(static_cast<size_t>(col.width), '-') << "|";
    (*_output_stream) << "\n";
}
} // namespace utl::table

// ========= header guard end =========

#endif
#endif






// ------------------------------
// --------- utl::timer ---------
// ------------------------------
#if !defined(UTL_PICK_MODULES) || defined(UTLMODULE_TIMER)
#ifndef UTLHEADERGUARD_TIMER
#define UTLHEADERGUARD_TIMER

// ======== header guard start ========

#include <chrono>
#include <ctime>
#include <string>

namespace utl::timer {

// _________ DEVELOPER DOCS _________

// Global-state timer with built-in formatting. Functions for local date and time.
//
// # ::start() #
// Starts the timer.
// Note: If start() wasn't called system will use uninitialized value as a start.
//
// # ::elapsed_ms(), ::elapsed_sec(), ::elapsed_min(), ::elapsed_hours() #
// Elapsed time as double.
//
// # ::elapsed_string_ms(), ::elapsed_string_sec(), ::elapsed_string_min(), ::elapsed_string_hours() #
// Elapsed time as std::string.
//
// # ::elapsed_string:fullform() #
// Elapsed time in format "%H hours %M min %S sec %MS ms".
//
// # ::datetime_string() #
// Current local date and time in format "%y-%m-%d %H:%M:%S".
//
// # ::datetime_string_id() #
// Current local date and time in format "%y-%m-%d-%H-%M-%S".
// Less readable that usual format, but can be used in filenames which prohibit ":" usage.

// _________ IMPLEMENTATION _________

// --- Internals ---
using _clock     = std::chrono::steady_clock;
using _chrono_ns = std::chrono::nanoseconds;

constexpr double NS_IN_MS = 1e6;

constexpr long long MS_IN_SEC  = 1000;
constexpr long long MS_IN_MIN  = 60 * MS_IN_SEC;
constexpr long long MS_IN_HOUR = 60 * MS_IN_MIN;

inline _clock::time_point _start_timepoint;

inline double _elapsed_time_as_ms() {
    const auto elapsed = std::chrono::duration_cast<_chrono_ns>(_clock::now() - _start_timepoint).count();
    return static_cast<double>(elapsed) / NS_IN_MS;
}

inline void start() { _start_timepoint = _clock::now(); }


// --- Elapsed time ---
inline double elapsed_ms() { return _elapsed_time_as_ms(); }
inline double elapsed_sec() { return _elapsed_time_as_ms() / static_cast<double>(MS_IN_SEC); }
inline double elapsed_min() { return _elapsed_time_as_ms() / static_cast<double>(MS_IN_MIN); }
inline double elapsed_hours() { return _elapsed_time_as_ms() / static_cast<double>(MS_IN_HOUR); }


// --- Elapsed string ---
inline std::string elapsed_string_ms() { return std::to_string(elapsed_ms()) + " ms"; }
inline std::string elapsed_string_sec() { return std::to_string(elapsed_sec()) + " sec"; }
inline std::string elapsed_string_min() { return std::to_string(elapsed_min()) + " min"; }
inline std::string elapsed_string_hours() { return std::to_string(elapsed_hours()) + " hours"; }


// --- Elapsed string (full form) ---
inline std::string elapsed_string_fullform() {
    long long unaccounted_ms = static_cast<long long>(_elapsed_time_as_ms());

    long long ms    = 0;
    long long min   = 0;
    long long sec   = 0;
    long long hours = 0;

    if (unaccounted_ms > MS_IN_HOUR) {
        hours += unaccounted_ms / MS_IN_HOUR;
        unaccounted_ms -= hours * MS_IN_HOUR;
    }

    if (unaccounted_ms > MS_IN_MIN) {
        min += unaccounted_ms / MS_IN_MIN;
        unaccounted_ms -= min * MS_IN_MIN;
    }

    if (unaccounted_ms > MS_IN_SEC) {
        sec += unaccounted_ms / MS_IN_SEC;
        unaccounted_ms -= sec * MS_IN_SEC;
    }

    ms = unaccounted_ms;

    return std::to_string(hours) + " hours " + std::to_string(min) + " min " + std::to_string(sec) + " sec " +
           std::to_string(ms) + " ms ";
}


// --- Localtime string ---
// SFINAE to select localtime_s() or localtime_r()
template <typename Arg_tm, typename Arg_time_t>
auto _available_localtime_impl(Arg_tm time_moment, Arg_time_t timer)
    -> decltype(localtime_s(std::forward<Arg_tm>(time_moment), std::forward<Arg_time_t>(timer))) {
    return localtime_s(std::forward<Arg_tm>(time_moment), std::forward<Arg_time_t>(timer));
}

template <typename Arg_tm, typename Arg_time_t>
auto _available_localtime_impl(Arg_tm time_moment, Arg_time_t timer)
    -> decltype(localtime_r(std::forward<Arg_time_t>(timer), std::forward<Arg_tm>(time_moment))) {
    return localtime_r(std::forward<Arg_time_t>(timer), std::forward<Arg_tm>(time_moment));
}

inline std::string _datetime_string_with_format(const char* format) {
    std::time_t timer = std::time(nullptr);
    std::tm     time_moment{};

    // Call localtime_s() or localtime_r() depending on which one is present
    _available_localtime_impl(&time_moment, &timer);

    // // Macro version, can be used instead of SFINAE resolution
    // // Get localtime safely (if possible)
    // #if defined(__unix__)
    // localtime_r(&timer, &time_moment);
    // #elif defined(_MSC_VER)
    // localtime_s(&time_moment, &timer);
    // #else
    // // static std::mutex mtx; // mutex can be used to make thread-safe version but add another dependency
    // // std::lock_guard<std::mutex> lock(mtx);
    // time_moment = *std::localtime(&timer);
    // #endif

    // Convert time to C-string
    char mbstr[100];
    std::strftime(mbstr, sizeof(mbstr), format, &time_moment);

    return std::string(mbstr);
}

inline std::string datetime_string() { return _datetime_string_with_format("%Y-%m-%d %H:%M:%S"); }

inline std::string datetime_string_id() { return _datetime_string_with_format("%Y-%m-%d-%H-%M-%S"); }

} // namespace utl::timer

// ========= header guard end =========

#endif
#endif





// -----------------------------------
// --------- utl::voidstream ---------
// -----------------------------------
#if !defined(UTL_PICK_MODULES) || defined(UTLMODULE_VOIDSTREAM)
#ifndef UTLHEADERGUARD_VOIDSTREAM
#define UTLHEADERGUARD_VOIDSTREAM

// ======== header guard start ========

#include <ostream>
#include <streambuf>

namespace utl::voidstream {

// _________ DEVELOPER DOCS _________

// "voidstream" that functions like std::ostream with no output.
// Can be passed to interfaces that use streams to silence their output.
//
// # ::vstreambuf #
// Stream buffer that overflows with no output, usage example:
//   > std::ofstream output_stream(&vstreambuf);
//   > output_stream << VALUE; // produces nothing
//
// # ::vout #
// Output stream that produces no output, usage example:
//   > vout << VALUE; // produces nothing

// _________ IMPLEMENTATION _________

class VoidStreamBuf : public std::streambuf {
public:
    inline int overflow(int c) { return c; }
};

class VoidStream : public std::ostream {
public:
    inline VoidStream() : std::ostream(&buffer) {}

private:
    VoidStreamBuf buffer;
};

inline VoidStreamBuf vstreambuf;
inline VoidStream    vout;
} // namespace utl::voidstream

// ========= header guard end =========

#endif
#endif






// ------------------------------
// --------- UTL_DEFINE ---------
// ------------------------------
#if !defined(UTL_PICK_MODULES) || defined(UTLMACRO_DEFINE)
#ifndef UTLHEADERGUARD_DEFINE
#define UTLHEADERGUARD_DEFINE

// ======== header guard start ========

#include <cctype>
#include <sstream>
#include <string>
#include <type_traits>

// _________ DEVELOPER DOCS _________

// Macros for automatic codegen, such as codegen for type traits, loops,
// operators, enum bitflags, enum <=> string conversions and etc.
//
// # UTL_DEFINE_CURRENT_OS_STRING #
// Evaluates to current platform name detected through compiler macros. Can evaluate to:
// "Windows64", "Windows32", "Windows (CYGWIN)", "Android", "Linux", "Unix-like OS", "MacOS", ""
//
// # UTL_DEFINE_IS_DEBUG #
// Resolves to 'true' in debug mode and to 'false' in release.
// Useful for some debug-only 'contexpr if ()' expressions.
//
// # UTL_DEFINE_REPEAT(N) #
// Equivalent to 'for(int i=0; i<N; ++i)'
//
// # UTL_DEFINE_VA_ARGS_COUNT(args...) #
// Counts how many comma-separated arguments are passed. Works in a true 'preprocessor' evaluating arguments
// as arbitrary text that doesn't have to be a valid identifier. Useful for a lot of variadic macros.
// Note: Since macro evaluates arbitrary preprocessor text, language construct with additional comma should
// be surrounded by parentheses aka '(std::pair<int, int>{4, 5})', since 'std::pair<int, int>{4, 5}' gets
// interpreted as 3 separate args 'std::pair<int', ' int>{4', ' 5}' by any variadic macro.
//
// # UTL_DEFINE_ENUM_WITH_STRING_CONVERSION(enum_name, enum_values...) #
// Create enum 'enum_name' with given enum_values and methods 'to_string()', 'from_string()' inside 'enum_name'
// namespace. Note: Enums can't be declared inside functions.
//
// # UTL_DEFINE_IS_FUNCTION_PRESENT(function_name, return_type, argumet_types...) #
// Creates integral constant named 'is_function_present_{function_name}' that
// returns in "::value" whether function with given name and mask is a valid identifier.
// Usefull for detecting existance of platform-specific methods.
// Note: This is kinda pushing the limits of compiler, but possible through SFINAE, which is the only
// mechanism that can perform that kind check. This is the reason implementing an inline FUNC_EXISTS(...)
// kind of macro is impossible and integral constant has to be created.
//
// # UTL_DEFINE_EXIT(message, exit_code) #
// Performs 'std::exit()' with some additional decorators. Displays formatted exit message, exit code,
// file, function and line of callsite to std::cerr. If ommited replaces message & exit_code with default values.

// _________ IMPLEMENTATION _________

// --- Detect OS ---
#if defined(_WIN64) // _WIN64 implies _WIN32 so it should be first
#define UTL_DEFINE_CURRENT_OS_STRING "Windows64"
#elif defined(_WIN32)
#define UTL_DEFINE_CURRENT_OS_STRING "Windows32"
#elif defined(__CYGWIN__) && !defined(_WIN32)
#define UTL_DEFINE_CURRENT_OS_STRING "Windows (CYGWIN)" // Cygwin POSIX under Microsoft Window
#elif defined(__ANDROID__)                              // __ANDROID__ implies __linux__ so it should be first
#define UTL_DEFINE_CURRENT_OS_STRING "Android"
#elif defined(__linux__)
#define UTL_DEFINE_CURRENT_OS_STRING "Linux"
#elif defined(unix) || defined(__unix__) || defined(__unix)
#define UTL_DEFINE_CURRENT_OS_STRING "Unix-like OS"
#elif defined(__APPLE__) && defined(__MACH__)
#define UTL_DEFINE_CURRENT_OS_STRING "MacOS" // Apple OSX and iOS (Darwin)
#else
#define UTL_DEFINE_CURRENT_OS_STRING ""
#endif

// --- Debug as bool constant ---
#ifdef _DEBUG
#define UTL_DEFINE_IS_DEBUG true
#else
#define UTL_DEFINE_IS_DEBUG false
#endif

// --- Repeat loop ---
#define UTL_DEFINE_REPEAT(repeats_)                                                                                    \
    for (std::remove_const<decltype(repeats_)>::type count_ = 0; count_ < repeats_; ++count_)

// --- Size of __VA_ARGS__ in variadic macros ---
#define _utl_expand_va_args(x_) x_ // a fix for MSVC bug not expanding __VA_ARGS__ properly

#define _utl_va_args_count_impl(x01_, x02_, x03_, x04_, x05_, x06_, x07_, x08_, x09_, x10_, x11_, x12_, x13_, x14_,    \
                                x15_, x16_, x17_, x18_, x19_, x20_, x21_, x22_, x23_, x24_, x25_, x26_, x27_, x28_,    \
                                x29_, x30_, x31_, x32_, x33_, x34_, x35_, x36_, x37_, x38_, x39_, x40_, x41_, x42_,    \
                                x43_, x44_, x45_, x46_, x47_, x48_, x49_, N_, ...)                                     \
    N_

#define UTL_DEFINE_VA_ARGS_COUNT(...)                                                                                  \
    _utl_expand_va_args(_utl_va_args_count_impl(__VA_ARGS__, 49, 48, 47, 46, 45, 44, 43, 42, 41, 40, 39, 38, 37, 36,   \
                                                35, 34, 33, 32, 31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19,    \
                                                18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0))

// --- Enum with string conversion ---
inline std::string _utl_trim_enum_string(const std::string& str) {
    std::string::const_iterator left_it = str.begin();
    while (left_it != str.end() && std::isspace(*left_it)) ++left_it;

    std::string::const_reverse_iterator right_it = str.rbegin();
    while (right_it.base() != left_it && std::isspace(*right_it)) ++right_it;

    return std::string(left_it, right_it.base()); // return string with whitespaces trimmed at both sides
}

inline void _utl_split_enum_args(const char* va_args, std::string* strings, int count) {
    std::stringstream ss(va_args);
    std::string       buffer;

    for (int i = 0; ss.good() && (i < count); ++i) {
        std::getline(ss, buffer, ',');
        strings[i] = _utl_trim_enum_string(buffer);
    }
};

#define UTL_DEFINE_ENUM_WITH_STRING_CONVERSION(enum_name_, ...)                                                        \
    namespace enum_name_ {                                                                                             \
    enum enum_name_ { __VA_ARGS__, _count };                                                                           \
                                                                                                                       \
    inline std::string _strings[_count];                                                                               \
                                                                                                                       \
    inline std::string to_string(enum_name_ enum_val) {                                                                \
        if (_strings[0].empty()) { _utl_split_enum_args(#__VA_ARGS__, _strings, _count); }                             \
        return _strings[enum_val];                                                                                     \
    }                                                                                                                  \
                                                                                                                       \
    inline enum_name_ from_string(const std::string& enum_str) {                                                       \
        if (_strings[0].empty()) { _utl_split_enum_args(#__VA_ARGS__, _strings, _count); }                             \
        for (int i = 0; i < _count; ++i) {                                                                             \
            if (_strings[i] == enum_str) { return static_cast<enum_name_>(i); }                                        \
        }                                                                                                              \
        return _count;                                                                                                 \
    }                                                                                                                  \
    }
    // IMPLEMENTATION COMMENTS:
    // We declare namespace with enum inside to simulate enum-class while having '_strings' array
    // and 'to_string()', 'from_string()' methods bundled with it.
    //
    // To count number of enum elements we add fake '_count' value at the end, which ends up being enum size
    //
    // '_strings' is declared compile-time, but gets filled through lazy evaluation upon first
    // 'to_string()' or 'from_string()' call. To fill it we interpret #__VA_ARGS__ as a single string
    // with some comma-separated identifiers. Those identifiers get split by commas, trimmed from
    // whitespaces and added to '_strings'
    //
    // Upon further calls (enum -> string) conversion is done though taking '_strings[enum_val]',
    // while (string -> enum) conversion requires searching through '_strings' to find enum index


// --- type trait: is_function_present ---
#define UTL_DEFINE_IS_FUNCTION_PRESENT(function_name_, return_type_, ...)                                              \
    template <typename ReturnType, typename... ArgTypes>                                                               \
    class _utl_is_function_present_impl_##function_name_ {                                                             \
    private:                                                                                                           \
        typedef char no[sizeof(ReturnType) + 1];                                                                       \
                                                                                                                       \
        template <typename... C>                                                                                       \
        static auto test(C... arg) -> decltype(function_name_(arg...));                                                \
                                                                                                                       \
        template <typename... C>                                                                                       \
        static no& test(...);                                                                                          \
                                                                                                                       \
    public:                                                                                                            \
        enum { value = (sizeof(test<ArgTypes...>(std::declval<ArgTypes>()...)) == sizeof(ReturnType)) };               \
    };                                                                                                                 \
                                                                                                                       \
    using is_function_present_##function_name_ =                                                                       \
        _utl_is_function_present_impl_##function_name_<return_type_, __VA_ARGS__>;
// TASK:
// We need to detect at compile time if function FUNC(ARGS...) exists.
// FUNC identifier isn't guaranteed to be declared.
//
// Ideal method would look like UTL_FUNC_EXISTS(FUNC, RETURN_TYPE, ARG_TYPES...) -> true/false
// This does not seem to be possible, we have to declare integral constant insted, see explanation below.
//
// WHY IS IT SO HARD:
// (1) Can this be done through preprocessor macros?
// No, preprocessor has no way to tell whether C++ identifier is defined or not.
//
// (2) Is there a compiler-specific way to do it?
// Doesn't seem to be the case.
//
// (3) Why not use some sort of template with FUNC as a parameter?
// Essentially we have to evaluate undeclared identifier, while compiler exits with error upon
// encountering anything undeclared. The only way to detect whether undeclared identifier exists
// or not seems to be through SFINAE.
//
// IMPLEMENTATION COMMENTS:
// We declate integral constant class with 2 functions 'test()', first one takes priority during overload
// resolution and compiles if FUNC(ARGS...) is defined, otherwise it's {Substitution Failure} which is
// {Is Not An Error} and second function compiles.
//
// To resolve which overload of 'test()' was selected we check the sizeof() return type, 2nd overload
// has a return type 'char[sizeof(ReturnType) + 1]' so it's always different from 1st overload.
// Resolution result (true/false) gets stored to '::value'.
//
// Note that we can't pass 'ReturnType' and 'ArgTypes' directly through '__VA_ARGS__' because
// to call function 'test(ARGS...)' in general case we have to 'std::declval<>()' all 'ARGS...'.
// To do so we can use variadic template syntax and then just forward '__VA_ARGS__' to the template
// through 'using is_function_present = is_function_present_impl<ReturnType, __VA_ARGS__>'.
//
// ALTERNATIVES: Perhaps some sort of tricky inline SFINAE can be done through C++14 generic lambdas.
//
// Note: Some versions of 'clangd' give a 'bugprone-sizeof-expression' warning for sizeof(*A),
// this is a false alarm.

// --- Exit with decorators ---
inline void _utl_exit_with_message(std::string_view file, int line, std::string_view func,
                                   std::string_view message = "<NO MESSAGE>", int code = 1) {
    constexpr int  HLINE_WIDTH  = 50;
    constexpr char HLINE_SYMBOL = '-';

    const std::string_view filename = file.substr(file.find_last_of("/\\") + 1);

    // Create some more space
    std::cerr << std::endl;

    // Horizontal line
    std::fill_n(std::ostreambuf_iterator<char>(std::cerr), HLINE_WIDTH, HLINE_SYMBOL);
    std::cerr << std::endl;

    // Text
    std::cerr << "Exit triggered on [" << filename << ":" << line << ", " << func << "()] with:" << std::endl
              << "Message => " << message << std::endl
              << "Code    => " << code << std::endl;

    // Horizontal line
    std::fill_n(std::ostreambuf_iterator<char>(std::cerr), HLINE_WIDTH, HLINE_SYMBOL);
    std::cerr << std::endl;

    std::exit(code);
}

#define UTL_DEFINE_EXIT(...) _utl_exit_with_message(__FILE__, __LINE__, __func__, __VA_ARGS__);

// TODO: UTL_DEFINE_BITFLAG_OPERATORS
// TODO: UTL_DEFINE_BITFLAG_COMPATIBILITY_CHECK

// ========= header guard end =========

#endif
#endif






// ---------------------------
// --------- UTL_LOG ---------
// ---------------------------
#if !defined(UTL_PICK_MODULES) || defined(UTLMACRO_LOG)
#ifndef UTLHEADERGUARD_LOG
#define UTLHEADERGUARD_LOG

// ======== header guard start ========

#include <iostream>
#include <ostream>
#include <string_view>

// _________ DEVELOPER DOCS _________

// # UTL_LOG_SET_OUTPUT() #
// Select ostream used by LOG macros.
//
// # UTL_LOG(), UTL_LOG_DEBUG() #
// Print message to selected ostream prefixed with [<filename>:<line> (<function>)].
// Accepts multiple args (with defined operator <<) that get concatenated into a single message.
// DEBUG version compiles to nothing in release.

// _________ IMPLEMENTATION _________

inline std::ostream* _utl_log_ostream = &std::cout;

#define UTL_LOG_SET_OUTPUT(new_stream_) _utl_log_ostream = &new_stream_;

template <typename... Args>
inline void _utl_log_print(std::string_view file, int line, std::string_view func, const Args&... args) {
    const std::string_view filename = file.substr(file.find_last_of("/\\") + 1);

    ///(*_utl_log_ostream) << "\033[31;1m"; // Supported by Linux and Windows10+, but prints to files, figure out a fix
    (*_utl_log_ostream) << "[" << filename << ":" << line << ", " << func << "()]";
    ///(*_utl_log_ostream) << "\033[0m";

    (*_utl_log_ostream) << " ";
    ((*_utl_log_ostream) << ... << args);
    (*_utl_log_ostream) << '\n';
}

#define UTL_LOG(...) _utl_log_print(__FILE__, __LINE__, __func__, __VA_ARGS__)

#ifdef _DEBUG
#define UTL_LOG_DEBUG(...) _utl_log_print(__FILE__, __LINE__, __func__, __VA_ARGS__)
#else
#define UTL_LOG_DEBUG(...)
#endif

// ========= header guard end =========

#endif
#endif






// --------------------------------
// --------- UTL_PROFILER ---------
// --------------------------------
#if !defined(UTL_PICK_MODULES) || defined(UTLMACRO_PROFILER)
#ifndef UTLHEADERGUARD_PROFILER
#define UTLHEADERGUARD_PROFILER

// ======== header guard start ========

#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <ios>
#include <iostream>
#include <map>
#include <ostream>
#include <sstream>
#include <string>
#include <string_view>

// _________ DEVELOPER DOCS _________

// Macros for quick code profiling.
//
// # UTL_PROFILER_REROUTE_TO_FILE(filepath) #
// Reroutes profiler summary from std::cout to a file with a given name.
//
// # UTL_PROFILER, UTL_PROFILER_LABELED() #
// Profiles the following expression/scope. If profiled scope was entered at any point of the program,
// upon exiting 'main()' the table with profiling results will be printed. Profiling results include:
// - Total program runtime
// - Total runtime of each profiled scope
// - % of total runtime taken by each profiled scope
// - Profiler labels (if using 'UTL_PROFILE_LABELED()', otherwise the label is set to "<NONE>")
// - Profiler locations: file, function, line

// _________ IMPLEMENTATION _________

inline std::string _utl_profiler_reroute_to_filepath = "";

#define UTL_PROFILER_REROUTE_TO_FILE(filepath_) _utl_profiler_reroute_to_filepath = filepath_;

using _utl_profiler_clock         = std::chrono::steady_clock;
using _utl_profiler_time_duration = _utl_profiler_clock::duration;
using _utl_profiler_time_point    = _utl_profiler_clock::time_point;

inline const _utl_profiler_time_point _utl_profiler_program_init_time_point = _utl_profiler_clock::now();
// automatically gets program launch time so we can compute total runtime later

inline std::string _utl_profiler_format_call_site(std::string_view file, int line, std::string_view func) {
    const std::string_view filename = file.substr(file.find_last_of("/\\") + 1);

    std::stringstream ss;
    ss << filename << ":" << line << ", " << func << "()";
    return ss.str();
}

struct _utl_profiler_record {
    std::string                 label;
    _utl_profiler_time_duration duration;
};

inline void _utl_profiler_atexit(); // predeclare, it needs '_utl_profiler' but used by '_utl_profiler'

class _utl_profiler {
private:
    std::string              call_site;
    std::string              label;
    _utl_profiler_time_point construction_time_point;

public:
    inline static std::map<std::string, _utl_profiler_record>
        records; // std::map because we do WANT sorting by call-site name

    operator bool() const { return true; } // needed so we can use 'if (auto x = _utl_profiler())' construct

public:
    _utl_profiler(std::string_view file, int line, std::string_view func, std::string_view label) {
        this->call_site               = _utl_profiler_format_call_site(file, line, func);
        this->label                   = label;
        this->construction_time_point = _utl_profiler_clock::now();

        // If profiler ever gets called => registed results output at std::exit()
        static bool first_call = true;
        if (first_call) {
            std::atexit(_utl_profiler_atexit);
            first_call = false;
        }
    }

    ~_utl_profiler() {
        const auto it                = _utl_profiler::records.find(this->call_site);
        const auto profiled_duration = _utl_profiler_clock::now() - this->construction_time_point;

        // Record with the same callsite exists => accumulate duration
        if (it != _utl_profiler::records.end()) {
            (*it).second.duration += profiled_duration;
        }
        // Otherwise => add new record with duration
        else {
            _utl_profiler::records.insert({
                std::string(this->call_site), _utl_profiler_record{this->label, profiled_duration}
            });
        }
    }
};

inline void _utl_profiler_atexit() {
    namespace chr = std::chrono;

    const auto total_runtime = _utl_profiler_clock::now() - _utl_profiler_program_init_time_point;
    // const auto total_runtime = std::chrono::duration_cast<std::chrono::milliseconds>();

    // Choose whether to print or reroute output to file
    std::ostream* ostr = &std::cout;
    std::ofstream output_file;

    if (!_utl_profiler_reroute_to_filepath.empty()) {
        output_file.open(_utl_profiler_reroute_to_filepath);
        ostr = &output_file;
    }

    // Convenience functions
    const auto duration_to_sec = [](_utl_profiler_time_duration duration) -> double {
        return chr::duration_cast<chr::nanoseconds>(duration).count() / 1e9;
    };

    const auto duration_percentage = [&](double duration_sec) -> double {
        const double total_runtime_sec = duration_to_sec(total_runtime);
        return (duration_sec / total_runtime_sec) * 100.;
    };

    const auto float_printed_size = [](double value, std::streamsize precision, decltype(std::fixed) format,
                                       std::string_view postfix) -> std::streamsize {
        std::stringstream ss;
        ss << std::setprecision(precision) << format << value << postfix;
        return ss.str().size(); // can be done faster but we don't really care here
    };

    const auto repeat_hline_symbol = [](std::streamsize repeats) -> std::string {
        return std::string(static_cast<size_t>(repeats), '-');
    };

    constexpr std::streamsize  duration_precision = 2;
    constexpr auto             duration_format    = std::fixed;
    constexpr std::string_view duration_postfix   = " s";

    constexpr std::streamsize  percentage_precision = 1;
    constexpr auto             percentage_format    = std::fixed;
    constexpr std::string_view percentage_postfix   = "%";

    // Collect max length of each column (for proper formatting)
    constexpr std::string_view column_name_call_site  = "Call Site";
    constexpr std::string_view column_name_label      = "Label";
    constexpr std::string_view column_name_duration   = "Time";
    constexpr std::string_view column_name_percentage = "Time %";

    std::streamsize max_length_call_site  = column_name_call_site.size();
    std::streamsize max_length_label      = column_name_label.size();
    std::streamsize max_length_duration   = column_name_duration.size();
    std::streamsize max_length_percentage = column_name_percentage.size();

    for (const auto& record : _utl_profiler::records) {
        const auto&  call_site    = record.first;
        const auto&  label        = record.second.label;
        const double duration_sec = duration_to_sec(record.second.duration);

        // 'Call Site' column
        const std::streamsize length_call_site = call_site.size();
        if (max_length_call_site < length_call_site) max_length_call_site = length_call_site;

        // 'Label' column
        const std::streamsize length_label = label.size();
        if (max_length_label < length_label) max_length_label = length_label;

        // 'Time' column
        const std::streamsize length_duration =
            float_printed_size(duration_sec, duration_precision, duration_format, duration_postfix);
        if (max_length_duration < length_duration) max_length_duration = length_duration;

        // 'Time %' column
        const auto            percentage = duration_percentage(duration_sec);
        const std::streamsize length_percentage =
            float_printed_size(percentage, percentage_precision, percentage_format, percentage_postfix);
        if (max_length_percentage < length_percentage) max_length_percentage = length_percentage;
    }

    // Print formatted profiler header
    constexpr std::string_view HEADER_TEXT = " UTL PROFILING RESULTS ";

    const std::streamsize total_table_length = sizeof("| ") - 1 + max_length_call_site + sizeof(" | ") - 1 +
                                               max_length_label + sizeof(" | ") - 1 + max_length_duration +
                                               sizeof(" | ") - 1 + max_length_percentage + sizeof(" |") -
                                               1; // -1 because sizeof(char[]) accounts for invisible '\0' at the end

    const std::streamsize header_text_length = HEADER_TEXT.size();
    const std::streamsize header_left_pad    = (total_table_length - header_text_length) / 2;
    const std::streamsize header_right_pad   = total_table_length - header_text_length - header_left_pad;

    *ostr << "\n"
          << repeat_hline_symbol(header_left_pad + 1) << HEADER_TEXT << repeat_hline_symbol(header_right_pad + 1)
          << '\n'
          // + 1 makes header hline extend 1 character past the table on both sides
          << "\n"
          << " Total runtime -> " << std::setprecision(duration_precision) << duration_format
          << duration_to_sec(total_runtime) << " sec\n"
          << "\n";

    // Print formatted table header
    *ostr << " | " << std::setw(max_length_call_site) << column_name_call_site << " | " << std::setw(max_length_label)
          << column_name_label << " | " << std::setw(max_length_duration) << column_name_duration << " | "
          << std::setw(max_length_percentage) << column_name_percentage << " |\n";

    *ostr << " |"
          << repeat_hline_symbol(max_length_call_site + 2) // add 2 to account for delimers not having spaces in hline
          << "|" << repeat_hline_symbol(max_length_label + 2) << "|" << repeat_hline_symbol(max_length_duration + 2)
          << "|" << repeat_hline_symbol(max_length_percentage + 2) << "|\n";

    *ostr << std::setfill(' '); // reset the fill so we don't mess with table contents


    // Print formatted table contents
    for (const auto& record : _utl_profiler::records) {
        const auto&  call_site    = record.first;
        const auto&  label        = record.second.label;
        const double duration_sec = duration_to_sec(record.second.duration);
        const double percentage   = duration_percentage(duration_sec);

        // Joint floats with their postfixes into a single string so they are properly handled by std::setw()
        // (which only affects the first value leading to a table misaligned by postfix size)
        std::stringstream ss_duration;
        ss_duration << std::setprecision(duration_precision) << duration_format << duration_sec << duration_postfix;

        std::stringstream ss_percentage;
        ss_percentage << std::setprecision(percentage_precision) << percentage_format << percentage
                      << percentage_postfix;

        *ostr << " | " << std::setw(max_length_call_site) << call_site << " | " << std::setw(max_length_label) << label
              << " | " << std::setw(max_length_duration) << ss_duration.str() << " | "
              << std::setw(max_length_percentage) << ss_percentage.str() << " |\n";
    }
}

#define _utl_profiler_concat_tokens(a, b) a##b
#define _utl_profiler_concat_tokens_wrapper(a, b) _utl_profiler_concat_tokens(a, b)
#define _utl_profiler_add_line_number_to_variable_name(varname_) _utl_profiler_concat_tokens_wrapper(varname_, __LINE__)
// This macro creates token 'varname_##__LINE__' from 'varname_'.
//
// The reason we can't just write it as is, is that function-macros only expands their macro-arguments
// if neither the stringizing operator # nor the token-pasting operator ## are applied to the arguments
// inside the macro body.
//
// Which means in a simple 'varname_##__LINE__' macro '__LINE__' doesn't expand to it's value.
//
// We can get around this fact by introducing indirection,
// '__LINE__' gets expanded in '_utl_profiler_concat_tokens_wrapper()'
// and then tokenized and concatenated in '_utl_profiler_concat_tokens()'

#define UTL_PROFILER_LABELED(label_)                                                                                   \
    if (auto _utl_profiler_add_line_number_to_variable_name(profiler_) =                                               \
            _utl_profiler(__FILE__, __LINE__, __func__, label_))
// We add line number to profiler variable name to prevent nested profiler scopes from shadowing
// each other's 'profiler_' variable. While such shadowing has no effect on an actual behavior,
// it does cause a warning from most compilers.
//
// Such solution isn't perfect (2 scopes can be nested on the same line in deffirent files), but it
// is good enough and there is no way to get a better 'unique varname' without resorting to
// specific compiler extensions.

#define UTL_PROFILER UTL_PROFILER_LABELED("<NONE>")

// ========= header guard end =========

#endif
#endif


// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
// $$$$$$$$$ Module Format Blueprint $$$$$$$$$
// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

// ----------------------------
// --------- utl::XXX ---------
// ----------------------------
#if !defined(UTL_PICK_MODULES) || defined(UTLMODULE_XXX)
#ifndef UTLHEADERGUARD_XXX
#define UTLHEADERGUARD_XXX

// ======== header guard start ========

// NOTE: INCLUDES

namespace utl::xxx {

// _________ DEVELOPER DOCS _________

// NOTE: DOCS

// _________ IMPLEMENTATION _________

// NOTE: IMPL

} // namespace utl::xxx

// ========= header guard end =========

#endif
#endif

// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
// $$$$$$$$$ Macro-Module Format Blueprint $$$$$$$$$
// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

// ---------------------------
// --------- UTL_XXX ---------
// ---------------------------
#if !defined(UTL_PICK_MODULES) || defined(UTLMACRO_XXX)
#ifndef UTLHEADERGUARD_XXX
#define UTLHEADERGUARD_XXX

// ======== header guard start ========

// NOTE: INCLUDES

// _________ DEVELOPER DOCS _________

// NOTE: DOCS

// _________ IMPLEMENTATION _________

// NOTE: IMPL

// ========= header guard end =========

#endif
#endif