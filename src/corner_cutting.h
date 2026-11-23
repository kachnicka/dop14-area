#pragma once

// Corner cutting, as the listing in Kacerik and Bittner 2024, Appendix A,
// with its edge-overlap constant halved (Eq. 3 and 4). T rounds as written
// (the build uses -ffp-contract=off). The formula is exact on tight slabs: its
// error is the fit's rounding, amplified, plus its own.

#include "dop14.h"

#include <algorithm>
#include <array>
#include <numbers>

// corner-cut depths and edge-overlap terms, for the tests
template<typename T>
struct AreaTerms {
    std::array<T, 8>  d;
    std::array<T, 12> s;
};

namespace corner_cutting {

// (3 - sqrt(3)) / 2 per corner, (2 - sqrt(3)) / 4 per edge overlap (the
// listing has twice that), each rounded once to T
template<typename T>
inline constexpr T cCut { static_cast<T>((3.L - std::numbers::sqrt3_v<fext>) / 2.L) };
template<typename T>
inline constexpr T cOverlap { static_cast<T>((2.L - std::numbers::sqrt3_v<fext>) / 4.L) };
static_assert(cCut<f32> == .6339745962155614f && cCut<f64> == .6339745962155614);
static_assert(cOverlap<f32> == .06698729810778067f && cOverlap<f64> == .06698729810778067);

}

template<typename T>
[[nodiscard]] T CornerCut(Dop14<T> const& dop, AreaTerms<T>* terms = nullptr)
{
    Vec3<T> const diag { .x = (dop[1] - dop[0]), .y = (dop[3] - dop[2]), .z = (dop[5] - dop[4]) };
    T const       aabbSA = T(2) * (diag.x * diag.y + diag.y * diag.z + diag.z * diag.x);

    T d[8];
    d[0] = dop[6] - DotN3(Vec3<T> { .x = dop[0], .y = dop[2], .z = dop[4] });
    d[1] = DotN3(Vec3<T> { .x = dop[1], .y = dop[3], .z = dop[5] }) - dop[7];
    d[2] = dop[8] - DotN4(Vec3<T> { .x = dop[0], .y = dop[2], .z = dop[5] });
    d[3] = DotN4(Vec3<T> { .x = dop[1], .y = dop[3], .z = dop[4] }) - dop[9];
    d[4] = dop[10] - DotN5(Vec3<T> { .x = dop[0], .y = dop[3], .z = dop[4] });
    d[5] = DotN5(Vec3<T> { .x = dop[1], .y = dop[2], .z = dop[5] }) - dop[11];
    d[6] = dop[12] - DotN6(Vec3<T> { .x = dop[0], .y = dop[3], .z = dop[5] });
    d[7] = DotN6(Vec3<T> { .x = dop[1], .y = dop[2], .z = dop[4] }) - dop[13];

    T tetrahedraCorners = d[0] * d[0] + d[1] * d[1] + d[2] * d[2] + d[3] * d[3] + d[4] * d[4] + d[5] * d[5] + d[6] * d[6] + d[7] * d[7];
    tetrahedraCorners *= corner_cutting::cCut<T>;

    T s[12];
    s[0]  = std::max(T(0), d[0] + d[7] - diag.x);
    s[1]  = std::max(T(0), d[1] + d[6] - diag.x);
    s[2]  = std::max(T(0), d[2] + d[5] - diag.x);
    s[3]  = std::max(T(0), d[3] + d[4] - diag.x);
    s[4]  = std::max(T(0), d[0] + d[4] - diag.y);
    s[5]  = std::max(T(0), d[1] + d[5] - diag.y);
    s[6]  = std::max(T(0), d[2] + d[6] - diag.y);
    s[7]  = std::max(T(0), d[3] + d[7] - diag.y);
    s[8]  = std::max(T(0), d[0] + d[2] - diag.z);
    s[9]  = std::max(T(0), d[1] + d[3] - diag.z);
    s[10] = std::max(T(0), d[4] + d[6] - diag.z);
    s[11] = std::max(T(0), d[5] + d[7] - diag.z);

    T tetrahedraCorrection = s[0] * s[0] + s[1] * s[1] + s[2] * s[2] + s[3] * s[3]
        + s[4] * s[4] + s[5] * s[5] + s[6] * s[6] + s[7] * s[7]
        + s[8] * s[8] + s[9] * s[9] + s[10] * s[10] + s[11] * s[11];
    tetrahedraCorrection *= corner_cutting::cOverlap<T>;

    if (terms) {
        std::copy_n(d, 8, terms->d.begin());
        std::copy_n(s, 12, terms->s.begin());
    }
    return aabbSA - tetrahedraCorners + tetrahedraCorrection;
}
