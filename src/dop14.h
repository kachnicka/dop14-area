#pragma once

// The 14-DOP and its fit, in f32 or f64:
// diagonal sums left to right, slabs by std::min and std::max.

#include "types.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <vector>

template<typename T>
struct Vec3 {
    T x;
    T y;
    T z;
};

inline constexpr u32 cSlabs { 7 };

// [min, max] interleaved over x, y, z, x+y+z, x+y-z, x-y+z, x-y-z (raw sums)
template<typename T>
using Dop14 = std::array<T, 2 * cSlabs>;

inline constexpr std::array<Vec3<i32>, cSlabs> cNormals { {
    { .x = 1, .y = 0, .z = 0 },
    { .x = 0, .y = 1, .z = 0 },
    { .x = 0, .y = 0, .z = 1 },
    { .x = 1, .y = 1, .z = 1 },
    { .x = 1, .y = 1, .z = -1 },
    { .x = 1, .y = -1, .z = 1 },
    { .x = 1, .y = -1, .z = -1 },
} };

template<typename T>
[[nodiscard]] constexpr T DotN0(Vec3<T> v)
{
    return v.x;
}
template<typename T>
[[nodiscard]] constexpr T DotN1(Vec3<T> v)
{
    return v.y;
}
template<typename T>
[[nodiscard]] constexpr T DotN2(Vec3<T> v)
{
    return v.z;
}
template<typename T>
[[nodiscard]] constexpr T DotN3(Vec3<T> v)
{
    return (v.x + v.y + v.z);
}
template<typename T>
[[nodiscard]] constexpr T DotN4(Vec3<T> v)
{
    return (v.x + v.y - v.z);
}
template<typename T>
[[nodiscard]] constexpr T DotN5(Vec3<T> v)
{
    return (v.x - v.y + v.z);
}
template<typename T>
[[nodiscard]] constexpr T DotN6(Vec3<T> v)
{
    return (v.x - v.y - v.z);
}

template<typename T>
[[nodiscard]] constexpr std::array<T, cSlabs> Project(Vec3<T> v)
{
    return { DotN0(v), DotN1(v), DotN2(v), DotN3(v), DotN4(v), DotN5(v), DotN6(v) };
}

template<typename T>
[[nodiscard]] constexpr Dop14<T> Init(Vec3<T> v)
{
    auto const p { Project(v) };
    Dop14<T>   dop {};
    for (u32 i { 0 }; i < cSlabs; i++) {
        dop[2 * i]     = p[i];
        dop[2 * i + 1] = p[i];
    }
    return dop;
}

template<typename T>
[[nodiscard]] constexpr Dop14<T> Init(Vec3<T> v0, Vec3<T> v1, Vec3<T> v2)
{
    auto const p0 { Project(v0) };
    auto const p1 { Project(v1) };
    auto const p2 { Project(v2) };
    Dop14<T>   dop {};
    for (u32 i { 0 }; i < cSlabs; i++) {
        dop[2 * i]     = std::min(p0[i], std::min(p1[i], p2[i]));
        dop[2 * i + 1] = std::max(p0[i], std::max(p1[i], p2[i]));
    }
    return dop;
}

template<typename T>
constexpr void Fit(Dop14<T>& dop, Vec3<T> v)
{
    auto const p { Project(v) };
    for (u32 i { 0 }; i < cSlabs; i++) {
        dop[2 * i]     = std::min(dop[2 * i], p[i]);
        dop[2 * i + 1] = std::max(dop[2 * i + 1], p[i]);
    }
}

template<typename T>
constexpr void Merge(Dop14<T>& dop, Dop14<T> const& other)
{
    for (u32 i { 0 }; i < cSlabs; i++) {
        dop[2 * i]     = std::min(dop[2 * i], other[2 * i]);
        dop[2 * i + 1] = std::max(dop[2 * i + 1], other[2 * i + 1]);
    }
}

// Init on the first point, Fit on the rest
template<typename T>
[[nodiscard]] constexpr Dop14<T> Fit(std::vector<Vec3<T>> const& points)
{
    assert(!points.empty());
    Dop14<T> dop { Init(points[0]) };
    for (u64 i { 1 }; i < points.size(); i++)
        Fit(dop, points[i]);
    return dop;
}

// exact: every f32 is an f64
[[nodiscard]] constexpr Dop14<f64> Widen(Dop14<f32> const& dop)
{
    Dop14<f64> result {};
    for (u32 i { 0 }; i < 2 * cSlabs; i++)
        result[i] = static_cast<f64>(dop[i]);
    return result;
}
