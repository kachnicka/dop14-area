#pragma once

// Face clipping (Kacerik and Bittner 2024, section 3.2): each of the 14
// planes starts from its section of the box and is clipped by the 12 planes
// not parallel to it; the faces' areas sum to the surface area. Exact up to
// rounding.

#include "dop14.h"

#include <array>
#include <cassert>
#include <cmath>
#include <numbers>
#include <utility>

namespace face_clipping {

template<typename T>
[[nodiscard]] constexpr T& at(Vec3<T>& v, u32 axis)
{
    return axis == 0 ? v.x : (axis == 1 ? v.y : v.z);
}

template<typename T>
[[nodiscard]] constexpr T at(Vec3<T> const& v, u32 axis)
{
    return axis == 0 ? v.x : (axis == 1 ? v.y : v.z);
}

// 4 to start, at most one more per clip; rounding can break convexity, hence the margin
template<typename T>
struct Polygon {
    std::array<Vec3<T>, 32> v;
    u32                     count;
};

// keeps n . x >= bound (min side) or n . x <= bound (max side)
template<typename T>
void clip(Polygon<T> const& in, Polygon<T>& out, Vec3<i32> n, T bound, bool maxSide)
{
    auto const outside { [&](Vec3<T> const& v) {
        T const d { static_cast<T>(n.x) * v.x + static_cast<T>(n.y) * v.y + static_cast<T>(n.z) * v.z };
        return maxSide ? d - bound : bound - d;
    } };
    out.count = 0;
    for (u32 i { 0 }; i < in.count; i++) {
        Vec3<T> const& p { in.v[i] };
        Vec3<T> const& q { in.v[(i + 1) % in.count] };
        T const        dp { outside(p) };
        T const        dq { outside(q) };
        if (dp <= T(0) && out.count < out.v.size())
            out.v[out.count++] = p;
        if (((dp < T(0) && dq > T(0)) || (dp > T(0) && dq < T(0))) && out.count < out.v.size()) {
            T const t { dp / (dp - dq) };
            out.v[out.count++] = { .x = p.x + t * (q.x - p.x), .y = p.y + t * (q.y - p.y), .z = p.z + t * (q.z - p.z) };
        }
    }
    assert(out.count < out.v.size());
}

}

template<typename T>
[[nodiscard]] T FaceClip(Dop14<T> const& slab)
{
    using namespace face_clipping;
    for (u32 i { 0 }; i < cSlabs; i++) {
        if (slab[2 * i] > slab[2 * i + 1])
            return T(0);
    }

    // twice the face areas: axis faces as they are, diagonal faces projected
    // onto xy, sqrt(3) too small
    T          twiceAxis { 0 };
    T          twiceDiag { 0 };
    Polygon<T> buffers[2] {};
    for (u32 f { 0 }; f < 14; f++) {
        u32 const       s { f / 2 };
        Vec3<i32> const n { cNormals[s] };
        // the plane n . x = slab[f], solved for axis k over the box's
        // bounds in a and b; n_k is +-1
        u32 const   k { s < 3 ? s : 2 };
        u32 const   a { (k + 1) % 3 };
        u32 const   b { (k + 2) % 3 };
        T const     na { static_cast<T>(at(n, a)) };
        T const     nb { static_cast<T>(at(n, b)) };
        T const     nk { static_cast<T>(at(n, k)) };
        T const     ua[4] { slab[2 * a], slab[2 * a + 1], slab[2 * a + 1], slab[2 * a] };
        T const     ub[4] { slab[2 * b], slab[2 * b], slab[2 * b + 1], slab[2 * b + 1] };
        Polygon<T>* p { &buffers[0] };
        Polygon<T>* q { &buffers[1] };
        for (u32 j { 0 }; j < 4; j++) {
            Vec3<T> v {};
            at(v, a) = ua[j];
            at(v, b) = ub[j];
            at(v, k) = (slab[f] - na * ua[j] - nb * ub[j]) * nk;
            p->v[j]  = v;
        }
        p->count = 4;

        for (u32 g { 0 }; g < 14 && p->count > 0; g++) {
            if (g / 2 == s)
                continue;
            clip(*p, *q, cNormals[g / 2], slab[g], (g & 1) != 0);
            std::swap(p, q);
        }
        if (p->count < 3)
            continue;

        T sum { 0 };
        for (u32 j { 1 }; j + 1 < p->count; j++) {
            T const ea { at(p->v[j], a) - at(p->v[0], a) };
            T const eb { at(p->v[j], b) - at(p->v[0], b) };
            T const fa { at(p->v[j + 1], a) - at(p->v[0], a) };
            T const fb { at(p->v[j + 1], b) - at(p->v[0], b) };
            sum += ea * fb - eb * fa;
        }
        (s < 3 ? twiceAxis : twiceDiag) += std::abs(sum);
    }
    return (twiceAxis + std::numbers::sqrt3_v<T> * twiceDiag) / T(2);
}
