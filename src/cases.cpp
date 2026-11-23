#include "cases.h"

#include <cmath>
#include <limits>
#include <numbers>

std::vector<CuratedCase> CuratedCases()
{
    std::vector<CuratedCase> cases;
    auto                     add { [&cases](std::string_view name, std::vector<Vec3<f32>> const& points, fext closedForm) {
        cases.push_back({ .name = name, .dop = Fit(points), .closedForm = closedForm });
    } };
    fext const sqrt3 { std::numbers::sqrt3_v<fext> };

    add("unit cube",
        {
            { .x = 0.f, .y = 0.f, .z = 0.f },
            { .x = 1.f, .y = 0.f, .z = 0.f },
            { .x = 0.f, .y = 1.f, .z = 0.f },
            { .x = 1.f, .y = 1.f, .z = 0.f },
            { .x = 0.f, .y = 0.f, .z = 1.f },
            { .x = 1.f, .y = 0.f, .z = 1.f },
            { .x = 0.f, .y = 1.f, .z = 1.f },
            { .x = 1.f, .y = 1.f, .z = 1.f },
        },
        6.L);

    // edge midpoints of [-1, 1]^3: eight depth-1 cuts that never meet
    std::vector<Vec3<f32>> cuboctahedron;
    for (f32 sx : { -1.f, 1.f }) {
        for (f32 sy : { -1.f, 1.f }) {
            for (f32 sz : { -1.f, 1.f }) {
                cuboctahedron.push_back({ .x = sx, .y = sy, .z = 0.f });
                cuboctahedron.push_back({ .x = sx, .y = 0.f, .z = sz });
                cuboctahedron.push_back({ .x = 0.f, .y = sy, .z = sz });
            }
        }
    }
    add("cuboctahedron", cuboctahedron, 12.L + 4.L * sqrt3);

    // alternating cube corners: four cuts of depth 2 meet at the edge midpoints
    add("regular tetrahedron",
        {
            { .x = 1.f, .y = 1.f, .z = 1.f },
            { .x = 1.f, .y = -1.f, .z = -1.f },
            { .x = -1.f, .y = 1.f, .z = -1.f },
            { .x = -1.f, .y = -1.f, .z = 1.f },
        },
        8.L * sqrt3);

    // both sides count
    add("flat unit square",
        {
            { .x = 0.f, .y = 0.f, .z = 0.f },
            { .x = 1.f, .y = 0.f, .z = 0.f },
            { .x = 0.f, .y = 0.f, .z = 1.f },
            { .x = 1.f, .y = 0.f, .z = 1.f },
        },
        2.L);

    add("thin diagonal box",
        {
            { .x = 0.f, .y = 0.f, .z = -.05f },
            { .x = 10.f, .y = 10.f, .z = -.05f },
            { .x = 0.f, .y = 0.f, .z = .05f },
            { .x = 10.f, .y = 10.f, .z = .05f },
            { .x = 0.f, .y = .2f, .z = -.05f },
            { .x = 10.f, .y = 10.2f, .z = -.05f },
            { .x = 0.f, .y = .2f, .z = .05f },
            { .x = 10.f, .y = 10.2f, .z = .05f },
        },
        std::numeric_limits<fext>::quiet_NaN());

    return cases;
}

namespace {

using V = Vec3<f64>;

class Rng {
    u64 state;

    [[nodiscard]] static u64 splitMix64(u64& s)
    {
        u64 z { s += 0x9e3779b97f4a7c15ull };
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ull;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebull;
        return z ^ (z >> 31);
    }

public:
    Rng(u64 seed, Gen gen, u64 index)
    {
        u64 s { seed };
        s     = splitMix64(s) ^ static_cast<u64>(gen);
        s     = splitMix64(s) ^ index;
        state = splitMix64(s);
    }

    [[nodiscard]] u64 next() { return splitMix64(state); }
    // [0, 1) from the top 53 bits
    [[nodiscard]] f64 uniform() { return static_cast<f64>(next() >> 11) * 0x1p-53; }
    [[nodiscard]] f64 uniform(f64 lo, f64 hi) { return lo + (hi - lo) * uniform(); }
    [[nodiscard]] f64 logUniform(f64 lo, f64 hi) { return lo * std::pow(hi / lo, uniform()); }
    [[nodiscard]] u32 below(u32 n) { return static_cast<u32>(next() % n); }
};

// triangles: Init of three per triangle, then merged; otherwise one point set
struct Shape {
    std::vector<V> points;
    bool           triangles { false };
};

[[nodiscard]] V extents(Rng& rng)
{
    return { .x = rng.logUniform(1e-2, 10.), .y = rng.logUniform(1e-2, 10.), .z = rng.logUniform(1e-2, 10.) };
}

[[nodiscard]] V center(Rng& rng)
{
    return { .x = rng.uniform(-1., 1.), .y = rng.uniform(-1., 1.), .z = rng.uniform(-1., 1.) };
}

[[nodiscard]] V uniformIn(Rng& rng, V c, V e)
{
    return { .x = c.x + (rng.uniform() - .5) * e.x, .y = c.y + (rng.uniform() - .5) * e.y, .z = c.z + (rng.uniform() - .5) * e.z };
}

[[nodiscard]] f64& component(V& v, u32 axis)
{
    return axis == 0 ? v.x : (axis == 1 ? v.y : v.z);
}

[[nodiscard]] f64 component(V const& v, u32 axis)
{
    return axis == 0 ? v.x : (axis == 1 ? v.y : v.z);
}

[[nodiscard]] Shape pointsInBox(Rng& rng, V e, u32 count)
{
    V const c { center(rng) };
    Shape   shape;
    for (u32 i { 0 }; i < count; i++)
        shape.points.push_back(uniformIn(rng, c, e));
    return shape;
}

// points on one of the seven DOP planes; f32 rounding leaves diagonal ones a few ulps off
[[nodiscard]] Shape flat(Rng& rng)
{
    u32 const axis { rng.below(cSlabs) };
    V const   c { center(rng) };
    V const   e { extents(rng) };
    u32 const count { 3 + rng.below(8) };
    Shape     shape;
    for (u32 i { 0 }; i < count; i++) {
        V p { uniformIn(rng, c, e) };
        if (axis < 3) {
            component(p, axis) = component(c, axis);
        } else {
            V const   n { .x = static_cast<f64>(cNormals[axis].x), .y = static_cast<f64>(cNormals[axis].y), .z = static_cast<f64>(cNormals[axis].z) };
            f64 const k { n.x * c.x + n.y * c.y + n.z * c.z };
            p.z = (k - n.x * p.x - n.y * p.y) * n.z;
        }
        shape.points.push_back(p);
    }
    return shape;
}

// the corners of a randomly oriented box (uniform rotation, Shoemake)
[[nodiscard]] Shape box(Rng& rng)
{
    f64 const u1 { rng.uniform() };
    f64 const a2 { 2. * std::numbers::pi * rng.uniform() };
    f64 const a3 { 2. * std::numbers::pi * rng.uniform() };
    f64 const x { std::sqrt(1. - u1) * std::sin(a2) };
    f64 const y { std::sqrt(1. - u1) * std::cos(a2) };
    f64 const z { std::sqrt(u1) * std::sin(a3) };
    f64 const w { std::sqrt(u1) * std::cos(a3) };
    V const   r0 { .x = 1. - 2. * (y * y + z * z), .y = 2. * (x * y - z * w), .z = 2. * (x * z + y * w) };
    V const   r1 { .x = 2. * (x * y + z * w), .y = 1. - 2. * (x * x + z * z), .z = 2. * (y * z - x * w) };
    V const   r2 { .x = 2. * (x * z - y * w), .y = 2. * (y * z + x * w), .z = 1. - 2. * (x * x + y * y) };

    V const c { center(rng) };
    V const e { extents(rng) };
    Shape   shape;
    for (u32 k { 0 }; k < 8; k++) {
        V const h { .x = (k & 1 ? .5 : -.5) * e.x, .y = (k & 2 ? .5 : -.5) * e.y, .z = (k & 4 ? .5 : -.5) * e.z };
        shape.points.push_back({
            .x = c.x + r0.x * h.x + r0.y * h.y + r0.z * h.z,
            .y = c.y + r1.x * h.x + r1.y * h.y + r1.z * h.z,
            .z = c.z + r2.x * h.x + r2.y * h.y + r2.z * h.z,
        });
    }
    return shape;
}

// 2 to 7 triangles inside one node box, each 5% to 50% of the node's size
[[nodiscard]] Shape merged(Rng& rng)
{
    V const   c { center(rng) };
    V const   e { extents(rng) };
    u32 const count { 2 + rng.below(6) };
    Shape     shape { .points = {}, .triangles = true };
    for (u32 t { 0 }; t < count; t++) {
        V const   anchor { uniformIn(rng, c, e) };
        f64 const size { rng.uniform(.05, .5) };
        for (u32 k { 0 }; k < 3; k++)
            shape.points.push_back(uniformIn(rng, anchor, { .x = size * e.x, .y = size * e.y, .z = size * e.z }));
    }
    return shape;
}

[[nodiscard]] Shape base(Gen gen, Rng& rng)
{
    switch (gen) {
    case Gen::eCloud:
        return pointsInBox(rng, extents(rng), 4 + rng.below(12));
    case Gen::eSliver: {
        V e { extents(rng) };
        component(e, rng.below(3)) = rng.logUniform(1e-5, 1e-2);
        return pointsInBox(rng, e, 4 + rng.below(12));
    }
    case Gen::eTriangle: {
        Shape shape { pointsInBox(rng, extents(rng), 3) };
        shape.triangles = true;
        return shape;
    }
    case Gen::eSegment:
        return pointsInBox(rng, extents(rng), 2);
    case Gen::eFlat:
        return flat(rng);
    case Gen::eBox:
        return box(rng);
    case Gen::eMerged:
        return merged(rng);
    case Gen::eFar:
    case Gen::eScale:
        break;
    }
    assert(false);
    return {};
}

template<typename T>
[[nodiscard]] Dop14<T> fit(std::vector<Vec3<T>> const& p, bool triangles)
{
    if (!triangles)
        return Fit(p);
    assert(p.size() % 3 == 0);
    Dop14<T> dop { Init(p[0], p[1], p[2]) };
    for (u64 i { 3 }; i < p.size(); i += 3)
        Merge(dop, Init(p[i], p[i + 1], p[i + 2]));
    return dop;
}

// the f32 points widen exactly, so both fits bound the same points
[[nodiscard]] Fits fit(Shape const& shape)
{
    std::vector<Vec3<f32>> p32;
    std::vector<Vec3<f64>> p64;
    for (V const& v : shape.points) {
        Vec3<f32> const p { .x = static_cast<f32>(v.x), .y = static_cast<f32>(v.y), .z = static_cast<f32>(v.z) };
        p32.push_back(p);
        p64.push_back({ .x = p.x, .y = p.y, .z = p.z });
    }
    return { .fit32 = fit(p32, shape.triangles), .fit64 = fit(p64, shape.triangles) };
}

}

Fits Generate(Gen gen, u64 seed, u64 index)
{
    Rng rng { seed, gen, index };
    if (gen != Gen::eFar && gen != Gen::eScale)
        return fit(base(gen, rng));

    Gen const inner { static_cast<Gen>(rng.below(static_cast<u32>(Gen::eFar))) };
    Shape     shape { base(inner, rng) };
    if (gen == Gen::eFar) {
        // a uniform direction, 1e2 to 1e5 away
        f64 const r { rng.logUniform(1e2, 1e5) };
        f64 const z { rng.uniform(-1., 1.) };
        f64 const phi { 2. * std::numbers::pi * rng.uniform() };
        f64 const s { std::sqrt(1. - z * z) };
        V const   offset { .x = r * s * std::cos(phi), .y = r * s * std::sin(phi), .z = r * z };
        for (V& v : shape.points)
            v = { .x = v.x + offset.x, .y = v.y + offset.y, .z = v.z + offset.z };
    } else {
        f64 const f { rng.logUniform(1e-4, 1e6) };
        for (V& v : shape.points)
            v = { .x = f * v.x, .y = f * v.y, .z = f * v.z };
    }
    return fit(shape);
}
