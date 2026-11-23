// Corner cutting and face clipping in f32 against cddlib in floating-point
// (ddf_, double, inside cddgmp) and rational mode, and the brute-force
// reference, all on the f32 fits; cdd and the references read them widened.

#include "timing.h"
#include "cases.h"
#include "corner_cutting.h"
#include "face_clipping.h"
#include "refs.h"

#include <algorithm>
#include <bit>
#include <chrono>
#include <cmath>
#include <numbers>
#include <print>
#include <string>
#include <string_view>
#include <vector>

#include <gmp.h>

extern "C" {
#include <setoper.h>

#include <cdd.h>
}

namespace {

constexpr std::array cTimed { Gen::eTriangle, Gen::eMerged, Gen::eCloud, Gen::eBox };
constexpr u64        cCases { 500 };

// One DOP per iteration: a loop the compiler vectorizes across DOPs times a
// batch instead, and whether it does depends on the compiler.
#if defined(__clang__)
#    define DOP14_NO_VECTORIZE _Pragma("clang loop vectorize(disable) interleave(disable)")
#else
#    define DOP14_NO_VECTORIZE _Pragma("GCC novector")
#endif

// ns per call, repeating the set for at least minSeconds
template<typename D, typename F>
[[nodiscard]] f64 nsPerCall(std::vector<D> const& dops, F&& f, f64 minSeconds)
{
    auto const start { std::chrono::steady_clock::now() };
    u64        calls { 0 };
    f64        elapsed { 0. };
    do {
        DOP14_NO_VECTORIZE
        for (D const& d : dops)
            f(d);
        calls += dops.size();
        elapsed = std::chrono::duration<f64>(std::chrono::steady_clock::now() - start).count();
    } while (elapsed < minSeconds);
    return 1e9 * elapsed / static_cast<f64>(calls);
}

// results go into an integer: a floating-point accumulator stays NaN once
// float cdd returns one, and each later x87 add on a NaN takes a slow
// microcode path that would inflate every time measured after it
void consume(u64& sink, f64 v)
{
    sink ^= std::bit_cast<u64>(v);
}

[[nodiscard]] std::string duration(f64 ns)
{
    return ns < 1e3 ? std::format("{:.1f} ns", ns) : std::format("{:.0f} us", ns / 1e3);
}

// the same half-spaces as RefCdd, in double
[[nodiscard]] ddf_MatrixPtr floatMatrix(Dop14<f64> const& dop)
{
    ddf_MatrixPtr m { ddf_CreateMatrix(14, 4) };
    m->representation = ddf_Inequality;
    m->numbtype       = ddf_Real;
    for (u32 i { 0 }; i < cSlabs; i++) {
        Vec3<i32> const n { cNormals[i] };
        myfloat* const  lo { m->matrix[2 * i] };
        myfloat* const  hi { m->matrix[2 * i + 1] };
        lo[0][0] = -dop[2 * i];
        hi[0][0] = dop[2 * i + 1];
        lo[1][0] = n.x;
        lo[2][0] = n.y;
        lo[3][0] = n.z;
        hi[1][0] = -n.x;
        hi[2][0] = -n.y;
        hi[3][0] = -n.z;
    }
    return m;
}

[[nodiscard]] u64 floatVertices(Dop14<f64> const& dop)
{
    ddf_MatrixPtr    m { floatMatrix(dop) };
    ddf_ErrorType    err { ddf_NoError };
    ddf_PolyhedraPtr p { ddf_DDMatrix2Poly(m, &err) };
    u64              count { 0 };
    if (p) {
        if (ddf_MatrixPtr g { ddf_CopyGenerators(p) }) {
            count = static_cast<u64>(g->rowsize);
            ddf_FreeMatrix(g);
        }
        ddf_FreePolyhedra(p);
    }
    ddf_FreeMatrix(m);
    return count;
}

// RefCdd's assembly in double; ok is false where cdd fails
[[nodiscard]] f64 floatArea(Dop14<f64> const& dop, bool& ok)
{
    ddf_MatrixPtr    m { floatMatrix(dop) };
    ddf_ErrorType    err { ddf_NoError };
    ddf_PolyhedraPtr p { ddf_DDMatrix2Poly(m, &err) };
    ddf_MatrixPtr    g { p ? ddf_CopyGenerators(p) : nullptr };
    ddf_SetFamilyPtr inc { p ? ddf_CopyIncidence(p) : nullptr };
    ok = err == ddf_NoError && g && inc && g->rowsize > 0;
    f64 area { 0. };
    if (ok) {
        std::vector<std::array<f64, 3>> v(static_cast<u64>(g->rowsize));
        for (u64 k { 0 }; k < v.size(); k++) {
            f64 const w { g->matrix[k][0][0] };
            ok &= w > 0.;
            for (u32 c { 0 }; c < 3; c++)
                v[k][c] = g->matrix[k][c + 1][0] / w;
        }
        std::vector<u64> ids;
        for (u32 r { 0 }; r < 14; r++) {
            ids.clear();
            for (u64 k { 0 }; k < v.size(); k++) {
                if (set_member(static_cast<long>(r + 1), inc->set[k]))
                    ids.push_back(k);
            }
            if (ids.size() < 3)
                continue;
            u32 const slab { r / 2 };
            u32 const drop { slab < 3 ? slab : 2 };
            u32 const a { (drop + 1) % 3 };
            u32 const b { (drop + 2) % 3 };
            f64       ca { 0. };
            f64       cb { 0. };
            for (u64 const id : ids) {
                ca += v[id][a];
                cb += v[id][b];
            }
            ca /= static_cast<f64>(ids.size());
            cb /= static_cast<f64>(ids.size());
            std::sort(ids.begin(), ids.end(), [&](u64 i, u64 j) { return std::atan2(v[i][b] - cb, v[i][a] - ca) < std::atan2(v[j][b] - cb, v[j][a] - ca); });
            f64 s { 0. };
            for (u64 n { 0 }; n < ids.size(); n++) {
                auto const& P { v[ids[n]] };
                auto const& Q { v[ids[(n + 1) % ids.size()]] };
                s += (P[a] - ca) * (Q[b] - cb) - (Q[a] - ca) * (P[b] - cb);
            }
            area += std::abs(s) / 2. * (slab < 3 ? 1. : std::numbers::sqrt3);
        }
    }
    if (inc)
        ddf_FreeSetFamily(inc);
    if (g)
        ddf_FreeMatrix(g);
    if (p)
        ddf_FreePolyhedra(p);
    ddf_FreeMatrix(m);
    return area;
}

[[nodiscard]] u64 rationalVertices(Dop14<f64> const& dop)
{
    dd_MatrixPtr m { dd_CreateMatrix(14, 4) };
    m->representation = dd_Inequality;
    m->numbtype       = dd_Rational;
    for (u32 i { 0 }; i < cSlabs; i++) {
        Vec3<i32> const n { cNormals[i] };
        mytype* const   lo { m->matrix[2 * i] };
        mytype* const   hi { m->matrix[2 * i + 1] };
        mpq_set_d(lo[0], -dop[2 * i]);
        mpq_set_d(hi[0], dop[2 * i + 1]);
        mpq_set_si(lo[1], n.x, 1);
        mpq_set_si(lo[2], n.y, 1);
        mpq_set_si(lo[3], n.z, 1);
        mpq_set_si(hi[1], -n.x, 1);
        mpq_set_si(hi[2], -n.y, 1);
        mpq_set_si(hi[3], -n.z, 1);
    }
    dd_ErrorType    err { dd_NoError };
    dd_PolyhedraPtr p { dd_DDMatrix2Poly(m, &err) };
    u64             count { 0 };
    if (p) {
        if (dd_MatrixPtr g { dd_CopyGenerators(p) }) {
            count = static_cast<u64>(g->rowsize);
            dd_FreeMatrix(g);
        }
        dd_FreePolyhedra(p);
    }
    dd_FreeMatrix(m);
    return count;
}

[[nodiscard]] std::vector<Dop14<f32>> casesOf(Gen gen, u64 seed, u64 count)
{
    std::vector<Dop14<f32>> dops;
    for (u64 i { 0 }; i < count; i++)
        dops.push_back(Generate(gen, seed, i).fit32);
    return dops;
}

}

void PrintTiming(u64 seed)
{
    ddf_set_global_constants();
    // RefCdd sets the rational globals on its first call
    static_cast<void>(RefCdd(Widen(Generate(Gen::eBox, seed, 0).fit32)));

    struct Row {
        std::string_view name;
        f64              minSeconds;
        std::vector<f64> ns;
    };
    std::vector<Row> rows {
        { .name = "corner cutting  CornerCut<f32>", .minSeconds = .5, .ns = {} },
        { .name = "face clipping   FaceClip<f32>", .minSeconds = .2, .ns = {} },
        { .name = "cdd float       vertex enumeration", .minSeconds = .2, .ns = {} },
        { .name = "cdd float       enumeration and area", .minSeconds = .2, .ns = {} },
        { .name = "cdd rational    vertex enumeration", .minSeconds = .2, .ns = {} },
        { .name = "cdd rational    RefCdd", .minSeconds = .2, .ns = {} },
        { .name = "brute force     RefBrute", .minSeconds = .2, .ns = {} },
    };
    u64 sink { 0 };
    for (Gen const gen : cTimed) {
        auto const              dops { casesOf(gen, seed, cCases) };
        std::vector<Dop14<f64>> wide;
        for (Dop14<f32> const& d : dops)
            wide.push_back(Widen(d));
        bool ok { true };
        rows[0].ns.push_back(nsPerCall(dops, [&](Dop14<f32> const& d) { sink ^= std::bit_cast<u32>(CornerCut(d)); }, rows[0].minSeconds));
        rows[1].ns.push_back(nsPerCall(dops, [&](Dop14<f32> const& d) { sink ^= std::bit_cast<u32>(FaceClip(d)); }, rows[1].minSeconds));
        rows[2].ns.push_back(nsPerCall(wide, [&](Dop14<f64> const& d) { sink += floatVertices(d); }, rows[2].minSeconds));
        rows[3].ns.push_back(nsPerCall(wide, [&](Dop14<f64> const& d) { consume(sink, floatArea(d, ok)); }, rows[3].minSeconds));
        rows[4].ns.push_back(nsPerCall(wide, [&](Dop14<f64> const& d) { sink += rationalVertices(d); }, rows[4].minSeconds));
        rows[5].ns.push_back(nsPerCall(wide, [&](Dop14<f64> const& d) { consume(sink, static_cast<f64>(RefCdd(d).area)); }, rows[5].minSeconds));
        rows[6].ns.push_back(nsPerCall(wide, [&](Dop14<f64> const& d) { consume(sink, static_cast<f64>(RefBrute(d).area)); }, rows[6].minSeconds));
    }

    std::println("\ntime per DOP, one thread, {} cases per generator\n", cCases);
    std::string header { std::format("{:<40}", "method") };
    for (Gen const gen : cTimed)
        header += std::format("{:>12}", strGen[static_cast<u32>(gen)]);
    std::println("{}", header);
    for (Row const& row : rows) {
        std::string line { std::format("{:<40}", row.name) };
        for (f64 const ns : row.ns)
            line += std::format("{:>12}", duration(ns));
        std::println("{}", line);
    }
    for (auto const& [name, row] : { std::pair { "CornerCut<f32>", 0 }, std::pair { "FaceClip<f32>", 1 } }) {
        std::string line { std::format("{:<40}", std::format("cdd float enumeration / {}", name)) };
        for (u64 g { 0 }; g < cTimed.size(); g++)
            line += std::format("{:>12}", std::format("{:.0f}x", rows[2].ns[g] / rows[row].ns[g]));
        std::println("{}", line);
    }

    ddf_free_global_constants();
    // a volatile store keeps every timed result alive
    [[maybe_unused]] static u64 volatile kept;
    kept = sink;
}
