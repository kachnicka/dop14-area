#include "exact.h"
#include "refs.h"

#include <algorithm>
#include <memory>
#include <numeric>
#include <vector>

// exact.h brings gmp.h, which declares C++ operators, ahead of the extern "C"
// that setoper.h (no C++ guard of its own) needs
extern "C" {
#include <setoper.h>

#include <cdd.h>
}

#ifndef GMPRATIONAL
#    error "ref_cdd needs cddlib's GMP rational build (cddgmp, GMPRATIONAL)"
#endif

namespace {

using exact::Q;

// required around every cddlib call
struct CddGlobals {
    CddGlobals() { dd_set_global_constants(); }
    ~CddGlobals() { dd_free_global_constants(); }
    CddGlobals(CddGlobals const&)            = delete;
    CddGlobals& operator=(CddGlobals const&) = delete;
};

struct FreeMatrix {
    void operator()(dd_MatrixType* m) const { dd_FreeMatrix(m); }
};
struct FreePolyhedra {
    void operator()(dd_PolyhedraType* p) const { dd_FreePolyhedra(p); }
};
struct FreeSetFamily {
    void operator()(dd_SetFamilyType* f) const { dd_FreeSetFamily(f); }
};
using Matrix    = std::unique_ptr<dd_MatrixType, FreeMatrix>;
using Polyhedra = std::unique_ptr<dd_PolyhedraType, FreePolyhedra>;
using SetFamily = std::unique_ptr<dd_SetFamilyType, FreeSetFamily>;

using Vertex = std::array<Q, 3>;

// Exact twice-area of a convex polygon projected along axis `drop`: angular
// sort around the centroid with exact predicates, then shoelace.
void twiceProjectedArea(std::vector<Vertex> const& vertices, std::vector<u32> const& ids, u32 drop, Q& result)
{
    u32 const a { (drop + 1) % 3 };
    u32 const b { (drop + 2) % 3 };
    u64 const k { ids.size() };

    Q ca;
    Q cb;
    for (u32 const id : ids) {
        mpq_add(ca.v, ca.v, vertices[id][a].v);
        mpq_add(cb.v, cb.v, vertices[id][b].v);
    }
    Q count;
    mpq_set_ui(count.v, static_cast<unsigned long>(k), 1);
    mpq_div(ca.v, ca.v, count.v);
    mpq_div(cb.v, cb.v, count.v);

    std::vector<Q> da(k);
    std::vector<Q> db(k);
    for (u64 i { 0 }; i < k; i++) {
        mpq_sub(da[i].v, vertices[ids[i]][a].v, ca.v);
        mpq_sub(db[i].v, vertices[ids[i]][b].v, cb.v);
    }

    Q lhs;
    Q rhs;
    // half-plane first, then the turn
    auto const       half { [&](u64 i) { return (mpq_sgn(db[i].v) > 0 || (mpq_sgn(db[i].v) == 0 && mpq_sgn(da[i].v) > 0)) ? 0 : 1; } };
    std::vector<u64> order(k);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](u64 i, u64 j) {
        i32 const hi { half(i) };
        i32 const hj { half(j) };
        if (hi != hj)
            return hi < hj;
        mpq_mul(lhs.v, da[i].v, db[j].v);
        mpq_mul(rhs.v, db[i].v, da[j].v);
        return mpq_cmp(lhs.v, rhs.v) > 0;
    });

    mpq_set_ui(result.v, 0, 1);
    for (u64 n { 0 }; n < k; n++) {
        u64 const i { order[n] };
        u64 const j { order[(n + 1) % k] };
        mpq_mul(lhs.v, da[i].v, db[j].v);
        mpq_mul(rhs.v, db[i].v, da[j].v);
        mpq_sub(lhs.v, lhs.v, rhs.v);
        mpq_add(result.v, result.v, lhs.v);
    }
    mpq_abs(result.v, result.v);
}

}

RefArea RefCdd(Dop14<f64> const& dop)
{
    // once per process, freed at exit
    static CddGlobals const globals;

    // row [b, a] is b + a . x >= 0: n . x >= min, then n . x <= max
    Matrix const m { dd_CreateMatrix(14, 4) };
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

    RefArea         result { .area = 0.L, .vertices = 0, .ok = false };
    dd_ErrorType    err { dd_NoError };
    Polyhedra const poly { dd_DDMatrix2Poly(m.get(), &err) };
    if (err != dd_NoError || !poly)
        return result;
    Matrix const    generators { dd_CopyGenerators(poly.get()) };
    SetFamily const incidence { dd_CopyIncidence(poly.get()) };
    if (!generators || !incidence || incidence->famsize != generators->rowsize)
        return result;

    // generator rows: [1, x, y, z] a vertex, [0, ...] a ray or line
    std::vector<Vertex> vertices(static_cast<u64>(generators->rowsize));
    for (u64 k { 0 }; k < vertices.size(); k++) {
        mytype* const row { generators->matrix[k] };
        if (mpq_sgn(row[0]) <= 0 || set_member(static_cast<long>(k + 1), generators->linset))
            return result;
        for (u32 c { 0 }; c < 3; c++)
            mpq_div(vertices[k][c].v, row[c + 1], row[0]);
    }
    result.ok       = true;
    result.vertices = static_cast<u32>(vertices.size());
    if (vertices.empty())
        return result;

    // twice the areas, diagonal facets projected along z
    Q                axis2;
    Q                diag2;
    Q                facet2;
    std::vector<u32> ids;
    for (u32 r { 0 }; r < 14; r++) {
        ids.clear();
        for (u64 k { 0 }; k < vertices.size(); k++) {
            if (set_member(static_cast<long>(r + 1), incidence->set[k]))
                ids.push_back(static_cast<u32>(k));
        }
        if (ids.size() < 3)
            continue;
        u32 const slab { r / 2 };
        twiceProjectedArea(vertices, ids, slab < 3 ? slab : 2, facet2);
        mpq_add(slab < 3 ? axis2.v : diag2.v, slab < 3 ? axis2.v : diag2.v, facet2.v);
    }
    result.area = exact::Area(axis2.v, diag2.v);
    return result;
}
