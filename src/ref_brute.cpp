#include "exact.h"
#include "refs.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <vector>

namespace {

using exact::Q;
using exact::Z;

// planes a, b, c from distinct slabs meet at x with det x = beta_a (b x c) +
// beta_b (c x a) + beta_c (a x b), det = a . (b x c). Every nonzero det of a
// 3x3 matrix over {-1, 0, 1} divides 12, so 12 x has integer weights.
struct Triple {
    std::array<u32, 3>       slab;
    std::array<Vec3<i32>, 3> weight;
};

struct Triples {
    std::array<Triple, 35> list;
    u32                    count;
    bool                   integral;
};

[[nodiscard]] consteval Vec3<i32> cross(Vec3<i32> a, Vec3<i32> b)
{
    return { .x = a.y * b.z - a.z * b.y, .y = a.z * b.x - a.x * b.z, .z = a.x * b.y - a.y * b.x };
}

[[nodiscard]] consteval Triples solvable()
{
    Triples t { .list = {}, .count = 0, .integral = true };
    for (u32 i { 0 }; i < cSlabs; i++) {
        for (u32 j { i + 1 }; j < cSlabs; j++) {
            for (u32 k { j + 1 }; k < cSlabs; k++) {
                Vec3<i32> const a { cNormals[i] };
                Vec3<i32> const b { cNormals[j] };
                Vec3<i32> const c { cNormals[k] };
                Vec3<i32> const bc { cross(b, c) };
                i32 const       det { a.x * bc.x + a.y * bc.y + a.z * bc.z };
                if (det == 0)
                    continue;
                t.integral = t.integral && 12 % det == 0;
                i32 const  scale { 12 / det };
                auto const scaled { [scale](Vec3<i32> v) { return Vec3<i32> { .x = scale * v.x, .y = scale * v.y, .z = scale * v.z }; } };
                t.list[t.count++] = { .slab = { i, j, k }, .weight = { scaled(bc), scaled(cross(c, a)), scaled(cross(a, b)) } };
            }
        }
    }
    return t;
}

constexpr Triples cTriples { solvable() };
static_assert(cTriples.integral);

[[nodiscard]] i32 component(Vec3<i32> v, u32 axis)
{
    return axis == 0 ? v.x : (axis == 1 ? v.y : v.z);
}

// acc += w * b for a small integer w
void addMul(mpz_t acc, mpz_srcptr b, i32 w)
{
    if (w > 0)
        mpz_addmul_ui(acc, b, static_cast<unsigned long>(w));
    else if (w < 0)
        mpz_submul_ui(acc, b, static_cast<unsigned long>(-w));
}

// 12 x in the slabs' integer units, and its seven projections n . 12 x
struct Vertex {
    std::array<Z, 3>      x;
    std::array<Z, cSlabs> n;
};

}

RefArea RefBrute(Dop14<f64> const& dop)
{
    // every slab value is m 2^e with an integer m; in units of 2^e for the
    // smallest e, all are integers
    std::array<i64, 14> mantissa {};
    std::array<i32, 14> exponent {};
    i32                 unitExp { 0 };
    bool                any { false };
    for (u32 p { 0 }; p < 14; p++) {
        assert(std::isfinite(dop[p]));
        i32       e { 0 };
        f64 const f { std::frexp(dop[p], &e) };
        mantissa[p] = static_cast<i64>(std::ldexp(f, 53));
        exponent[p] = e - 53;
        if (mantissa[p] != 0) {
            unitExp = any ? std::min(unitExp, exponent[p]) : exponent[p];
            any     = true;
        }
    }
    // the planes, times 12
    std::array<Z, 14> beta;
    std::array<Z, 14> beta12;
    for (u32 p { 0 }; p < 14; p++) {
        mpz_set_si(beta[p].v, mantissa[p]);
        if (mantissa[p] != 0)
            mpz_mul_2exp(beta[p].v, beta[p].v, static_cast<mp_bitcnt_t>(exponent[p] - unitExp));
        mpz_mul_ui(beta12[p].v, beta[p].v, 12);
    }

    std::vector<Vertex> vertices;
    Vertex              v;
    for (u32 t { 0 }; t < cTriples.count; t++) {
        Triple const& triple { cTriples.list[t] };
        for (u32 sides { 0 }; sides < 8; sides++) {
            std::array<u32, 3> plane {};
            for (u32 j { 0 }; j < 3; j++)
                plane[j] = 2 * triple.slab[j] + ((sides >> j) & 1);
            for (u32 axis { 0 }; axis < 3; axis++) {
                mpz_set_ui(v.x[axis].v, 0);
                for (u32 j { 0 }; j < 3; j++)
                    addMul(v.x[axis].v, beta[plane[j]].v, component(triple.weight[j], axis));
            }
            bool inside { true };
            for (u32 s { 0 }; s < cSlabs && inside; s++) {
                mpz_set_ui(v.n[s].v, 0);
                for (u32 axis { 0 }; axis < 3; axis++)
                    addMul(v.n[s].v, v.x[axis].v, component(cNormals[s], axis));
                inside = mpz_cmp(v.n[s].v, beta12[2 * s].v) >= 0 && mpz_cmp(v.n[s].v, beta12[2 * s + 1].v) <= 0;
            }
            if (!inside)
                continue;
            bool const known { std::any_of(vertices.begin(), vertices.end(), [&v](Vertex const& u) {
                return mpz_cmp(u.x[0].v, v.x[0].v) == 0 && mpz_cmp(u.x[1].v, v.x[1].v) == 0 && mpz_cmp(u.x[2].v, v.x[2].v) == 0;
            }) };
            if (!known)
                vertices.push_back(std::move(v));
        }
    }

    // twice the face areas, diagonal faces projected along z: the polygon's
    // vertices are extreme points, so an angular sort around the
    // lexicographically smallest one is exact
    Z                twiceAxis;
    Z                twiceDiag;
    Z                term;
    Z                polygon;
    std::vector<u32> ids;
    std::vector<Z>   da;
    std::vector<Z>   db;
    for (u32 p { 0 }; p < 14; p++) {
        u32 const s { p / 2 };
        ids.clear();
        for (u32 i { 0 }; i < vertices.size(); i++) {
            if (mpz_cmp(vertices[i].n[s].v, beta12[p].v) == 0)
                ids.push_back(i);
        }
        if (ids.size() < 3)
            continue;
        u32 const drop { s < 3 ? s : 2 };
        u32 const a { (drop + 1) % 3 };
        u32 const b { (drop + 2) % 3 };
        auto      lexLess { [&](u32 i, u32 j) {
            i32 const ca { mpz_cmp(vertices[i].x[a].v, vertices[j].x[a].v) };
            return ca < 0 || (ca == 0 && mpz_cmp(vertices[i].x[b].v, vertices[j].x[b].v) < 0);
        } };
        std::iter_swap(ids.begin(), std::min_element(ids.begin(), ids.end(), lexLess));
        da.resize(ids.size());
        db.resize(ids.size());
        for (u32 i { 1 }; i < ids.size(); i++) {
            mpz_sub(da[i].v, vertices[ids[i]].x[a].v, vertices[ids[0]].x[a].v);
            mpz_sub(db[i].v, vertices[ids[i]].x[b].v, vertices[ids[0]].x[b].v);
        }
        std::vector<u32> order(ids.size() - 1);
        for (u32 i { 0 }; i < order.size(); i++)
            order[i] = i + 1;
        // the cross product of the edges to i and j: positive when j turns left
        auto const turn { [&](u32 i, u32 j) {
            mpz_mul(term.v, da[i].v, db[j].v);
            mpz_submul(term.v, db[i].v, da[j].v);
            return mpz_sgn(term.v);
        } };
        std::sort(order.begin(), order.end(), [&](u32 i, u32 j) { return turn(i, j) > 0; });
        mpz_set_ui(polygon.v, 0);
        for (u32 k { 0 }; k + 1 < order.size(); k++) {
            mpz_mul(term.v, da[order[k]].v, db[order[k + 1]].v);
            mpz_submul(term.v, db[order[k]].v, da[order[k + 1]].v);
            mpz_add(polygon.v, polygon.v, term.v);
        }
        mpz_add(s < 3 ? twiceAxis.v : twiceDiag.v, s < 3 ? twiceAxis.v : twiceDiag.v, polygon.v);
    }

    // 12 x counts units of 2^unitExp, so the twice-areas count 4^unitExp / 144
    Q axis;
    Q diag;
    for (auto [from, to] : { std::pair { &twiceAxis, &axis }, std::pair { &twiceDiag, &diag } }) {
        mpq_set_z(to->v, from->v);
        if (unitExp >= 0)
            mpq_mul_2exp(to->v, to->v, static_cast<mp_bitcnt_t>(2 * unitExp));
        else
            mpq_div_2exp(to->v, to->v, static_cast<mp_bitcnt_t>(-2 * unitExp));
        mpz_mul_ui(mpq_denref(to->v), mpq_denref(to->v), 144);
        mpq_canonicalize(to->v);
    }
    return { .area = exact::Area(axis.v, diag.v), .vertices = static_cast<u32>(vertices.size()), .ok = true };
}
