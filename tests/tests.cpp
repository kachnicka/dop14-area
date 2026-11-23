#include "cases.h"
#include "corner_cutting.h"
#include "face_clipping.h"
#include "refs.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <print>
#include <string_view>
#include <vector>

namespace {

enum class Overlap : u8 {
    eNone,
    ePairwise,
    eCornerInSeveral,
};

// the corners each s[i] of CornerCut pairs
constexpr std::array<std::array<u32, 2>, 12> cOverlapPairs { {
    { 0, 7 },
    { 1, 6 },
    { 2, 5 },
    { 3, 4 },
    { 0, 4 },
    { 1, 5 },
    { 2, 6 },
    { 3, 7 },
    { 0, 2 },
    { 1, 3 },
    { 4, 6 },
    { 5, 7 },
} };

template<typename T>
[[nodiscard]] constexpr Overlap classify(AreaTerms<T> const& terms)
{
    std::array<u32, 8> perCorner {};
    bool               any { false };
    for (u32 i { 0 }; i < 12; i++) {
        if (terms.s[i] > T(0)) {
            any = true;
            perCorner[cOverlapPairs[i][0]]++;
            perCorner[cOverlapPairs[i][1]]++;
        }
    }
    for (u32 const n : perCorner) {
        if (n >= 2)
            return Overlap::eCornerInSeveral;
    }
    return any ? Overlap::ePairwise : Overlap::eNone;
}

bool check(bool condition, std::string_view caseName, std::string_view what)
{
    if (!condition)
        std::println(stderr, "  '{}': {}", caseName, what);
    return condition;
}

[[nodiscard]] fext relError(fext value, fext reference)
{
    return std::abs(value - reference) / reference;
}

bool cornerCutClosedForms(auto areaOf, fext bound)
{
    bool ok { true };
    for (auto const& c : CuratedCases()) {
        if (std::isnan(c.closedForm))
            continue;
        fext const area { areaOf(c.dop) };
        fext const err { relError(area, c.closedForm) };
        ok &= check(err <= bound, c.name, std::format("area {} vs closed form {}: rel error {:.3e} > {:.0e}", area, c.closedForm, err, bound));
    }
    return ok;
}

bool cornerCutF64ClosedForms()
{
    return cornerCutClosedForms([](Dop14<f32> const& dop) { return static_cast<fext>(CornerCut(Widen(dop))); }, 1e-15L);
}

bool cornerCutF32ClosedForms()
{
    return cornerCutClosedForms([](Dop14<f32> const& dop) { return static_cast<fext>(CornerCut(dop)); }, 1e-6L);
}

// left to right: (1 + 2^-24) + 2^-24 is 1, 1 + (2^-24 + 2^-24) is not
bool projectSumsLeftToRight()
{
    f32 const  tiny { 0x1p-24f };
    auto const p { Project(Vec3<f32> { .x = 1.f, .y = tiny, .z = tiny }) };
    bool       ok { check(p[3] == 1.f, "n3", std::format("{:a} != 1", p[3])) };
    ok &= check(p[4] == 1.f - tiny, "n4", std::format("{:a} != 1 - 2^-24", p[4]));
    auto const q { Project(Vec3<f32> { .x = 1.f, .y = -tiny, .z = -tiny }) };
    ok &= check(q[6] == 1.f, "n6", std::format("{:a} != 1", q[6]));
    return ok;
}

bool faceClipClosedForms()
{
    bool ok { true };
    for (auto const& c : CuratedCases()) {
        if (std::isnan(c.closedForm))
            continue;
        fext const e64 { relError(FaceClip(Widen(c.dop)), c.closedForm) };
        fext const e32 { relError(FaceClip(c.dop), c.closedForm) };
        ok &= check(e64 <= 1e-15L, c.name, std::format("FaceClip<f64> rel error {:.3e}", e64));
        ok &= check(e32 <= 1e-6L, c.name, std::format("FaceClip<f32> rel error {:.3e}", e32));
    }
    return ok;
}

// both references are exact: only sqrt(3) and the long double round
bool refsClosedForms()
{
    bool ok { true };
    for (auto const& c : CuratedCases()) {
        Dop14<f64> const dop { Widen(c.dop) };
        RefArea const    cdd { RefCdd(dop) };
        RefArea const    brute { RefBrute(dop) };
        ok &= check(cdd.ok && cdd.vertices > 0, c.name, "cdd failed or found the case empty");
        ok &= check(RefsAgree(cdd, brute), c.name, std::format("cdd {} ({} vertices) vs brute force {} ({})", cdd.area, cdd.vertices, brute.area, brute.vertices));
        if (std::isnan(c.closedForm))
            continue;
        fext const err { relError(brute.area, c.closedForm) };
        ok &= check(err <= 1e-18L, c.name, std::format("area {} vs closed form {}: rel error {:.3e}", brute.area, c.closedForm, err));
    }
    return ok;
}

[[nodiscard]] Dop14<f32> box(Vec3<f32> lo, Vec3<f32> hi)
{
    std::vector<Vec3<f32>> corners;
    for (u32 k { 0 }; k < 8; k++)
        corners.push_back({ .x = (k & 1) ? hi.x : lo.x, .y = (k & 2) ? hi.y : lo.y, .z = (k & 4) ? hi.z : lo.z });
    return Fit(corners);
}

// position and scale cost the references nothing
bool refsTranslatedAndScaled()
{
    struct Case {
        std::string_view name;
        Dop14<f32>       dop;
        fext             area;
    };
    Case const cases[] {
        { .name = "cube at 1e5", .dop = box({ .x = 1e5f, .y = 1e5f, .z = 1e5f }, { .x = 1e5f + 1.f, .y = 1e5f + 1.f, .z = 1e5f + 1.f }), .area = 6.L },
        { .name = "cube at -3e4", .dop = box({ .x = -3e4f, .y = 0.f, .z = 2e4f }, { .x = -3e4f + 1.f, .y = 1.f, .z = 2e4f + 1.f }), .area = 6.L },
        { .name = "cube of 2^-60", .dop = box({ .x = 0.f, .y = 0.f, .z = 0.f }, { .x = 0x1p-60f, .y = 0x1p-60f, .z = 0x1p-60f }), .area = 6.L * 0x1p-120L },
    };
    bool ok { true };
    for (auto const& c : cases) {
        RefArea const cdd { RefCdd(Widen(c.dop)) };
        RefArea const brute { RefBrute(Widen(c.dop)) };
        ok &= check(cdd.area == c.area && brute.area == c.area, c.name, std::format("cdd {} and brute force {} vs {}", cdd.area, brute.area, c.area));
    }
    return ok;
}

// a point is nonempty as stored only if its diagonal sums are exact; an
// axis-aligned segment, since a general segment's DOP is a polygon or solid
bool refsDegenerate()
{
    struct Case {
        std::string_view name;
        Dop14<f32>       dop;
        u32              vertices;
    };
    f32 const                    tiny { 0x1p-24f };
    std::vector<Vec3<f32>> const segment { { .x = 0.f, .y = 1.f, .z = 2.f }, { .x = 3.f, .y = 1.f, .z = 2.f } };
    Dop14<f32>                   crossed { box({ .x = 0.f, .y = 0.f, .z = 0.f }, { .x = 1.f, .y = 1.f, .z = 1.f }) };
    crossed[6] = 2.5f;
    crossed[7] = 2.f;
    Case const cases[] {
        { .name = "exact point", .dop = Init(Vec3<f32> { .x = 1.f, .y = 2.f, .z = 4.f }), .vertices = 1 },
        { .name = "rounded point", .dop = Init(Vec3<f32> { .x = 1.f, .y = tiny, .z = tiny }), .vertices = 0 },
        { .name = "axis segment", .dop = Fit(segment), .vertices = 2 },
        { .name = "min above max", .dop = crossed, .vertices = 0 },
    };
    bool ok { true };
    for (auto const& c : cases) {
        for (RefArea const r : { RefCdd(Widen(c.dop)), RefBrute(Widen(c.dop)) }) {
            ok &= check(r.ok && r.vertices == c.vertices && r.area == 0.L, c.name,
                std::format("ok {} vertices {} area {}, expected {} vertices and no area", r.ok, r.vertices, r.area, c.vertices));
        }
    }
    return ok;
}

// 1000 cases per generator, both fits
bool refsAgreeOnGenerated()
{
    bool ok { true };
    u32  dops { 0 };
    u32  zeroArea { 0 };
    for (u32 g { 0 }; g < strGen.size(); g++) {
        for (u64 i { 0 }; i < 1000; i++) {
            Fits const fits { Generate(static_cast<Gen>(g), 1, i) };
            for (Dop14<f64> const& dop : { Widen(fits.fit32), fits.fit64 }) {
                RefArea const cdd { RefCdd(dop) };
                RefArea const brute { RefBrute(dop) };
                dops++;
                zeroArea += cdd.area > 0.L ? 0 : 1;
                ok &= check(RefsAgree(cdd, brute), strGen[g],
                    std::format("case {}: cdd {} ({} vertices, ok {}) vs brute force {} ({})", i, cdd.area, cdd.vertices, cdd.ok, brute.area, brute.vertices));
            }
        }
    }
    std::println("  {} DOPs, {} with a zero area", dops, zeroArea);
    return ok;
}

// integer points give tight slabs, where the formula is exact, overlaps
// included
bool cornerCutF64ExactOnTightDops()
{
    u64  state { 12345 };
    auto next { [&state](i32 lo, i32 hi) {
        state = state * 6364136223846793005ull + 1442695040888963407ull;
        return lo + static_cast<i32>((state >> 33) % static_cast<u64>(hi - lo + 1));
    } };
    auto point { [&next](i32 e) { return Vec3<f64> { .x = static_cast<f64>(next(-e, e)), .y = static_cast<f64>(next(-e, e)), .z = static_cast<f64>(next(-e, e)) }; } };

    bool ok { true };
    u32  several { 0 };
    u32  evaluated { 0 };
    fext maxRel { 0.L };
    for (u32 c { 0 }; c < 2000; c++) {
        i32 const  extent { next(1, 40) };
        Dop14<f64> dop {};
        if (c % 2 == 0) {
            std::vector<Vec3<f64>> points;
            for (i32 n { next(4, 15) }; n > 0; n--)
                points.push_back(point(extent));
            dop = Fit(points);
        } else {
            dop = Init(point(extent), point(extent), point(extent));
        }
        RefArea const r { RefCdd(dop) };
        if (!(r.area > 0.L))
            continue;
        AreaTerms<f64> terms {};
        f64 const      area { CornerCut(dop, &terms) };
        fext const     err { relError(area, r.area) };
        evaluated++;
        several += classify(terms) == Overlap::eCornerInSeveral ? 1 : 0;
        maxRel = std::max(maxRel, err);
        ok &= check(err <= 1e-14L, "tight DOP", std::format("case {}: CornerCut<f64> {} vs cdd {}: rel error {:.3e}", c, area, r.area, err));
    }
    ok &= check(several * 2 > evaluated, "tight DOP", std::format("only {} of {} cases with a corner in two or more overlaps", several, evaluated));
    std::println("  {} cases, {} with a corner in two or more overlaps, max rel error {:.2Le}", evaluated, several, maxRel);
    return ok;
}

// the fit and CornerCut round each diagonal sum in the same order, and rounding
// is monotonic: a corner-cut depth cannot come out negative
bool cornerCutDepthsNonnegative()
{
    bool ok { true };
    for (u32 g { 0 }; g < strGen.size(); g++) {
        for (u64 i { 0 }; i < 20000; i++) {
            Fits const     fits { Generate(static_cast<Gen>(g), 3, i) };
            AreaTerms<f32> t32 {};
            AreaTerms<f64> t64 {};
            static_cast<void>(CornerCut(fits.fit32, &t32));
            static_cast<void>(CornerCut(fits.fit64, &t64));
            for (u32 k { 0 }; k < 8; k++) {
                ok &= check(t32.d[k] >= 0.f, strGen[g], std::format("case {}: f32 depth {:a}", i, t32.d[k]));
                ok &= check(t64.d[k] >= 0., strGen[g], std::format("case {}: f64 depth {:a}", i, t64.d[k]));
            }
        }
    }
    return ok;
}

struct Test {
    std::string_view name;
    bool (*run)();
};

constexpr Test cTests[] {
    { .name = "cornercut_f64_closed_forms", .run = &cornerCutF64ClosedForms },
    { .name = "cornercut_f32_closed_forms", .run = &cornerCutF32ClosedForms },
    { .name = "project_sums_left_to_right", .run = &projectSumsLeftToRight },
    { .name = "faceclip_closed_forms", .run = &faceClipClosedForms },
    { .name = "refs_closed_forms", .run = &refsClosedForms },
    { .name = "refs_translated_and_scaled", .run = &refsTranslatedAndScaled },
    { .name = "refs_degenerate", .run = &refsDegenerate },
    { .name = "refs_agree_on_generated", .run = &refsAgreeOnGenerated },
    { .name = "cornercut_f64_exact_on_tight_dops", .run = &cornerCutF64ExactOnTightDops },
    { .name = "cornercut_depths_nonnegative", .run = &cornerCutDepthsNonnegative },
};

}

// dop14_tests [<name>...]; no argument runs all
int main(int argc, char** argv)
{
    bool ok { true };
    u32  ran { 0 };
    for (auto const& t : cTests) {
        bool selected { argc == 1 };
        for (i32 i { 1 }; i < argc; i++)
            selected |= std::string_view { argv[i] } == t.name;
        if (!selected)
            continue;
        ran++;
        bool const passed { t.run() };
        std::println("{} {}", passed ? "pass" : "FAIL", t.name);
        ok &= passed;
    }
    if (ran == 0) {
        std::println(stderr, "no test matches");
        return EXIT_FAILURE;
    }
    return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
