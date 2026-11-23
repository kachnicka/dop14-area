#pragma once

#include "dop14.h"

#include <array>
#include <string_view>
#include <vector>

// fitted in f32; the closed forms belong to integer points, whose slabs are exact
struct CuratedCase {
    std::string_view name;
    Dop14<f32>       dop;
    // NaN: no closed form
    fext closedForm;
};

[[nodiscard]] std::vector<CuratedCase> CuratedCases();

enum class Gen : u8 {
    eCloud,
    eSliver,
    eTriangle,
    eSegment,
    eFlat,
    eBox,
    eMerged,
    eFar,
    eScale,
};

inline constexpr std::array<std::string_view, 9> strGen {
    "cloud", "sliver", "triangle", "segment", "flat", "box", "merged", "far", "scale"
};

// one point set, stored as f32, fitted in each precision
struct Fits {
    Dop14<f32> fit32;
    Dop14<f64> fit64;
};

// Deterministic per (gen, seed, index), so any range of cases stands alone.
[[nodiscard]] Fits Generate(Gen gen, u64 seed, u64 index);
