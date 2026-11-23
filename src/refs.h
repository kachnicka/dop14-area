#pragma once

// Two exact areas of the stored slabs, read as exact half-spaces: rational
// throughout, so only sqrt(3) and the final long double round. A flat DOP
// counts both sides, as CornerCut's box term does; an empty one has area 0.
// f32 DOPs widen exactly.

#include "dop14.h"

struct RefArea {
    fext area;
    // of the stored polytope; 0 when it is empty
    u32 vertices;
    // false when cdd fails
    bool ok;
};

// cddlib in GMP rational mode: vertices and incidence by the double
// description method. One thread only: cddlib has global state.
[[nodiscard]] RefArea RefCdd(Dop14<f64> const& dop);

// every triple of planes from distinct slabs, solved in GMP integers; no cdd
[[nodiscard]] RefArea RefBrute(Dop14<f64> const& dop);

// the same rational areas give the same long double; a mismatch is a bug in
// this tool
[[nodiscard]] inline bool RefsAgree(RefArea const& cdd, RefArea const& brute)
{
    return cdd.ok && brute.ok && cdd.vertices == brute.vertices && cdd.area == brute.area;
}
