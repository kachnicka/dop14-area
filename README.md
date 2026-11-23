# dop14-area

Correctness verification and performance benchmark for the 14-DOP surface area algorithms from
M. Káčerik and J. Bittner, [SAH-Optimized k-DOP Hierarchies for Ray Tracing](https://doi.org/10.1145/3675391),
Proc. ACM Comput. Graph. Interact. Tech. 7(3), 2024: face clipping in section 3.2 and corner cutting in section 3.3.

## Erratum in the paper's listing

Corner cutting's Eq. 3 gives the overlap leg s = (√3(d₁ + d₂) − u)/2, with d measured along unit normals.
The listing measures depths along the raw diagonals (√3 times longer) and computes s = d₁ + d₂ − u, twice Eq. 3's value.
Eq. 4's correction (2 − √3)s² is therefore (2 − √3)/4 · s² ≈ 0.0670 s² in the listing's terms;
the listing in the paper's Appendix A multiplies by (2 − √3)/2 ≈ 0.1340 instead.
This repository uses (2 − √3)/4.

## Build

C++23, GCC or Clang, CMake 3.25 or newer, Ninja, and cddlib's GMP build
with GMP from the system packages:

- Arch: `pacman -S cddlib gmp`
- Debian/Ubuntu: `apt install libcdd-dev libgmp-dev` (not tested here)

CMake finds cddlib through the pkg-config module `cddgmp`, or else through
`find_library(cddgmp)`; the headers may sit in `cddlib/` (upstream) or `cdd/`.

```
cmake -B build -G Ninja
cmake --build build
ctest --test-dir build
```

## Usage

```
dop14-area [--count N] [--seed S] [--jobs N]
```

One run checks correctness, then measures speed, and prints both as plain-text tables.
It evaluates `--count` cases (default 10000) of each data [generator](#data) under `--seed` (default 1).
The exact references cost about 2 ms per case, so the correctness part forks `--jobs` worker processes
(default: one per CPU; cddlib keeps global state, so there are no threads).
The report does not depend on `--jobs`.

## Code structure

- `src/dop14.h`: the DOP type, its fit in f32 or f64 (diagonal sums left to right) and the merge.
- `src/face_clipping.h`: `FaceClip`, face clipping templated on the float type.
- `src/corner_cutting.h`: `CornerCut`, corner cutting templated on the float type.
- `src/ref_cdd.cpp`: the first exact reference.
  cddlib in GMP rational mode gives the exact vertices and vertex-facet incidence;
  each facet's polygon is ordered and summed in rationals.
- `src/ref_brute.cpp`: the second exact reference, without cdd: every triple of planes from distinct slabs
  is solved in GMP integers, the solutions inside all slabs are the vertices,
  and each face's vertices are ordered and summed exactly.
- `src/cases.cpp`: the curated cases with closed-form areas, which the tests check, and the seeded generators.

In both references, only the final √3 of the diagonal faces and the conversion to long double round (about 1e-19).
They take the stored slab planes as exact half-spaces.
A flat DOP counts both sides (twice its polygon), as the box term of `CornerCut` does;
lower-dimensional DOPs need no special handling, and slabs that enclose nothing give area 0.

The references check each other on every DOP: the same vertex count and the same area, bit for bit,
since both reduce the area to the same two rationals before rounding. A mismatch is a bug in this tool.

## Data

| Generator | Cases |
|---|---|
| `cloud` | 4 to 15 points, per-axis extent 0.01 to 10 |
| `sliver` | a cloud with one axis 1e-5 to 1e-2 |
| `triangle` | 3 points, fitted as one triangle |
| `segment` | 2 points |
| `flat` | 3 to 10 points on one of the seven DOP planes |
| `box` | the 8 corners of a randomly oriented box (deep diagonal cuts) |
| `merged` | the merge of 2 to 7 fitted triangles |
| `far` | any of the above, translated by 1e2 to 1e5 |
| `scale` | any of the above, scaled by 1e-4 to 1e6 |

Each case is a point set generated in f64 and stored as f32, as a vertex buffer holds it,
then fitted twice: in f32 and in f64. The f32 points widen exactly, so both DOPs bound the same points.
Extents are log-uniform, box centres within [−1, 1]³. A case depends only on
(seed, generator, index), and on the platform's `pow`, `sin` and `cos`.

Exit codes:
- 0 when the references agree on every DOP
- 1 when they disagree, or cdd or a worker process fails
- 2 for a bad command line

## Results

`dop14-area --count 100000` on an Intel Core i9-10900X, GCC 16.2, Release build, cddlib 0.94n, GMP 6.3:
900k cases, and the references agree exactly on all 1.8M DOPs.

Relative error against the exact area of the same fit (`<T>`: fitted and evaluated in T), max / p99.

| Generator | Cases | `CornerCut<f32>` | `CornerCut<f64>` | `FaceClip<f32>` | `FaceClip<f64>` |
|---|---:|---:|---:|---:|---:|
| cloud | 100000 | 7.5e-6 / 8.0e-7 | 1.3e-15 / 2.2e-16 | 8.1e-6 / 9.2e-7 | 6.1e-16 / 2.3e-16 |
| sliver | 100000 | 3.5e-5 / 2.0e-6 | 1.6e-15 / 1.9e-16 | 3.9e-5 / 2.1e-6 | 6.1e-16 / 2.2e-16 |
| triangle | 100000 | 3.6e-5 / 2.2e-6 | 5.8e-15 / 3.8e-16 | 5.3e-5 / 1.9e-6 | 1.8e-15 / 2.0e-16 |
| segment | 100000 | 4.0e-4 / 8.6e-6 | 4.5e-14 / 1.3e-15 | 3.3e-4 / 6.5e-6 | 4.2e-14 / 1.9e-16 |
| flat | 100000 | 1.1e-3 / 2.5e-5 | 2.8e-12 / 3.3e-14 | 7.6e-4 / 7.6e-6 | 1.6e-13 / 4.3e-15 |
| box | 100000 | 3.7e-6 / 6.6e-7 | 9.2e-15 / 5.4e-16 | 3.6e-6 / 6.0e-7 | 1.8e-15 / 2.8e-16 |
| merged | 100000 | 1.7e-5 / 1.1e-6 | 1.6e-15 / 2.5e-16 | 2.3e-5 / 1.1e-6 | 6.7e-16 / 2.3e-16 |
| far | 99957 | 12 / 6.4e-2 | 8.9e-13 / 5.5e-15 | 9.7 / 5.2e-2 | 2.7e-15 / 1.4e-16 |
| scale | 100000 | 3.9e-4 / 6.7e-6 | 2.6e-13 / 5.9e-15 | 8.4e-5 / 3.6e-6 | 3.2e-14 / 4.8e-16 |
| far, 1e2 to 1e3 | 33073 | 7.3e-2 / 1.6e-3 | 8.9e-13 / 6.1e-15 | 5.6e-2 / 1.2e-3 | 1.7e-15 / 1.4e-16 |
| far, 1e3 to 1e4 | 33424 | 0.24 / 1.6e-2 | 2.1e-13 / 5.3e-15 | 0.36 / 1.2e-2 | 2.7e-15 / 1.4e-16 |
| far, 1e4 to 1e5 | 33460 | 12 / 0.15 | 1.0e-13 / 5.4e-15 | 9.7 / 0.13 | 2.2e-16 / 1.4e-16 |

43 `far` cases have an exact area of zero and are left out: in 15, rounding
to f32 far from the origin snaps a small shape's points onto a line, and both
fits are segments; in the other 28, only the f32 fit degenerates, its rounded
diagonal sums leaving a segment (22), a point (1) or nothing (5).

Time per DOP, one thread, one DOP per loop iteration, on the f32 fits:

| Per DOP | triangle | merged | cloud | box |
|---|---:|---:|---:|---:|
| corner cutting: `CornerCut<f32>` | 15.6 ns | 15.5 ns | 15.6 ns | 15.6 ns |
| face clipping: `FaceClip<f32>` | 5 µs | 5 µs | 5 µs | 5 µs |
| cdd float: vertex enumeration | 44 µs | 51 µs | 50 µs | 57 µs |
| cdd float: enumeration and area | 58 µs | 70 µs | 70 µs | 80 µs |
| cdd rational: vertex enumeration | 623 µs | 706 µs | 701 µs | 783 µs |
| cdd rational: `RefCdd` | 811 µs | 940 µs | 919 µs | 1039 µs |
| brute force: `RefBrute` | 51 µs | 57 µs | 57 µs | 61 µs |
| cdd float enumeration / `CornerCut<f32>` | 2811× | 3330× | 3215× | 3665× |
| cdd float enumeration / `FaceClip<f32>` | 9× | 11× | 11× | 13× |

What the results say:

- Corner cutting is exact on tight slabs, edge-overlap corrections included, as the paper states:
  on DOPs of integer points, where every slab is exact,
  `CornerCut<f64>` matches cdd within 1e-15 although nearly all of them
  have a corner in two or more overlaps (the `cornercut_f64_exact_on_tight_dops` test).
  The same holds for the closed forms.
- Fitted and evaluated in f64, both methods are exact up to rounding, far from the origin included:
  `CornerCut<f64>` keeps its p99 at or below 3.3e-14 and its worst case at 2.8e-12,
  on flat DOPs, where the formula's terms cancel;
  `FaceClip<f64>` stays within 1.6e-13.
- In f32, the worst case near the origin ranges from about 4e-6 on boxes to 1.1e-3 on flat DOPs,
  and the error grows about tenfold per decade of distance from the origin
  (p99 1.6e-3, 1.6e-2 and 0.15 in the three far bands);
  beyond 1e4 the result can be off by more than the area itself.
  Both compute with coordinates at the DOP's position, where the f32 spacing is coarse against a small DOP;
  corner cutting also assumes tight slabs, which the f32 fit's rounded diagonal sums do not quite give.
- Corner cutting takes 16 ns, 2800 to 3700 times less than floating-point
  cdd's vertex enumeration and about 300 times less than face clipping;
  face clipping is 9 to 13 times faster than that enumeration, the paper's "about an order of magnitude".
- `RefBrute` is exact without cdd and 16 to 17 times faster than `RefCdd`.

## License

All code in this repository is MIT, see [LICENSE](LICENSE). The repository
contains no third-party source. A built binary links cddlib (GPL-2.0-or-later)
and GMP (LGPL-3.0-or-later or GPL-2.0-or-later), so the binary as a whole is
GPL-2.0-or-later; no source file changes license.
