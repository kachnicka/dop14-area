#include "cases.h"
#include "corner_cutting.h"
#include "face_clipping.h"
#include "refs.h"
#include "timing.h"

#include <algorithm>
#include <charconv>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <optional>
#include <print>
#include <span>
#include <string>
#include <string_view>
#include <type_traits>
#include <vector>

#include <sys/wait.h>
#include <unistd.h>

namespace {

constexpr i32 cExitOk { 0 };
constexpr i32 cExitFailed { 1 };
constexpr i32 cExitUsage { 2 };

void usage()
{
    std::println(stderr, "usage: dop14-area [--count N] [--seed S] [--jobs N]");
}

[[nodiscard]] fext relError(fext value, fext reference)
{
    return std::abs(value - reference) / reference;
}

[[nodiscard]] fext centreDistance(Dop14<f64> const& dop)
{
    fext const cx { (static_cast<fext>(dop[0]) + dop[1]) / 2.L };
    fext const cy { (static_cast<fext>(dop[2]) + dop[3]) / 2.L };
    fext const cz { (static_cast<fext>(dop[4]) + dop[5]) / 2.L };
    return std::sqrt(cx * cx + cy * cy + cz * cz);
}

// the measured methods, in report order; each on the fit of its own type
constexpr u32                                    cMethods { 4 };
inline constexpr std::array<std::string_view, 4> strMethod { "CornerCut<f32>", "CornerCut<f64>", "FaceClip<f32>", "FaceClip<f64>" };

using Errors = std::array<std::vector<f64>, cMethods>;

struct GenStats {
    u64 cases { 0 };
    // DOPs, two per case, where cdd fails or the references disagree
    u64 refsFailed { 0 };
    // cases left out: a fit with an exact area of zero
    u64 zeroArea { 0 };
    u64 nonFinite { 0 };
    // relative errors against the exact area of the same fit
    Errors errors;
    // far only: by the box centre's distance from the origin, below 1e3,
    // below 1e4, beyond
    std::array<Errors, 3> byDistance;
};

[[nodiscard]] fext exactArea(GenStats& s, Dop14<f64> const& dop)
{
    RefArea const cdd { RefCdd(dop) };
    s.refsFailed += RefsAgree(cdd, RefBrute(dop)) ? 0 : 1;
    return cdd.area;
}

void add(GenStats& s, Fits const& fits, bool far)
{
    s.cases++;
    fext const exact32 { exactArea(s, Widen(fits.fit32)) };
    fext const exact64 { exactArea(s, fits.fit64) };
    if (!(exact32 > 0.L && exact64 > 0.L)) {
        s.zeroArea++;
        return;
    }

    fext const areas[cMethods] { CornerCut(fits.fit32), CornerCut(fits.fit64), FaceClip(fits.fit32), FaceClip(fits.fit64) };
    fext const exact[cMethods] { exact32, exact64, exact32, exact64 };
    for (fext const a : areas) {
        if (!std::isfinite(a)) {
            s.nonFinite++;
            return;
        }
    }
    fext const distance { centreDistance(fits.fit64) };
    for (u32 m { 0 }; m < cMethods; m++) {
        f64 const err { static_cast<f64>(relError(areas[m], exact[m])) };
        s.errors[m].push_back(err);
        if (far)
            s.byDistance[distance < 1e3L ? 0 : (distance < 1e4L ? 1 : 2)][m].push_back(err);
    }
}

// part follows into's cases: appending keeps case order, so the result equals
// one pass
void merge(GenStats& into, GenStats const& part)
{
    into.cases += part.cases;
    into.refsFailed += part.refsFailed;
    into.zeroArea += part.zeroArea;
    into.nonFinite += part.nonFinite;
    auto const append { [](Errors& to, Errors const& from) {
        for (u32 m { 0 }; m < cMethods; m++)
            to[m].insert(to[m].end(), from[m].begin(), from[m].end());
    } };
    append(into.errors, part.errors);
    for (u32 d { 0 }; d < 3; d++)
        append(into.byDistance[d], part.byDistance[d]);
}

class Writer {
public:
    std::vector<u8> bytes;

    template<typename T>
    void put(T const& v)
    {
        static_assert(std::is_trivially_copyable_v<T>);
        auto const* p { reinterpret_cast<u8 const*>(&v) };
        bytes.insert(bytes.end(), p, p + sizeof(T));
    }
    template<typename T>
    void putVector(std::vector<T> const& v)
    {
        put(static_cast<u64>(v.size()));
        auto const* p { reinterpret_cast<u8 const*>(v.data()) };
        bytes.insert(bytes.end(), p, p + v.size() * sizeof(T));
    }
};

class Reader {
    std::span<u8 const> bytes;
    u64                 at { 0 };

public:
    bool ok { true };

    explicit Reader(std::span<u8 const> b)
        : bytes(b)
    {
    }

    template<typename T>
    void get(T& v)
    {
        static_assert(std::is_trivially_copyable_v<T>);
        ok = ok && at + sizeof(T) <= bytes.size();
        if (!ok)
            return;
        std::copy_n(bytes.data() + at, sizeof(T), reinterpret_cast<u8*>(&v));
        at += sizeof(T);
    }
    template<typename T>
    void getVector(std::vector<T>& v)
    {
        u64 n { 0 };
        get(n);
        ok = ok && n <= (bytes.size() - at) / sizeof(T);
        if (!ok)
            return;
        v.resize(n);
        std::copy_n(bytes.data() + at, n * sizeof(T), reinterpret_cast<u8*>(v.data()));
        at += n * sizeof(T);
    }
    [[nodiscard]] bool done() const { return ok && at == bytes.size(); }
};

// the fields in one order for both directions
template<typename Io, typename Stats>
void transfer(Io& io, Stats& s, auto value, auto vector)
{
    value(io, s.cases);
    value(io, s.refsFailed);
    value(io, s.zeroArea);
    value(io, s.nonFinite);
    for (u32 m { 0 }; m < cMethods; m++) {
        vector(io, s.errors[m]);
        for (u32 d { 0 }; d < 3; d++)
            vector(io, s.byDistance[d][m]);
    }
}

void serialize(Writer& w, GenStats const& s)
{
    transfer(w, s, [](Writer& io, auto const& v) { io.put(v); }, [](Writer& io, auto const& v) { io.putVector(v); });
}

void deserialize(Reader& r, GenStats& s)
{
    transfer(r, s, [](Reader& io, auto& v) { io.get(v); }, [](Reader& io, auto& v) { io.getVector(v); });
}

// "max / p99", or "-" without cases
[[nodiscard]] std::string maxAndP99(std::vector<f64> values)
{
    if (values.empty())
        return "-";
    f64 const worst { *std::max_element(values.begin(), values.end()) };
    u64 const rank { (values.size() * 99 + 99) / 100 - 1 };
    std::nth_element(values.begin(), values.begin() + static_cast<i64>(rank), values.end());
    return std::format("{:.1e} / {:.1e}", worst, values[rank]);
}

void printRow(std::string_view label, Errors const& errors)
{
    std::string line { std::format("{:<16} {:>8}", label, errors[0].size()) };
    for (std::vector<f64> const& e : errors)
        line += std::format(" {:>19}", maxAndP99(e));
    std::println("{}", line);
}

template<typename T>
[[nodiscard]] std::optional<T> parseNumber(std::string_view s)
{
    T value {};
    auto const [end, ec] { std::from_chars(s.data(), s.data() + s.size(), value) };
    if (ec != std::errc {} || end != s.data() + s.size())
        return std::nullopt;
    return value;
}

struct Options {
    u64 count { 10'000 };
    u64 seed { 1 };
    u32 jobs { static_cast<u32>(std::max(1L, sysconf(_SC_NPROCESSORS_ONLN))) };
};

[[nodiscard]] std::optional<Options> parseOptions(std::span<char* const> args)
{
    Options options;
    bool    valid { args.size() % 2 == 0 };
    for (u64 i { 0 }; valid && i < args.size(); i += 2) {
        std::string_view const flag { args[i] };
        std::string_view const value { args[i + 1] };
        if (flag == "--count") {
            auto const v { parseNumber<u64>(value) };
            valid         = v && *v > 0;
            options.count = v.value_or(options.count);
        } else if (flag == "--seed") {
            auto const v { parseNumber<u64>(value) };
            valid        = v.has_value();
            options.seed = v.value_or(options.seed);
        } else if (flag == "--jobs") {
            auto const v { parseNumber<u32>(value) };
            valid        = v && *v > 0;
            options.jobs = v.value_or(options.jobs);
        } else {
            valid = false;
        }
    }
    if (!valid)
        return std::nullopt;
    return options;
}

// the same contiguous slice of every generator's cases
[[nodiscard]] std::vector<GenStats> evaluateSlice(Options const& options, u32 job)
{
    u64 const             begin { options.count * job / options.jobs };
    u64 const             end { options.count * (job + 1) / options.jobs };
    std::vector<GenStats> stats(strGen.size());
    for (u32 g { 0 }; g < strGen.size(); g++) {
        Gen const gen { static_cast<Gen>(g) };
        for (u64 i { begin }; i < end; i++)
            add(stats[g], Generate(gen, options.seed, i), gen == Gen::eFar);
    }
    return stats;
}

[[nodiscard]] bool writeAll(i32 fd, std::span<u8 const> bytes)
{
    while (!bytes.empty()) {
        ssize_t const n { write(fd, bytes.data(), bytes.size()) };
        if (n <= 0)
            return false;
        bytes = bytes.subspan(static_cast<u64>(n));
    }
    return true;
}

// cddlib has global state, so workers are processes; each pipes its
// statistics back
[[nodiscard]] std::optional<std::vector<GenStats>> evaluateForked(Options const& options)
{
    struct Worker {
        pid_t pid;
        i32   fd;
    };
    std::fflush(stdout);
    std::fflush(stderr);
    std::vector<Worker> workers;
    for (u32 job { 0 }; job < options.jobs; job++) {
        i32 fds[2];
        if (pipe(fds) != 0) {
            std::println(stderr, "dop14-area: pipe for worker {} failed", job);
            return std::nullopt;
        }
        pid_t const pid { fork() };
        if (pid < 0) {
            std::println(stderr, "dop14-area: fork of worker {} failed", job);
            return std::nullopt;
        }
        if (pid == 0) {
            close(fds[0]);
            for (Worker const& w : workers)
                close(w.fd);
            Writer writer;
            for (GenStats const& s : evaluateSlice(options, job))
                serialize(writer, s);
            _exit(writeAll(fds[1], writer.bytes) ? EXIT_SUCCESS : EXIT_FAILURE);
        }
        close(fds[1]);
        workers.push_back({ .pid = pid, .fd = fds[0] });
    }

    std::vector<GenStats> stats(strGen.size());
    bool                  ok { true };
    std::vector<u8>       bytes;
    for (u32 job { 0 }; job < options.jobs; job++) {
        bytes.clear();
        u8 chunk[1 << 16];
        for (ssize_t n; (n = read(workers[job].fd, chunk, sizeof(chunk))) > 0;)
            bytes.insert(bytes.end(), chunk, chunk + n);
        close(workers[job].fd);
        i32 status { 0 };
        waitpid(workers[job].pid, &status, 0);

        std::vector<GenStats> part(strGen.size());
        Reader                reader { bytes };
        for (GenStats& s : part)
            deserialize(reader, s);
        if (!WIFEXITED(status) || WEXITSTATUS(status) != EXIT_SUCCESS || !reader.done()) {
            std::println(stderr, "dop14-area: worker {} failed (status {:#x}, {} bytes)", job, status, bytes.size());
            ok = false;
            continue;
        }
        for (u64 g { 0 }; g < stats.size(); g++)
            merge(stats[g], part[g]);
    }
    if (!ok)
        return std::nullopt;
    return stats;
}

// false when the references disagree
[[nodiscard]] bool printCorrectness(std::vector<GenStats> const& stats, Options const& options, f64 seconds)
{
    u64 refsFailed { 0 };
    u64 zeroArea { 0 };
    u64 nonFinite { 0 };
    for (GenStats const& s : stats) {
        refsFailed += s.refsFailed;
        zeroArea += s.zeroArea;
        nonFinite += s.nonFinite;
    }
    u64 const dops { 2 * options.count * strGen.size() };
    std::println("correctness: seed {}, {} cases per generator, {} job(s), {:.1f} s\n", options.seed, options.count, options.jobs, seconds);
    if (refsFailed == 0)
        std::println("references: cdd and brute force agree exactly on all {} DOPs (both fits of every case)\n", dops);
    else
        std::println("references: disagree or cdd failed on {} of {} DOPs: a bug in this tool\n", refsFailed, dops);

    std::println("relative error against the exact area of the same fit, max / p99; <T> fits and evaluates in T");
    std::println("{:<16} {:>8} {:>19} {:>19} {:>19} {:>19}", "generator", "cases", strMethod[0], strMethod[1], strMethod[2], strMethod[3]);
    for (u32 g { 0 }; g < stats.size(); g++)
        printRow(strGen[g], stats[g].errors);
    constexpr std::array<std::string_view, 3> cDistances { "far, 1e2 to 1e3", "far, 1e3 to 1e4", "far, 1e4 to 1e5" };
    for (u32 d { 0 }; d < 3; d++)
        printRow(cDistances[d], stats[static_cast<u32>(Gen::eFar)].byDistance[d]);
    if (zeroArea > 0)
        std::println("cases with a zero exact area, left out above: {}", zeroArea);
    if (nonFinite > 0)
        std::println("cases with a non-finite result, left out above: {}", nonFinite);
    return refsFailed == 0;
}

}

int main(int argc, char** argv)
{
    auto const options { parseOptions({ argv + 1, static_cast<u64>(argc - 1) }) };
    if (!options) {
        usage();
        return cExitUsage;
    }

    auto const            start { std::chrono::steady_clock::now() };
    std::vector<GenStats> stats;
    if (options->jobs == 1) {
        stats = evaluateSlice(*options, 0);
    } else {
        auto forked { evaluateForked(*options) };
        if (!forked)
            return cExitFailed;
        stats = std::move(*forked);
    }
    std::chrono::duration<f64> const elapsed { std::chrono::steady_clock::now() - start };
    bool const                       ok { printCorrectness(stats, *options, elapsed.count()) };

    PrintTiming(options->seed);
    return ok ? cExitOk : cExitFailed;
}
